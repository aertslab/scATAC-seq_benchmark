import platform
import os
import sys
from ConfigParser import RawConfigParser
from base64 import b64decode
from datetime import datetime
from email.mime.text import MIMEText
from smtplib import SMTP

from MySQLdb import connect


JOBS_IN_QUEUE_STATEMENT = r"""
  SELECT count(*)
  FROM jobQueue
  WHERE jobStatusCode = 1;
"""

WAITING_TIME_STATEMENT = r"""
  SELECT max(timestampdiff( MINUTE, jobRequestTime, now() ))
  FROM jobQueue
  WHERE jobStatusCode = 1;
"""


def create_db_connection(cfg):
    connection = connect(host=cfg.get('dbserver', 'servername'),
                         port=cfg.getint('dbserver', 'port'),
                         user=cfg.get('dbserver', 'username'),
                         passwd=b64decode(cfg.get('dbserver', 'password')),
                         db=cfg.get('dbserver', 'database'))

    # Automatically reconnect to the MySQL server if the connection goes away.
    connection.ping(True)

    return connection


def send_mail(cfg, message):
    subject = cfg.get('compserver', 'email_subject')
    sender = cfg.get('compserver', 'webserver_email_address')
    recipients = map(lambda s: s.strip(), cfg.get('compserver', 'recipient_email_addresses').split(';'))
    mail_server_host = (cfg.get('compserver', 'mail_server_host').strip()
                        if cfg.has_option('compserver', 'mail_server_host')
                        else None)
    mail_server_port = (int(cfg.get('compserver', 'mail_server_port').strip())
                        if cfg.has_option('compserver', 'mail_server_port')
                        else 25)

    for recipient in recipients:
        msg = MIMEText(message)
        msg['Subject'] = subject
        msg['From'] = sender
        msg['To'] = recipient

        s = SMTP(host=mail_server_host, port=mail_server_port)
        s.sendmail(sender, [msg['To']], msg.as_string())
        s.quit()


def query_jobs_count(connection):
    cursor = connection.cursor()
    cursor.execute(JOBS_IN_QUEUE_STATEMENT)
    count = cursor.fetchone()[0]
    cursor.close()

    return count


def query_waiting_time(connection):
    cursor = connection.cursor()
    cursor.execute(WAITING_TIME_STATEMENT)
    waiting_time = cursor.fetchone()[0]
    cursor.close()

    return waiting_time if waiting_time is not None else 0


def main():
    if len(sys.argv) != 2:
        print >> sys.stderr, '\nUsage: iregulon-check-availability check-availability.ini\n\n' + \
                             'Example configfile:\n-------------------\n\n' + \
                             '[dbserver]\n' + \
                             'servername=localhost\n' + \
                             'port=3306\n' + \
                             'username=iregulondaemon\n' + \
                             'password=<base64_encoded_password>\n' + \
                             'database=iregulon_daemon\n' + \
                             '[compserver]\n' + \
                             'email_subject=Problem with iRegulon daemon on ' + platform.node() + ' ...\n' + \
                             'webserver_email_address=noreply@example.com\n' + \
                             'recipient_email_addresses=yourname@example.com\n' + \
                             'max_waiting_time_in_min=20\n' + \
                             'max_jobs_in_queue=5\n' + \
                             'daemon_pid_file=<path_to_iregulon-daemon.pid>\n'
        sys.exit(1)

    ini_filename = sys.argv[1]

    if not os.path.exists(ini_filename):
        print >> sys.stderr, "Config file '{0:s}' does not exist.".format(ini_filename)
        sys.exit(1)
    if not os.path.isfile(ini_filename):
        print >> sys.stderr, "'{0:s}' is not a config file.".format(ini_filename)
        sys.exit(1)

    cfg = RawConfigParser()
    cfg.read(ini_filename)

    daemon_pid_file = os.path.expandvars(cfg.get('compserver', 'daemon_pid_file'))

    if os.path.exists(daemon_pid_file) and os.path.isfile(daemon_pid_file):
        with open(daemon_pid_file, 'r') as fh:
            daemon_pid = fh.readline().rstrip()

        daemon_pid_cmdline_file = '/proc/' + daemon_pid + '/cmdline'

        if not os.path.exists('/proc/' + daemon_pid):
            send_mail(cfg,
                      "The iRegulon daemon with PID={0:s} is not running anymore on {1:s}.".format(daemon_pid,
                                                                                                   platform.node()))
            print datetime.now()
            print "The iRegulon daemon with PID={0:s} is not running anymore on {1:s}.".format(daemon_pid,
                                                                                               platform.node())
        elif os.path.exists(daemon_pid_cmdline_file):
            with open(daemon_pid_cmdline_file, 'r') as fh:
                cmdline = fh.read()

            # Check the PID is indeed from the iregulon-daemon.
            if 'iregulon-daemon' in cmdline:
                # Check the number of jobs and the waiting time in the queue.
                max_jobs_in_queue = cfg.getint('compserver', 'max_jobs_in_queue')
                max_waiting_time_in_min = cfg.getint('compserver', 'max_waiting_time_in_min')

                connection = create_db_connection(cfg)
                jobs_in_queue = query_jobs_count(connection)
                waiting_time = query_waiting_time(connection)
                if jobs_in_queue > max_jobs_in_queue or waiting_time >= max_waiting_time_in_min:
                    send_mail(cfg,
                              "There are {0:d} jobs in the job queue of the iRegulon daemon on {1:s}\nand the waiting time for a requested job exceeds {2:d} minutes.".format(
                                  jobs_in_queue, platform.node(), waiting_time))

                print datetime.now()
                print "# jobs in the queue                      = {0:d}".format(jobs_in_queue)
                print "Highest waiting time for a requested job = {0:d} min".format(waiting_time)

                connection.close()
            else:
                send_mail(cfg,
                          "PID {0:s} is in use by another process.\nDelete the PID file '{1:s}' and start the iRegulon daemon on {2:s} again.".format(
                              daemon_pid, daemon_pid_file, platform.node()))
                print datetime.now()
                print "PID {0:s} is in use by another process.\nDelete the PID file '{1:s}' and start the iRegulon daemon on {2:s} again.".format(
                    daemon_pid, daemon_pid_file, platform.node())

    else:
        print datetime.now()
        print "The iRegulon daemon is not started on {0:s}.\nPID file '{1:s}' could not be found.".format(
            platform.node(), daemon_pid_file)


if __name__ == "__main__":
    main()
