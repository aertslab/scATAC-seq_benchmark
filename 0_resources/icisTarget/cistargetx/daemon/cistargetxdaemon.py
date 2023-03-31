import StringIO
import base64
import datetime
import logging
import os
import paramiko
import shutil
import signal
import socket
import stat
import sys
import tempfile
import time
import traceback
import urllib2
from ConfigParser import RawConfigParser

from job import JobQueue, JobRunningInOtherInstanceError
from jobcontext import JobContext


from ftplib import all_errors as FTP_all_errors
try:
    from ftplib import FTP_TSL as FTP
except ImportError:
    from ftplib import FTP

INI_TEMPLATE = os.path.join(os.path.dirname(__file__), 'cistargetxdaemon-template.ini')


def upload_files(cfg, job_unique_id, results_folder, subdirs_list, filenames_list):
    upload_results = cfg.get('compserver', 'upload_results')

    if upload_results == 'local':
        local_reports_dir = os.path.expandvars(cfg.get('local', 'folder'))

        if not os.path.isdir(local_reports_dir):
            logging.error("Local reports directory {0:s} does not exist.".format(local_reports_dir))
            raise AssertionError("Report results could not be copied.")

        dest_folder = os.path.join(local_reports_dir, str(job_unique_id))
        logging.info("Copy report results to {0:s}.".format(dest_folder))

        shutil.copytree(results_folder, dest_folder)
    elif upload_results == 'ftpserver':
        logging.info("Upload report results via FTP.")
        if cfg.getboolean('ftpserver', 'upload_archive_only'):
            try:
                ftp = FTP()
                ftp.connect(cfg.get('ftpserver', 'servername'),
                            port=(cfg.getint('ftpserver', 'port')
                                  if cfg.has_option('ftpserver', 'port')
                                  else 21))

                ftp.login(cfg.get('ftpserver', 'username'), base64.b64decode(cfg.get('ftpserver', 'password')))

                if cfg.get('ftpserver', 'folder'):
                    ftp.cwd(os.path.expandvars(cfg.get('ftpserver', 'folder')))

                # Upload report zip archive only.
                with open(os.path.join(results_folder, 'archive.zip'), 'rb') as input_fh:
                    ftp.storbinary('STOR {0:s}-archive.zip'.format(str(job_unique_id)), input_fh)

                ftp.quit()
            except FTP_all_errors, e:
                logging.error("FTP error: {0:s}.".format(e))
                raise AssertionError("Something went wrong while uploading the report results.")

            # Extract report zip archive on webserver.
            req = urllib2.Request(url='{0:s}/unzip_report.php'.format(cfg.get('compserver', 'webserver_url')),
                                  data='job_uid={0:s}'.format(str(job_unique_id)))

            try:
                response = urllib2.urlopen(req)
                response.close()
            except urllib2.URLError as e:
                if hasattr(e, 'reason'):
                    logging.error("HTTP error: {0:s}.".format(e.reason))
                    raise AssertionError("Something went wrong while extracting the report results.")
                elif hasattr(e, 'code'):
                    logging.error("HTTP error: {0:s}.".format(e.code))
                    raise AssertionError("Something went wrong while extracting the report results.")
        else:
            try:
                ftp = FTP()
                ftp.connect(cfg.get('ftpserver', 'servername'),
                            port=(cfg.getint('ftpserver', 'port')
                                  if cfg.has_option('ftpserver', 'port')
                                  else 21))

                ftp.login(cfg.get('ftpserver', 'username'), base64.b64decode(cfg.get('ftpserver', 'password')))

                if cfg.get('ftpserver', 'folder'):
                    ftp.cwd(os.path.expandvars(cfg.get('ftpserver', 'folder')))

                ftp.mkd(str(job_unique_id))
                ftp.cwd(str(job_unique_id))

                # First make all subdirectories.
                for subdir in subdirs_list:
                    ftp.mkd(subdir)

                # Upload all files.
                for filename in filenames_list:
                    with open(os.path.join(results_folder, filename), 'rb') as input_fh:
                        ftp.storbinary('STOR {0:s}'.format(filename), input_fh)

                ftp.quit()
            except FTP_all_errors, e:
                logging.error("FTP error: {0:s}.".format(e))
                raise AssertionError("Something went wrong while uploading the report results.")
    elif upload_results == 'sshserver':
        ssh_hostname = cfg.get('sshserver', 'hostname')
        ssh_port = int(cfg.get('sshserver', 'port'))
        ssh_username = cfg.get('sshserver', 'username')
        ssh_private_key_file = os.path.expandvars(cfg.get('sshserver', 'private_key_file'))
        reports_folder = cfg.get('sshserver', 'reports_folder')

        if reports_folder.endswith('/'):
            reports_folder = reports_folder[:-1]

        job_results_zip_file = os.path.join(results_folder, 'archive.zip')
        job_results_folder = os.path.join(reports_folder, job_unique_id)

        logging.info("Upload report results via SSH.")

        # Transfer zip archive to webserver with scp.
        os.system("scp -i {0:s} {1:s} {2:s}@{3:s}:{4:s}-archive.zip".format(ssh_private_key_file,
                                                                            job_results_zip_file,
                                                                            ssh_username,
                                                                            ssh_hostname,
                                                                            job_results_folder))

        try:
            ssh = paramiko.SSHClient()
            ssh.load_system_host_keys()

            ssh_private_key = paramiko.RSAKey.from_private_key_file(filename=ssh_private_key_file)

            ssh.connect(hostname=ssh_hostname,
                        port=ssh_port,
                        username=ssh_username,
                        pkey=ssh_private_key,
                        look_for_keys=False)

            ssh_transport = ssh.get_transport()

            # Send keep alive packets every 10 seconds.
            ssh_transport.set_keepalive(10)

            # Transfer zip archive to webserver.
            #sftp = paramiko.SFTPClient.from_transport(ssh_transport)
            #sftp.put(job_results_zip_file, job_results_folder + '-archive.zip')

            extract_results_zip_script = 'mkdir {0:s}-temp;\n'.format(job_results_folder) + \
                                         'unzip -qq -o -d {0:s}-temp {0:s}-archive.zip;\n'.format(job_results_folder) + \
                                         'mv {0:s}-temp/icistarget {0:s};\n'.format(job_results_folder) + \
                                         'mv {0:s}-archive.zip {0:s}/archive.zip;\n'.format(job_results_folder) + \
                                         'rmdir {0:s}-temp;\n'.format(job_results_folder)

            ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command(extract_results_zip_script,
                                                                 bufsize=-1,
                                                                 timeout=None,
                                                                 get_pty=True)

            # Read stdout so the commands are executed.
            ssh_stdout_str = ssh_stdout.read()

            # Close SSH connection.
            ssh.close()
        except paramiko.BadHostKeyException, e:
            logging.error("SSH error: {0:s}.".format(e))
            raise AssertionError("Something when wrong while uploading the report results.")
        except paramiko.AuthenticationException, e:
            logging.error("SSH error: {0:s}.".format(e))
            raise AssertionError("Something when wrong while uploading the report results.")
        except paramiko.SSHException, e:
            logging.error("SSH error: {0:s}.".format(e))
            raise AssertionError("Something when wrong while uploading the report results.")
        except socket.error, e:
            logging.error("SSH error: {0:s}.".format(e))
            raise AssertionError("Something when wrong while uploading the report results.")

    logging.info("Uploading report results finished.")


def display_usage(output_fh=sys.stderr, cmd=sys.argv[0]):
    print >> output_fh, "Usage: python {0:s} <ini-file>".format(cmd)
    print >> output_fh
    print >> output_fh, "Example ini file:"
    print >> output_fh

    with open(INI_TEMPLATE, 'r') as input_fh:
        print >> output_fh, input_fh.read()


def parse_arguments(arguments):
    if len(arguments) == 1:
        ini_filename = arguments[0]
    else:
        display_usage(output_fh=sys.stderr)
        sys.exit(2)

    if not os.path.isfile(ini_filename):
        print >> sys.stderr, "'{0:s}' file doesn't exist.".format(ini_filename)
        display_usage(output_fh=sys.stderr)
        sys.exit(1)

    cfg = RawConfigParser()
    cfg.read(ini_filename)

    id_description_filename = (os.path.expandvars(cfg.get('compserver', 'id_description_filename'))
                               if cfg.has_option('compserver', 'id_description_filename') else None)
    id_location_bed_filename = (os.path.expandvars(cfg.get('compserver', 'id_location_bed_filename'))
                                if cfg.has_option('compserver', 'id_location_bed_filename') else None)
    track_cluster_filename = None
    stamp_command = cfg.get('STAMP', 'command')
    output_folder = os.path.expandvars(cfg.get('compserver', 'results_folder'))

    two_bit_to_fa_command = os.path.expandvars(cfg.get('twoBitToFa', 'command'))
    genome_2bit_filename = os.path.expandvars(cfg.get('twoBitToFa', 'genome_2bit_filename'))
    cluster_buster_command = os.path.expandvars(cfg.get('clusterBuster', 'command'))

    motif_annotations_filename = (os.path.expandvars(cfg.get('compserver', 'motif_annotations_filename'))
                                  if cfg.has_option('compserver', 'motif_annotations_filename') else None)

    track_annotations_filename = (os.path.expandvars(cfg.get('compserver', 'track_annotations_filename'))
                                  if cfg.has_option('compserver', 'track_annotations_filename') else None)

    db_cache_size = (cfg.getint('compserver', 'db_cache_size')
                     if cfg.has_option('compserver', 'db_cache_size') else 0)

    if cfg.has_option('compserver', 'rankings_metadata'):
        context = JobContext.create_from_zip_archive("regions",
                                                     id_description_filename,
                                                     id_location_bed_filename,
                                                     os.path.expandvars(cfg.get('compserver', 'rankings_metadata')),
                                                     os.path.expandvars(cfg.get('compserver', 'rankings_folder')),
                                                     dict(cfg.items('dbmapping')),
                                                     dict(cfg.items('dbnames')),
                                                     output_folder,
                                                     cfg.get('compserver', 'webserver_url').rstrip('/'),
                                                     cfg.get('compserver', 'webserver_upload_folder').strip('/'),
                                                     (cfg.get('compserver', 'webserver_email_address').strip()
                                                      if cfg.has_option('compserver', 'webserver_email_address')
                                                      else None),
                                                     (cfg.get('compserver', 'mail_server_host').strip()
                                                      if cfg.has_option('compserver', 'mail_server_host')
                                                      else None),
                                                     (int(cfg.get('compserver', 'mail_server_port').strip())
                                                      if cfg.has_option('compserver', 'mail_server_port')
                                                      else 25),
                                                     two_bit_to_fa_command,
                                                     genome_2bit_filename,
                                                     cluster_buster_command,
                                                     track_cluster_filename,
                                                     stamp_command,
                                                     db_cache_size,
                                                     motif_annotations_filename=motif_annotations_filename,
                                                     track_annotations_filename=track_annotations_filename)
    else:
        logo_folder = (os.path.expandvars(cfg.get('compserver', 'logo_folder'))
                       if cfg.has_option('compserver', 'logo_folder')
                       else None)
        logo_extension = (cfg.get('compserver', 'logo_extension')
                          if cfg.has_option('compserver', 'logo_extension')
                          else None)
        stamp_motif_database_filename = os.path.expandvars(cfg.get('STAMP', 'transfac_motif_database'))
        stamp_score_distribution_filename = os.path.expandvars(cfg.get('STAMP', 'score_database'))
        pwm_folder = os.path.expandvars(cfg.get('clusterBuster', 'pwm_folder'))
        pwm_extension = os.path.expandvars(cfg.get('clusterBuster', 'pwm_extension'))
        context = JobContext.create_from_files("regions",
                                               id_description_filename,
                                               id_location_bed_filename,
                                               logo_folder,
                                               logo_extension,
                                               os.path.expandvars(cfg.get('compserver', 'rankings_folder')),
                                               dict(cfg.items('dbmapping')),
                                               dict(cfg.items('dbnames')),
                                               output_folder,
                                               cfg.get('compserver', 'webserver_url').rstrip('/'),
                                               cfg.get('compserver', 'webserver_upload_folder').strip('/'),
                                               (cfg.get('compserver', 'webserver_email_address').strip()
                                                if cfg.has_option('compserver', 'webserver_email_address')
                                                else None),
                                               (cfg.get('compserver', 'mail_server_host').strip()
                                                if cfg.has_option('compserver', 'mail_server_host')
                                                else None),
                                               (int(cfg.get('compserver', 'mail_server_port').strip())
                                                if cfg.has_option('compserver', 'mail_server_port')
                                                else 25),
                                               two_bit_to_fa_command,
                                               genome_2bit_filename,
                                               cluster_buster_command,
                                               pwm_folder,
                                               pwm_extension,
                                               track_cluster_filename,
                                               stamp_command,
                                               stamp_motif_database_filename,
                                               stamp_score_distribution_filename,
                                               db_cache_size,
                                               motif_annotations_filename=motif_annotations_filename,
                                               track_annotations_filename=track_annotations_filename)

    return cfg, context


def init_logging(log_filename=None):
    output_fh = open(log_filename, 'a') if log_filename else sys.stderr
    if log_filename:
        # Make log file readable for everybody.
        os.chmod(log_filename, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH)
    logging.basicConfig(stream=output_fh, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)
    # Disable most logging output for the paramiko module.
    logging.getLogger('paramiko').setLevel(logging.WARNING)


def dump_stacktrace2log(job, administrator_email_address=None):
    exc_type, exc_value, exc_traceback = sys.exc_info()
    output = StringIO.StringIO()
    traceback.print_exception(exc_type, exc_value, exc_traceback, file=output)
    logging.error("error while running job {0:s}: {1:s}".format(job.name, output.getvalue()))

    if administrator_email_address:
        subject = 'Server error occurred for {0:s} (ID = {1:s})'.format(job.name, job.unique_id)
        message = (
            'A server error occurred for {0:s} (ID = {1:s})\n\n'.format(job.name, job.unique_id) +
            output.getvalue() + '\n\n' +
            ('If necessary, contact the submitter of the job: {0:s}\n\n'.format(job.email_address)
             if job.email_address
             else '')
        )
        job.mail_administrator(administrator_email_address, subject, message)


def dump_errormessage(job, message):
    try:
        job.error(message)
    except:
        pass

    try:
        job.mail_error()
    except:
        pass


def check_signal(signum, frame):
    if signum == signal.SIGUSR1:
        global halt_processing_jobs

        if halt_processing_jobs:
            halt_processing_jobs = False
            logging.info("Resume processing jobs.")
        else:
            halt_processing_jobs = True
            logging.info("Halt processing jobs after the current one is finished.")
    elif signum == signal.SIGUSR2:
        global stop_processing_jobs_and_exit
        stop_processing_jobs_and_exit = True

        logging.info("Stop processing jobs after the current one is finished and exit.")


def main():
    cfg, context = parse_arguments(sys.argv[1:])

    log_directory = (os.path.expandvars(cfg.get('compserver', 'log_dir'))
                     if cfg.has_option('compserver', 'log_dir')
                     else None)
    log_filename = (tempfile.mkstemp(suffix='.log',
                                     prefix=datetime.datetime.now().strftime('i-cistarget-daemon.%Y-%m-%d__%Hh%Mm%Ss.'),
                                     dir=log_directory)[1]
                    if log_directory
                    else None
                    )
    init_logging(log_filename)

    # Set signal handlers to control:
    #   - halting/resuming processing jobs: send SIGUSR1 signal
    #   - stopping the i-cisTarget daemon:  send SIGUSR2 signal
    # after running jobs are completed.
    signal.signal(signal.SIGUSR1, check_signal)
    signal.signal(signal.SIGUSR2, check_signal)

    global stop_processing_jobs_and_exit
    stop_processing_jobs_and_exit = False

    global halt_processing_jobs
    halt_processing_jobs = False

    # Connection to database is kept open all the time.
    # Opening and reopening a connection is too slow.
    with JobQueue(context,
                  cfg.get('dbserver', 'servername'),
                  cfg.getint('dbserver', 'port'),
                  cfg.get('dbserver', 'username'),
                  base64.b64decode(cfg.get('dbserver', 'password')),
                  cfg.get('dbserver', 'database')) as queue:
        for job in queue:
            if stop_processing_jobs_and_exit:
                logging.info("Stop processing jobs as requested and exit.")
                sys.exit(0)
            elif halt_processing_jobs is False:
                if not job:
                    time.sleep(3)
                else:
                    try:
                        # Check if all requested databases can be handled by this daemon.
                        if set(job.db_ids).issubset(context.db_id2name):
                            job.start()
                            logging.info("running job for {0:s} (ID = {1:s})".format(job.name, job.unique_id))
                            results_folder, subdirs_list, filenames_list = job.do_analysis()
                            upload_files(cfg, job.unique_id, results_folder, subdirs_list, filenames_list)

                            if os.path.isdir(results_folder):
                                shutil.rmtree(results_folder)

                            job.finish()
                            job.mail_success()
                            logging.info("finished job for {0:s} (ID = {1:s})".format(job.name, job.unique_id))
                        else:
                            queue.skipjob(job.job_id)
                    except AssertionError as e:
                        # First print information to the log, then insert error message into the database.
                        dump_stacktrace2log(job)
                        dump_errormessage(job, e.message)
                    except JobRunningInOtherInstanceError:
                        # The job was already running in the meantime by another daemon process.
                        pass
                    except:
                        administrator_email_address = cfg.get('compserver', 'administrator_email_address').strip()
                        # First print information to the log, then insert error message into the database.
                        dump_stacktrace2log(job, administrator_email_address)
                        dump_errormessage(job,
                                          ("A server error occurred during the analysis of your query. "
                                           "Please, send a mail with your job ID {0:s} to {1:s}.".format(
                                              job.unique_id,
                                              administrator_email_address))
                                          )

            time.sleep(1)


if __name__ == "__main__":
    main()
