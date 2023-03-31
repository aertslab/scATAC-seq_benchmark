import StringIO
import base64
import datetime
import logging
import os
import signal
import stat
import sys
import tempfile
import time
import traceback
from ConfigParser import RawConfigParser as ConfigParser
from multiprocessing import Process

if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

from Job import JobRunningInOtherInstanceError, UnsupportedRankingsDatabaseError
from JobContext import JobContext
from JobQueue import JobQueue


INI_TEMPLATE = os.path.join(os.path.dirname(__file__), 'template.ini')
DDL_FILENAME = os.path.join(os.path.dirname(__file__), 'iregulon.sql')

QUERY_NOMENCLATURE_CODE_STATEMENT = r"""
  SELECT code, name
  FROM geneNomenclature;
"""

# Maximum error message length which can be written to errorMessage (varchar(40960)
# field in iRegulon daemon MySQL database.
MAX_ERROR_MESSAGE_LENGTH = 4096


def display_usage(output_fh=sys.stderr, cmd=sys.argv[0]):
    print >> output_fh, "Usage: python {0:s} <ini-file>".format(cmd)
    print >> output_fh, "Usage: python {0:s} <ini-file> init".format(cmd)
    print >> output_fh
    print >> output_fh, "Example ini file:"
    print >> output_fh

    with open(INI_TEMPLATE, 'r') as input_fh:
        print >> output_fh, input_fh.read()


def parse_arguments(arguments):
    if len(arguments) == 1:
        mode = "run"
    elif len(arguments) == 2 and arguments[1] == "init":
        mode = "init"
    else:
        display_usage(output_fh=sys.stderr)
        sys.exit(2)

    ini_filename = arguments[0]

    if not os.path.isfile(ini_filename):
        print >> sys.stderr, "'{0:s}' file doesn't exist.".format(ini_filename)
        display_usage(output_fh=sys.stderr)
        sys.exit(1)

    cfg = ConfigParser()
    cfg.read(ini_filename)

    return mode, cfg, JobContext.create_from_cfg_file(cfg)


def init_logging(log_filename=None):
    output_fh = open(log_filename, 'a') if log_filename else sys.stderr
    if log_filename:
        # Make log file readable for everybody.
        os.chmod(log_filename, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH)
    logging.basicConfig(stream=output_fh, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)
    logging.captureWarnings(True)


def dump_stacktrace2log(job=None):
    exc_type, exc_value, exc_traceback = sys.exc_info()
    output = StringIO.StringIO()
    traceback.print_exception(exc_type, exc_value, exc_traceback, file=output)

    if job:
        logging.error("error while running job {0:s}: {1:s}".format(job.name, output.getvalue()))


def dump_errormessage(job, message):
    try:
        # Shorten error message if necessary, so it can be written to errorMessage
        # field in MySQL iRegulon daemon.
        job.error(message[0:MAX_ERROR_MESSAGE_LENGTH])
    except:
        pass


def load_ddl():
    with open(DDL_FILENAME, 'r') as input_fh:
        return input_fh.read()


def init_database(connection):
    cursor = connection.cursor()
    cursor.execute(load_ddl())
    connection.commit()
    cursor.close()


def load_nomenclature2code(connection):
    cursor = connection.cursor()
    cursor.execute(QUERY_NOMENCLATURE_CODE_STATEMENT)
    result = dict((nomenclature, code) for code, nomenclature in cursor)
    cursor.close()
    return result


def check_signal(signum, frame):
    if signum == signal.SIGUSR1:
        global halt_processing_jobs

        if halt_processing_jobs:
            halt_processing_jobs = False
            logging.info("Resume processing jobs.")
        else:
            halt_processing_jobs = True
            logging.info("Halt processing jobs after the current ones are finished.")
    elif signum == signal.SIGUSR2:
        global stop_processing_jobs_and_exit
        stop_processing_jobs_and_exit = True

        logging.info("Stop processing jobs after the current ones are finished and exit.")


def single_core_daemon(cfg, context):
    # Connection to database is kept open all the time.
    # Opening and reopening a connection is too slow.
    while True:
        try:
            with JobQueue(context,
                          cfg.get('dbserver', 'servername'),
                          cfg.getint('dbserver', 'port'),
                          cfg.get('dbserver', 'username'),
                          base64.b64decode(cfg.get('dbserver', 'password')),
                          cfg.get('dbserver', 'database')) as queue:
                for job in queue:
                    if stop_processing_jobs_and_exit:
                        sys.exit(0)
                    elif halt_processing_jobs is False:
                        if not job:
                            time.sleep(5)
                        else:
                            process_job(job)

                    time.sleep(1)
        except:
            dump_stacktrace2log()


def multiple_core_daemon(cfg, context, n_cores, niceness):
    # Connection to database is kept open all the time.
    # Opening and reopening a connection is too slow.
    n_consumers = n_cores - 1
    running_processes = []

    # One of the imported python modules sets the CPU affinity to one core.
    # If we do not recent the CPU affinity all processes (python and cbust) wil run
    # on the same core, defeating the whole purpose of the multiprocessing module.
    #
    # Most likely one of the imported modules in linked with OpenBLAS:
    #     http://stackoverflow.com/questions/23537716/importing-scipy-breaks-multiprocessing-support-in-python
    #
    # If we ever upgrade to Python 3, the call to taskset can be replaced with:
    #     os.sched_setaffinity(pid, mask)
    DEVNULL = open(os.devnull, 'w')
    subprocess.call(['taskset','-p', '0xFFFFFFFF', str(os.getpid())], stdout=DEVNULL)
    DEVNULL.close()

    while True:
        try:
            with JobQueue(context,
                          cfg.get('dbserver', 'servername'),
                          cfg.getint('dbserver', 'port'),
                          cfg.get('dbserver', 'username'),
                          base64.b64decode(cfg.get('dbserver', 'password')),
                          cfg.get('dbserver', 'database')) as queue:
                while True:
                    running_processes = filter(lambda p: p.is_alive(), running_processes)
                    idle_processes = n_consumers - len(running_processes)

                    if running_processes == 0:
                        if stop_processing_jobs_and_exit:
                            sys.exit(0)
                    elif halt_processing_jobs is False:
                        if idle_processes:
                            for job in queue.fetch_jobs(idle_processes):
                                job.connection = None
                                p = Process(target=process_job, args=(job, niceness, True))
                                running_processes.append(p)
                                p.start()

                                time.sleep(1)

                    time.sleep(3)
        except:
            dump_stacktrace2log()

            for p in running_processes:
                if p.is_alive:
                    p.terminate()


def process_job(job, niceness=1, new_connection=False):
    os.nice(niceness)

    try:
        if new_connection: job.reset()
        job.start()
        logging.info("running job: {0:s} (ID = {1:d})".format(job.name, job.job_id))
        logging.info("client info: ID = {0:d} submitted from {1:s} [{2:s}]".format(job.job_id, job.ip, job.user_agent))
        job.do_analysis()
        job.finish()
        logging.info("finished job: {0:s} (ID = {1:d})".format(job.name, job.job_id))
    except AssertionError as e:
        # First print information to the log, then insert error message into the database.
        dump_stacktrace2log(job)
        dump_errormessage(job, e.message)
    except UnsupportedRankingsDatabaseError:
        # Do not try to run this job when the database is not supported (needs to processed by the other daemon).
        pass
    except JobRunningInOtherInstanceError:
        # The job was already running in the meantime by another worker thread.
        pass
    except:
        # First print information to the log, then insert error message into the database.
        dump_stacktrace2log(job)
        dump_errormessage(job, "An error occurred during analysis of your query.")


def main():
    mode, cfg, context = parse_arguments(sys.argv[1:])

    if mode == "run":
        log_directory = (os.path.expandvars(cfg.get('compserver', 'log_dir'))
                         if cfg.has_option('compserver', 'log_dir')
                         else None)
        log_filename = (tempfile.mkstemp(suffix='.log',
                                         prefix=datetime.datetime.now().strftime('iregulon-daemon.%Y-%m-%d__%Hh%Mm%Ss.'),
                                         dir=log_directory)[1]
                        if log_directory
                        else None
                        )
        init_logging(log_filename)

        # Multiprocessing parameters ...
        n_cores = (cfg.getint('compserver', 'n_cores')
                   if cfg.has_option('compserver', 'n_cores')
                   else 1)
        niceness = (cfg.getint('compserver', 'process_niceness')
                    if cfg.has_option('compserver', 'process_niceness')
                    else 19)

        # Set signal handlers to control:
        #   - halting/resuming processing jobs: send SIGUSR1 signal
        #   - stopping the iRegulon daemon:     send SIGUSR2 signal
        # after running jobs are completed.
        signal.signal(signal.SIGUSR1, check_signal)
        signal.signal(signal.SIGUSR2, check_signal)

        global stop_processing_jobs_and_exit
        stop_processing_jobs_and_exit = False

        global halt_processing_jobs
        halt_processing_jobs = False

        if n_cores == 1:
            single_core_daemon(cfg, context)
        else:
            multiple_core_daemon(cfg, context, n_cores, niceness)
    elif mode == "init":
        connection = JobQueue._create_db_connection(cfg.get('dbserver', 'servername'),
                                                    cfg.getint('dbserver', 'port'),
                                                    cfg.get('dbserver', 'username'),
                                                    base64.b64decode(cfg.get('dbserver', 'password')),
                                                    cfg.get('dbserver', 'database'))
        init_database(connection)
        connection.close()


if __name__ == "__main__":
    main()
