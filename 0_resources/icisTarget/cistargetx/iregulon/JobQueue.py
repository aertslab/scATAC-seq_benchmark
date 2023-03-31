import MySQLdb

from Job import Job
from cistargetx.cistargetdb.bunch import Bunch


JOBS_QUERY_STATEMENT = r"""
  SELECT ID, name, nomenclatureCode, rankingDatabaseCode, conversionDelineation, conversionUpstreamRegionInBp,
         conversionDownstreamRegionInBp, conversionFractionOfOverlap, rankThreshold, AUCRankThresholdAsPercentage,
         NESThreshold, minOrthologousIdentity, maxMotifSimilarityFDR, geneIDs, ip, user_agent
  FROM jobQueue
  WHERE jobStatusCode = 1
  ORDER BY jobRequestTime;
"""

JOBS_COUNT_STATEMENT = r"""
  SELECT count(*)
  FROM jobQueue
  WHERE jobStatusCode = 1;
"""


class JobQueue(object):
    """ Implements iteratable, iterator and context manager protocol. """

    @staticmethod
    def _create_db_connection(servername, port, username, password, database):
        connection = MySQLdb.connect(host=servername,
                                     port=port,
                                     user=username,
                                     passwd=password,
                                     db=database)

        # Automatically reconnect to the MySQL server if the connection goes away.
        connection.ping(True)

        return connection

    def __init__(self, context, servername, port, username, password, database):
        self.servername = servername
        self.port = port
        self.username = username
        self.password = password
        self.database = database
        self.context = context

        self.configuration = Bunch(servername=servername,
                                   port=port,
                                   username=username,
                                   password=password,
                                   database=database)

    def __enter__(self):
        self.connection = JobQueue._create_db_connection(self.servername,
                                                         self.port,
                                                         self.username,
                                                         self.password,
                                                         self.database)
        return self

    def __exit__(self, type, value, traceback):
        self.connection.close()

        return False

    def __iter__(self):
        return self

    def fetch_jobs(self, max_jobs=1):
        cursor = self.connection.cursor()
        cursor.execute(JOBS_QUERY_STATEMENT)
        rows = cursor.fetchmany(size=max_jobs)

        if rows:
            jobs = [Job(self.connection, self.configuration, self.context, *row) for row in rows]
        else:
            jobs = []

        cursor.close()

        return jobs

    def next(self):
        jobs = self.fetch_jobs(max_jobs=1)

        return jobs[0] if jobs else None

    def __len__(self):
        cursor = self.connection.cursor()
        cursor.execute(JOBS_COUNT_STATEMENT)

        count = int(cursor.fetchone()[0])

        cursor.close()

        return count
