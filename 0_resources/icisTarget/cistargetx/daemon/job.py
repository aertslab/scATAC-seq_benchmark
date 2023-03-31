import Queue
import hashlib
import httplib
import itertools
import logging
import operator
import os
import re
import shutil
import smtplib
import urllib
from email.mime.text import MIMEText
from urlparse import urlparse

import MySQLdb

from cistargetx.common.rankingsdatabase import load_databases
from cistargetx.conversion.tools import gene_ids2regions, bed2regions, get_fraction_of_mapped_gene_ids, GeneIDType, \
    input_bed2input_bed_and_overlap_icistarget_regions
from cistargetx.conversion.featurelist import FeatureList
from cistargetx.recoveryanalysis.recoverycurves import RecoveryCurves
from cistargetx.motifscanning.motifscanner import scan
from cistargetx.recoveryanalysis.reportgeneration import generate_report as generate_recovery_analysis_report, \
    generate_batch_report, get_all_subdirs_and_filenames_for_directory
from cistargetx.motifscanning.reportgeneration import generate_report as generate_motif_scanning_report
from cistargetx.common.signatureformats import GeneSignature, LEGACY, GMT_COMMA, GMT_SEMICOLON, GMT_TAB, \
    GeneSignatureFileFormat


RCA_JOBS_QUERY_STATEMENT = r"""
  SELECT jobRequestTime, jobID, jobName, jobStatus, queryOriginalType, queryOriginalData, queryType,
         queryData, conversionParameters, analysisParameters, emailAddress, trueCRMlocationsBED
  FROM recoveryAnalysisQueue
  WHERE jobStatus = 'requested'
  ORDER BY jobRequestTime;
"""

RCA_START_JOB_STATEMENT = r"""
  UPDATE recoveryAnalysisQueue
  SET jobUID = '{0:s}', jobStatus = 'running', jobStartTime = now()
  WHERE jobID = {1:d} AND jobStatus = 'requested';
"""

RCA_STOP_JOB_STATEMENT = r"""
  UPDATE recoveryAnalysisQueue
  SET jobStatus = 'finished', jobFinishTime = now(),
      queryType = "{0:s}", queryData = "{1:s}", emailed = 1
  WHERE jobID = {2:d};
"""

RCA_ERROR_JOB_STATEMENT = r"""
  UPDATE recoveryAnalysisQueue
  SET jobStatus = 'error', jobFinishTime = now(), errorMessage = %s
  WHERE jobID = %s;
"""

RCA_DUMP_RESULTS_STATEMENT = r"""
  INSERT INTO recoveryAnalysisResults(jobID, featureID, candidateTargetIDs, topRankedIDs)
  VALUES({0:d}, "{1:s}", "{2:s}", "{3:s}");
"""

SCR_START_JOB_STATEMENT = r"""
  UPDATE clusterBusterScanQueue
  SET jobUID = '{0:s}', jobStatus = 'running', jobStartTime = now()
  WHERE jobID = {1:d} AND jobStatus = 'requested';
"""

SCR_STOP_JOB_STATEMENT = r"""
  UPDATE clusterBusterScanQueue
  SET jobStatus = 'finished', jobFinishTime = now(), emailed = 1
  WHERE jobID = {0:d};
"""

SCR_ERROR_JOB_STATEMENT = r"""
  UPDATE clusterBusterScanQueue
  SET jobStatus = 'error', jobFinishTime = now(), errorMessage = '{1:s}'
  WHERE jobID = {0:d};
"""

SCR_JOBS_QUERY_STATEMENT = r"""
  SELECT jobRequestTime, jobID, parentID, jobName, jobStatus, queryType, queryData, parameters, emailAddress
  FROM clusterBusterScanQueue
  WHERE jobStatus = 'requested'
  ORDER BY jobRequestTime;
"""

RCA_PARENT_ANALYSISPARAMETERS_SCR_JOBS_QUERY_STATEMENT = r"""
  SELECT analysisParameters
   FROM recoveryAnalysisQueue
   WHERE jobID = {0:d};
"""

REQUIRED_FILES = map(lambda f: os.path.join(os.path.dirname(__file__), f), ['1024px.css', 'background.gif'])

INPUT_MAPPED_TO_ICISTARGET_REGIONS_BED_FILENAME = 'input_mapped_to_icistarget_regions.bed'

FRACTION_OF_MAPPED_GENE_IDS = 0.80
MAX_NUMBER_OF_SIGNATURES = 250


class JobRunningInOtherInstanceError(Exception):
    pass


class JobStatus:
    REQUESTED = 'requested'
    RUNNING = 'running'
    DONE = 'finished'
    FAILED = 'error'


class Job(object):
    @staticmethod
    def _parse_parameters(parameters):
        if parameters:
            return dict([tuple(element.split('='))
                         for element in re.split('[\n\r]+', parameters)
                         if len(element.split('=')) == 2])
        else:
            return dict()

    def __init__(self, connection, context, request_time, job_id, name, status, email_address):
        self.connection = connection
        self.context = context
        self.request_time = request_time
        self.job_id = job_id
        self.name = unicode(name, errors='replace').encode('ascii', 'replace')
        self.status = status
        self.email_address = email_address

    def __str__(self):
        return "Job #{0:d} '{1:s}'".format(self.job_id, self.name)

    @property
    def unique_id(self):
        return str(self.job_id)

    def init_folder(self, results_folder):
        if not os.path.exists(results_folder):
            os.mkdir(results_folder)

        for filename in REQUIRED_FILES:
            shutil.copy(filename, results_folder)

    def start(self):
        pass

    def finish(self):
        pass

    def error(self, message=""):
        pass

    def do_analysis(self):
        pass

    @property
    def url(self):
        return "{0:s}/{1:s}/{2:s}/report.html".format(self.context.webserver_url,
                                                      self.context.webserver_upload_folder,
                                                      self.unique_id)

    def mail_error(self):
        self.mail_user_via_webserver()

    def mail_success(self):
        self.mail_user_via_webserver()

    def mail_user(self, msg):
        if self.email_address and self.context.webserver_email_address:
            s = smtplib.SMTP('localhost')
            s.sendmail(self.context.webserver_email_address, [self.email_address], msg.as_string())
            s.quit()

    def mail_user_via_webserver(self):
        pass

    def _mail_user_via_webserver(self, queue_type):
        if self.email_address:
            params = urllib.urlencode({'job_id': self.job_id, 'queue_type': queue_type})
            headers = {"Content-type": "application/x-www-form-urlencoded", "Accept": "text/plain"}
            parsed_url = urlparse(self.context.webserver_url)
            HTTPOrHTTPSConnection = (httplib.HTTPSConnection
                                     if parsed_url.scheme == 'https'
                                     else httplib.HTTPConnection)
            connection = HTTPOrHTTPSConnection(parsed_url.hostname)
            connection.request("POST", parsed_url.path + "/mail.php", params, headers)
            response = connection.getresponse()

            if response.status != 200:
                logging.error('Sending e-mail via PHP server resulted in an error: {0:s}.'.format(response.reason))
            else:
                logging.info('Sending e-mail via PHP server: response from server = "{0:s}".'.format(response.read()))

            connection.close()

    def mail_administrator(self, administrator_email_address, subject, message):
        if administrator_email_address and self.context.webserver_email_address and self.context.mail_server_host:
            msg = MIMEText(message, 'plain')
            msg['Subject'] = subject
            msg['From'] = self.context.webserver_email_address
            msg['To'] = administrator_email_address

            s = smtplib.SMTP(host=self.context.mail_server_host, port=self.context.mail_server_port)
            s.sendmail(msg['From'], [msg['To']], msg.as_string())
            s.quit()


def _strip_blanks(string):
    string = re.sub('[ \t]+$', '', string)

    return re.sub('^[ \t]+', '', string)


def _is_not_empty(string):
    return len(string) > 0


def memoize(fnc):
    cache = {}

    def _helper(*args):
        if args not in cache:
            cache[args] = fnc(*args)

        return cache[args]

    return _helper


PARSER2FORMAT = {
    GMT_TAB: 'gmt_tab',
    GMT_COMMA: 'gmt_comma',
    GMT_SEMICOLON: 'gmt_semicolon',
    LEGACY: 'legacy'}

FORMAT2PARSER = {
    'gmt_tab': GMT_TAB,
    'gmt_comma': GMT_COMMA,
    'gmt_semicolon': GMT_SEMICOLON,
    'legacy': LEGACY
}


class RecoveryAnalysisJob(Job):
    def __init__(self, connection, context, request_time, job_id, name, status, query_original_type, query_original_data,
                 query_type, query_data, conversion_parameters, analysis_parameters, email_address, true_crm_locations):
        super(RecoveryAnalysisJob, self).__init__(connection, context, request_time, job_id, name, status, email_address)
        self.query_original_type = query_original_type
        self.query_original_data = query_original_data
        self.query_data = filter(_is_not_empty, query_data.split('\n')) if query_data else None
        self.conversion_parameters = Job._parse_parameters(conversion_parameters)
        self.analysis_parameters = Job._parse_parameters(analysis_parameters)
        self.true_crm_locations = true_crm_locations

    def mail_user_via_webserver(self):
        self._mail_user_via_webserver('recoveryanalysis')

    # ###################################################################################################################
    # Database related operations ...
    ####################################################################################################################

    def start(self):
        cursor = self.connection.cursor()
        cursor.execute(RCA_START_JOB_STATEMENT.format(self.unique_id, self.job_id))
        self.connection.commit()
        job_already_running = cursor.rowcount == 0
        cursor.close()

        if job_already_running:
            # Skip this job if it was already started in the meanwhile by another daemon process.
            raise JobRunningInOtherInstanceError

    def finish(self):
        cursor = self.connection.cursor()
        cursor.execute(
            RCA_STOP_JOB_STATEMENT.format(self.query_type_as_str, self.converted_query_data_as_str, self.job_id))
        self.connection.commit()
        cursor.close()

    def error(self, message):
        cursor = self.connection.cursor()
        cursor.execute(RCA_ERROR_JOB_STATEMENT, (message, self.job_id))
        self.connection.commit()
        cursor.close()

        results_folder = os.path.join(self.context.results_folder, self.unique_id)

        if os.path.exists(results_folder):
            shutil.rmtree(results_folder)

    ####################################################################################################################
    # Primary properties ...
    ####################################################################################################################

    @property
    def assembly(self):
        return self.analysis_parameters['assembly']

    @property
    def gene_annotation_version(self):
        return self.analysis_parameters['gene_annotation_version']

    @property
    def delineation(self):
        return self.conversion_parameters['delineation']

    @property
    def overlap_fraction(self):
        if 'overlap_fraction' not in self.conversion_parameters or self.conversion_parameters['overlap_fraction'] == '':
            return None

        try:
            return float(self.conversion_parameters['overlap_fraction'])
        except ValueError:
            raise AssertionError("Invalid value for minimum fraction of overlap.")

    @property
    def db_ids(self):
        return self.analysis_parameters['ranking_dbs'].split(';')

    @property
    def control_db_id(self):
        return self.analysis_parameters['control_db']

    @property
    def nes_threshold(self):
        try:
            return float(re.sub(',', '.', self.analysis_parameters['nes_threshold']))
        except ValueError:
            raise AssertionError("Invalid value for NES threshold.")

    @property
    def enrichment_within_db(self):
        return self.analysis_parameters['enrichment_within_db'].lower() == "true"

    @property
    def auc_threshold(self):
        try:
            return float(re.sub(',', '.', self.analysis_parameters['auc_threshold']))
        except ValueError:
            raise AssertionError("Invalid value for AUC threshold.")

    @property
    def rank_threshold(self):
        try:
            return int(self.analysis_parameters['rank_threshold'])
        except ValueError:
            raise AssertionError("Invalid value for rank threshold.")

    @property
    def feature_masks(self):
        value = self.analysis_parameters['feature_masks']

        if not value or value == '*':
            return None

        return value.split(';')

    @property
    def feature_combinations(self):
        combinations = self.analysis_parameters['feature_combinations']

        return combinations.split(';') if combinations else None

    @property
    def combination_method(self):
        return self.analysis_parameters['combination_method']

    @property
    def true_regions(self):
        if not self.true_crm_locations:
            return []
        else:
            return bed2regions(self.assembly, self.gene_annotation_version, self.true_crm_locations, fraction=self.overlap_fraction)

    ####################################################################################################################
    # Derived properties ...
    ####################################################################################################################

    @property
    def unique_id(self):
        return hashlib.sha1("recovery_analysis_" + str(self.job_id)).hexdigest()

    @property
    @memoize
    def query_type(self):
        """
        Returns tuple (format, nomenclature, fraction).

          - Format can be: bed, gmt_tab, gmt_comma, gmt_semicolon, legacy.
          - Nomenclature can be: region, fbgn, cg, symbol or None (for bed format).
          - Fraction can be None.
        """
        supplied_info = _strip_blanks(self.query_original_type)
        format_type = self._get_format(supplied_info)

        if format_type == "bed":
            return format_type, "region", None

        nomenclature, fraction = self._get_nomenclature(format_type, supplied_info)

        return format_type, nomenclature, fraction

    def _parser2format(self, parser):
        return PARSER2FORMAT.get(parser)

    def _format2parser(self, format_type):
        return FORMAT2PARSER.get(format_type)

    @property
    @memoize
    def _data_lines(self):
        return filter(_is_not_empty, map(_strip_blanks, re.split('[\n\r]+', self.query_original_data)))

    def _get_format(self, supplied_info):
        if supplied_info == "bed":
            return "bed"
        else:
            return self._auto_detect_format()

    def _auto_detect_format(self):
        if FeatureList.is_bed_file(self.query_original_data):
            return 'bed'

        return self._parser2format(GeneSignatureFileFormat.guess_file_format(self._data_lines))

    def _get_gene_ids(self, format_type):
        def combine(s1, s2):
            return s1 | s2

        return reduce(combine,
                      (signature.identifiers for signature in GeneSignature._load_from_lines(
                          self._data_lines,
                          self.name,
                          gene_signature_file_format=self._format2parser(format_type))),
                      set())

    def _get_nomenclature(self, format_type, supplied_info):
        gene_ids = self._get_gene_ids(format_type)

        if supplied_info in ("symbol", "fbgn", "cg", "hgnc_symbol", "mgi_symbol"):
            fraction = get_fraction_of_mapped_gene_ids(self.assembly, self.gene_annotation_version, gene_ids, supplied_info, self.delineation)

            return supplied_info, fraction
        elif supplied_info == "region":
            return supplied_info, 1.0
        else:
            return self._auto_detect_nomenclature(gene_ids)

    def _auto_detect_nomenclature(self, gene_ids):
        if self._all_cgs(gene_ids):
            fraction = get_fraction_of_mapped_gene_ids(self.assembly, self.gene_annotation_version, gene_ids, "cg", self.delineation)

            return "cg", fraction
        elif self._all_fbgns(gene_ids):
            fraction = get_fraction_of_mapped_gene_ids(self.assembly, self.gene_annotation_version, gene_ids, "fbgn", self.delineation)

            return "fbgn", fraction
        else:
            gene_id_types = (GeneIDType.FBGN, GeneIDType.SYMBOL, GeneIDType.CG, GeneIDType.HGNC_SYMBOL, GeneIDType.MGI_SYMBOL)
            fractions = [(gene_id_type, get_fraction_of_mapped_gene_ids(self.assembly, self.gene_annotation_version, gene_ids, gene_id_type, self.delineation))
                         for gene_id_type in gene_id_types]
            best_gene_id_type, best_fraction = sorted(fractions, key=operator.itemgetter(1), reverse=True)[0]

            assert best_fraction > FRACTION_OF_MAPPED_GENE_IDS, "Could not detect type of supplied gene IDs."

            return best_gene_id_type, best_fraction

    def _all_cgs(self, gene_ids):
        return all(gene_id.startswith("CG") for gene_id in gene_ids)

    def _all_fbgns(self, gene_ids):
        return all(gene_id.startswith("FBgn") for gene_id in gene_ids)

    @property
    def query_match_fraction(self):
        """ Returns the matched number of IDs. If the query type could not be derived None is returned. """
        try:
            return self.query_type[2]
        except AssertionError:
            return None
        except KeyError:
            return None

    @property
    def query_type_as_str(self):
        """
        Returns the derived query type that is stored in the database and displayed in the report.
        Current supported output is "bed", "region", "fbgn", "cg" and "symbol". Returning "auto" would
        not be informative. If the query type could not be derived None is returned.
        """
        try:
            format_type, nomenclature = self.query_type[:2]
            if format_type == 'bed':
                return 'bed'
            else:
                return nomenclature
        except AssertionError:
            return None
        except KeyError:
            return None

    @property
    @memoize
    def original_query_data(self):
        """ Returns [(name, {IDs})]. If BED file was supplied this BED file is already converted to region IDs. """
        format_type, nomenclature, fraction = self.query_type

        if format_type == 'bed':
            regions = bed2regions(self.assembly, self.gene_annotation_version, self.query_original_data, fraction=self.overlap_fraction)

            assert len(regions) > 0, (
                "Not a single peak in the supplied BED file could be mapped to a i-cisTarget region. "
                "Please, check carefully that the correct species and input type were selected."
            )

            return [(self.name, regions)]
        else:
            assert fraction > FRACTION_OF_MAPPED_GENE_IDS, (
                "Less than {0:.0f}% of genes could be mapped. ".format(FRACTION_OF_MAPPED_GENE_IDS * 100) +
                "Please, check carefully that the correct species and input type were selected."
            )

            def _generate_tuples():
                for signature in GeneSignature._load_from_lines(self._data_lines, self.name,
                                                                gene_signature_file_format=self._format2parser(
                                                                        format_type)):
                    assert len(signature.identifiers) > 1, (
                        "An empty or singleton signature was supplied. i-cisTarget only works on sets of gene IDs."
                    )

                    yield signature.name, signature.identifiers

            return [t for t in _generate_tuples()]

    #@property
    #@memoize
    def write_bed_original_query_data_overlap_with_icistarget_regions(self):
        """ Returns [(name, {IDs})]. If BED file was supplied this BED file is already converted to region IDs. """
        format_type, nomenclature, fraction = self.query_type

        if format_type == 'bed':
            input_mapped_to_icistarget_regions_bed_filename = os.path.join(self.context.results_folder,
                                                                           self.unique_id,
                                                                           INPUT_MAPPED_TO_ICISTARGET_REGIONS_BED_FILENAME)

            with open(input_mapped_to_icistarget_regions_bed_filename, 'w') as fh:
                for input_bed_region_overlap_icistarget_region in input_bed2input_bed_and_overlap_icistarget_regions(
                        self.assembly,
                        self.gene_annotation_version,
                        self.query_original_data,
                        fraction=self.overlap_fraction):
                    fh.write(input_bed_region_overlap_icistarget_region + "\n")

    @property
    @memoize
    def converted_query_data(self):
        """ Returns [(name, {regionIDs})]. """
        format_type, nomenclature = self.query_type[:2]

        if format_type == 'bed' or nomenclature == 'region':
            return self.original_query_data

        def _generate_tuples():
            for name, gene_ids in self.original_query_data:
                region_ids = gene_ids2regions(self.assembly,
                                              self.gene_annotation_version,
                                              gene_ids,
                                              gene_id_type=nomenclature,
                                              delineation=self.delineation,
                                              fraction=self.overlap_fraction)

                assert len(region_ids) > 0, "Not a single valid gene identifier was supplied."

                yield name, region_ids

        return [t for t in _generate_tuples()]

    @property
    @memoize
    def first_signature_ids(self):
        """
        This property returns the signature IDs of the first supplied signature as [IDs].
        If multiple signatures or a BED file was supplied then None is returned.
        """
        format_type, nomenclature = self.query_type[:2]

        if format_type == 'bed' or nomenclature not in ('fbgn', 'symbol') or len(self.original_query_data) > 1:
            return None
        else:
            return self.original_query_data[0][1]

    @property
    def converted_query_data_as_str(self):
        """
        The final result of conversion is a list of regions IDs when a BED file or one signature was
        supplied. Because further analysis must still be supported this must be stored in the database
        as a string where IDs are separated by a new line character.
        """
        if len(self.converted_query_data) == 1:
            return "\n".join(self.converted_query_data[0][1])
        else:
            return None

    ####################################################################################################################
    # Actual analysis method ...
    ####################################################################################################################

    def do_analysis(self):
        context = self.context

        # Setting up results folder.
        results_folder = os.path.join(context.results_folder, self.unique_id)
        self.init_folder(results_folder)

        self.write_bed_original_query_data_overlap_with_icistarget_regions()

        # Calculating all necessary recovery curves.
        name2curves = dict()

        for name, signature in self.converted_query_data:
            name2curves[name] = [RecoveryCurves(db,
                                                self.auc_threshold,
                                                self.rank_threshold,
                                                self.nes_threshold)
                                 for db in self._databases_for_signature(signature)]

        # Write out report and all necessary auxiliary files.
        if len(name2curves) == 1:
            curves = name2curves.values()[0]
            generate_recovery_analysis_report(curves, context, results_folder,
                                              title=self.name, elements_to_highlight=self.first_signature_ids,
                                              template_variables=self._template_variables(),
                                              true_target_ids=self.true_regions)

            # Dump results to database for further analysis by user.
            self._write_results_to_db(curves)
        elif len(name2curves) <= MAX_NUMBER_OF_SIGNATURES:
            generate_batch_report(name2curves, results_folder, context, template_variables=self._template_variables())
        else:
            raise AssertionError("You provided more signatures than the maximum allowed number ({0:d}).".format(
                MAX_NUMBER_OF_SIGNATURES))

        # Return results folder with all subdirs and filenames (results_folder, subdirs_list, filenames_list) tuple.
        return get_all_subdirs_and_filenames_for_directory(results_folder)

    def _write_results_to_db(self, curves):
        cursor = self.connection.cursor()

        def convert(rank_ID_list):
            return "\n".join(map(operator.itemgetter(1), rank_ID_list))

        iterators = map(RecoveryCurves.iterate_enriched_features, curves)
        for feature_id, auc, nes, curve in itertools.chain(*iterators):
            cursor.execute(RCA_DUMP_RESULTS_STATEMENT.format(self.job_id, feature_id,
                                                             convert(curve.get_candidate_target_gene_ids(feature_id)),
                                                             convert(curve.get_top_ranked_gene_ids(feature_id))))
        self.connection.commit()
        cursor.close()

    def _databases_for_signature(self, signature):
        return load_databases([self.context.get_full_path_for_db_id(db_id) for db_id in self.db_ids],
                              gene_ids=signature,
                              feature_masks=None,
                              features_to_combine=self.feature_combinations,
                              enrichment_within_db=self.enrichment_within_db,
                              control_db=self.context.get_full_path_for_db_id(self.control_db_id),
                              combination_method=self.combination_method,
                              name_translation=self.context.get_description_for_db_filename,
                              db_cache_size=self.context.db_cache_size)

    def _template_variables(self):
        return {
            'job_id': self.unique_id,
            'job_id_number': self.job_id,
            'databases': list({'id': db_id, 'name': db_name}
                              for db_id, db_name in self.context.db_id2name.iteritems()
                              if db_id in self.db_ids),
            'webserver_url': self.context.webserver_url,
            'webserver_reports_folder': self.context.webserver_upload_folder,
            'input_query_type': self.query_type_as_str,
            'input_mapped_to_icistarget_regions_bed': INPUT_MAPPED_TO_ICISTARGET_REGIONS_BED_FILENAME if self.query_type_as_str else None,
            'input_query_match': format(self.query_match_fraction, '.3f') if self.query_match_fraction else None,
            'overlap_fraction': self.overlap_fraction,
            'assembly': self.assembly
        }


class ScanRegionsJob(Job):
    def __init__(self, connection, context, request_time, job_id, parent_job_id, name, status, query_type, query_data,
                 analysis_parameters, email_address):
        super(ScanRegionsJob, self).__init__(connection, context, request_time, job_id, name, status, email_address)
        self.parent_job_id = parent_job_id
        self.query_type = query_type
        self.query_data = query_data
        self.analysis_parameters = Job._parse_parameters(analysis_parameters)

        cursor = self.connection.cursor()
        cursor.execute(RCA_PARENT_ANALYSISPARAMETERS_SCR_JOBS_QUERY_STATEMENT.format(self.parent_job_id))
        rca_analysis_parameters = cursor.fetchone()[0]
        self.analysis_parameters_RCA = Job._parse_parameters(rca_analysis_parameters)
        self.connection.commit()
        cursor.close()

    @property
    def db_ids(self):
        return self.analysis_parameters_RCA['ranking_dbs'].split(';')

    def mail_user_via_webserver(self):
        self._mail_user_via_webserver('cbustscan')

    def start(self):
        cursor = self.connection.cursor()
        cursor.execute(SCR_START_JOB_STATEMENT.format(self.unique_id, self.job_id))
        self.connection.commit()
        job_already_running = cursor.rowcount == 0
        cursor.close()

        if job_already_running:
            # Skip this job if it was already started in the meanwhile by another daemon process.
            raise JobRunningInOtherInstanceError

    def finish(self):
        cursor = self.connection.cursor()
        cursor.execute(SCR_STOP_JOB_STATEMENT.format(self.job_id))
        self.connection.commit()
        cursor.close()

    def error(self, message):
        cursor = self.connection.cursor()
        cursor.execute(SCR_ERROR_JOB_STATEMENT.format(self.job_id, message))
        self.connection.commit()
        cursor.close()

        results_folder = os.path.join(self.context.results_folder, self.unique_id)
        if os.path.exists(results_folder): shutil.rmtree(results_folder)

    @property
    def assembly(self):
        return self.analysis_parameters['assembly']

    @property
    def crm_type(self):
        return self.analysis_parameters['crm_type']

    @property
    def motif_ids(self):
        return self.analysis_parameters['motif_ids'].split(';')

    @property
    def region_ids(self):
        return set(map(lambda s: s.strip(), re.split('[\r\n]+', self.query_data)))

    @property
    def unique_id(self):
        return hashlib.sha1("motif_scan_" + str(self.job_id)).hexdigest()

    @property
    def parent_unique_id(self):
        return hashlib.sha1("recovery_analysis_" + str(self.parent_job_id)).hexdigest()

    @property
    def template_variables(self):
        return {
            'job_id': self.unique_id,
            'parent_job_id': self.parent_unique_id,
            'job_id_number': self.job_id,
            'webserver_url': self.context.webserver_url,
            'webserver_reports_folder': self.context.webserver_upload_folder,
            'assembly': self.assembly
        }

    def do_analysis(self):
        locations = self.context.gene_or_region_id_locations.filter_by_name(self.region_ids)

        if self.crm_type == 'hetero':
            motifs = self.context.motif_collection.reduce(self.motif_ids)
            feature_iterator = scan(locations, motifs, self.context.cluster_buster_command,
                                    self.context.two_bit_to_fa_command, self.context.genome_2bit_filename)
        else:
            feature_iterators = []
            for motif_id in self.motif_ids:
                motifs = self.context.motif_collection.reduce([motif_id])
                feature_iterators.append(
                    scan(locations, motifs, self.context.cluster_buster_command, self.context.two_bit_to_fa_command,
                         self.context.genome_2bit_filename))
            feature_iterator = itertools.chain(*feature_iterators)

        results_folder = os.path.join(self.context.results_folder, self.unique_id)
        self.init_folder(results_folder)

        generate_motif_scanning_report(results_folder, feature_iterator, self.context.motif_scan_template_filename,
                                       self.template_variables)

        # Return results folder with all subdirs and filenames (results_folder, subdirs_list, filenames_list) tuple.
        return get_all_subdirs_and_filenames_for_directory(results_folder)


class JobQueue:
    """ Implements iterator protocol and context manager protocol. """

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
        self.queue = Queue.PriorityQueue()

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

    def _fetch_recovery_analysis_job(self, job_id=None):
        cursor = self.connection.cursor()
        cursor.execute(RCA_JOBS_QUERY_STATEMENT)
        row = cursor.fetchone()
        # Skip jobs not suited for the current daemon.
        if job_id:
            while row is not None and int(row[1]) <= int(job_id):
                row = cursor.fetchone()
        job = RecoveryAnalysisJob(self.connection, self.context, *row) if row else None
        self.connection.commit()
        cursor.close()

        return job

    def _fetch_motif_scanning_job(self, job_id=None):
        cursor = self.connection.cursor()
        cursor.execute(SCR_JOBS_QUERY_STATEMENT)
        row = cursor.fetchone()
        # Skip jobs not suited for the current daemon.
        if job_id:
            while row is not None and int(row[1]) <= int(job_id):
                row = cursor.fetchone()
        job = ScanRegionsJob(self.connection, self.context, *row) if row else None
        self.connection.commit()
        cursor.close()

        return job

    def next(self):
        if not self.queue.empty():
            return self.queue.get()[1]

        rca_job = self._fetch_recovery_analysis_job()
        scr_job = self._fetch_motif_scanning_job()

        if not rca_job:
            return scr_job
        elif not scr_job:
            return rca_job
        elif rca_job.request_time > scr_job.request_time:
            self.queue.put((rca_job.request_time, rca_job))
            return scr_job
        else:
            self.queue.put((scr_job.request_time, scr_job))
            return rca_job

    def skipjob(self, job_id):
        rca_job = self._fetch_recovery_analysis_job(job_id)
        scr_job = self._fetch_motif_scanning_job(job_id)

        if scr_job:
            self.queue.put((scr_job.request_time, scr_job))
        if rca_job:
            self.queue.put((rca_job.request_time, rca_job))
