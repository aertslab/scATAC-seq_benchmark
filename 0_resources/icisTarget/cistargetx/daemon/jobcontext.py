import os

from cistargetx.motifscanning.motifscanner import MotifCollection
from cistargetx.recoveryanalysis.reportcontext import ReportContext

TEMPLATE_WEB_REGIONS_MAIN = os.path.join(os.path.dirname(__file__), 'web-regions-main.vm')
TEMPLATE_WEB_REGIONS_CANDIDATE_TARGETS = os.path.join(os.path.dirname(__file__), 'web-regions-candidate-targets.vm')
TEMPLATE_WEB_REGIONS_TOP_RANKED_TARGETS = os.path.join(os.path.dirname(__file__), 'web-regions-top-ranked-targets.vm')
TEMPLATE_WEB_REGIONS_MOTIF_SCAN = os.path.join(os.path.dirname(__file__), 'web-regions-motif-scan.vm')
TEMPLATE_WEB_REGIONS_BATCH = os.path.join(os.path.dirname(__file__), 'web-regions-batch-report.vm')


class JobContext(ReportContext):
    @staticmethod
    def create_from_zip_archive(id_type,
                                id_description_filename,
                                id_location_bed_filename,
                                rankings_metadata_zip_filename,
                                rankings_folder,
                                db_id2filename,
                                db_id2name,
                                results_folder,
                                webserver_url,
                                webserver_upload_folder,
                                webserver_email_address,
                                mail_server_host,
                                mail_server_port,
                                two_bit_to_fa_command,
                                genome_2bit_filename,
                                cluster_buster_command,
                                track_cluster_filename=None,
                                stamp_command=None,
                                db_cache_size=0,
                                motif_annotations_filename=None,
                                track_annotations_filename=None):
        assert id_type in ['regions', 'genes']

        if id_type == 'regions':
            main_template_filename = TEMPLATE_WEB_REGIONS_MAIN
            candidate_targets_template_filename = TEMPLATE_WEB_REGIONS_CANDIDATE_TARGETS
            top_ranked_targets_template_filename = TEMPLATE_WEB_REGIONS_TOP_RANKED_TARGETS
            motif_scan_template_filename = TEMPLATE_WEB_REGIONS_MOTIF_SCAN
            batch_report_template_filename = TEMPLATE_WEB_REGIONS_BATCH
        else:
            raise ValueError

        gene_or_region_id2description, gene_or_region_id_locations = \
            ReportContext._read_gene_or_region_meta_data(id_description_filename,
                                                         id_location_bed_filename)
        feature_id2metadata = ReportContext._read_feature_meta_data_from_zip_archive(rankings_metadata_zip_filename,
                                                                                     motif_annotations_filename)

        if track_annotations_filename:
            feature_id2metadata.update(
                ReportContext._read_feature_track_annotations_from_file(track_annotations_filename,
                                                                        motif_annotations_filename))

        motif_collection = MotifCollection.load_from_zip_archive(rankings_metadata_zip_filename)

        return JobContext(rankings_folder,
                          db_id2filename,
                          db_id2name,
                          gene_or_region_id2description,
                          feature_id2metadata,
                          track_cluster_filename,
                          stamp_command,
                          None,
                          None,
                          results_folder,
                          webserver_url,
                          webserver_upload_folder,
                          webserver_email_address,
                          mail_server_host,
                          mail_server_port,
                          id_type,
                          gene_or_region_id_locations,
                          two_bit_to_fa_command,
                          genome_2bit_filename,
                          cluster_buster_command,
                          motif_collection,
                          main_template_filename,
                          candidate_targets_template_filename,
                          top_ranked_targets_template_filename,
                          motif_scan_template_filename,
                          db_cache_size,
                          rankings_metadata_zip_filename,
                          batch_report_template_filename)

    @staticmethod
    def create_from_files(id_type,
                          id_description_filename,
                          id_location_bed_filename,
                          logo_folder, logo_extension,
                          rankings_folder,
                          db_id2filename,
                          db_id2name,
                          results_folder,
                          webserver_url,
                          webserver_upload_folder,
                          webserver_email_address,
                          mail_server_host,
                          mail_server_port,
                          two_bit_to_fa_command,
                          genome_2bit_filename,
                          cluster_buster_command,
                          pwm_folder,
                          pwm_extension,
                          track_cluster_filename=None,
                          stamp_command=None,
                          stamp_motif_database_filename=None,
                          stamp_score_distribution_filename=None,
                          db_cache_size=0,
                          motif_annotations_filename=None,
                          track_annotations_filename=None):
        assert id_type in ['regions', 'genes']

        if id_type == 'regions':
            main_template_filename = TEMPLATE_WEB_REGIONS_MAIN
            candidate_targets_template_filename = TEMPLATE_WEB_REGIONS_CANDIDATE_TARGETS
            top_ranked_targets_template_filename = TEMPLATE_WEB_REGIONS_TOP_RANKED_TARGETS
            motif_scan_template_filename = TEMPLATE_WEB_REGIONS_MOTIF_SCAN
            batch_report_template_filename = TEMPLATE_WEB_REGIONS_BATCH
        else:
            raise ValueError

        gene_or_region_id2description, gene_or_region_id_locations = \
            ReportContext._read_gene_or_region_meta_data(id_description_filename,
                                                         id_location_bed_filename)
        feature_id2metadata = ReportContext._read_feature_meta_data_from_folder(logo_folder,
                                                                                logo_extension,
                                                                                motif_annotations_filename)

        if track_annotations_filename:
            feature_id2metadata.update(
                ReportContext._read_feature_track_annotations_from_file(track_annotations_filename,
                                                                        motif_annotations_filename))

        motif_collection = MotifCollection.load_from_folder(pwm_folder, pwm_extension)

        return JobContext(rankings_folder,
                          db_id2filename,
                          db_id2name,
                          gene_or_region_id2description,
                          feature_id2metadata,
                          track_cluster_filename,
                          stamp_command,
                          stamp_motif_database_filename,
                          stamp_score_distribution_filename,
                          results_folder,
                          webserver_url,
                          webserver_upload_folder,
                          webserver_email_address,
                          mail_server_host,
                          mail_server_port,
                          id_type,
                          gene_or_region_id_locations,
                          two_bit_to_fa_command,
                          genome_2bit_filename,
                          cluster_buster_command,
                          motif_collection,
                          main_template_filename,
                          candidate_targets_template_filename,
                          top_ranked_targets_template_filename,
                          motif_scan_template_filename,
                          db_cache_size,
                          batch_report_template_filename=batch_report_template_filename)

    def __init__(self,
                 rankings_folder,
                 db_id2filename,
                 db_id2name,
                 gene_or_region_id2description,
                 feature_id2metadata,
                 track_cluster_filename,
                 stamp_command,
                 stamp_motif_database_filename,
                 stamp_score_distribution_filename,
                 results_folder,
                 webserver_url,
                 webserver_upload_folder,
                 webserver_email_address,
                 mail_server_host,
                 mail_server_port,
                 ylabel,
                 gene_or_region_id_locations,
                 two_bit_to_fa_command,
                 genome_2bit_filename,
                 cluster_buster_command,
                 motif_collection,
                 main_template_filename,
                 candidate_targets_template_filename,
                 top_ranked_targets_template_filename,
                 motif_scan_template_filename,
                 db_cache_size,
                 rankings_metadata_zip_filename=None,
                 batch_report_template_filename=None):

        ReportContext.__init__(self,
                               ylabel,
                               gene_or_region_id2description,
                               None,
                               gene_or_region_id_locations,
                               feature_id2metadata,
                               main_template_filename,
                               candidate_targets_template_filename,
                               top_ranked_targets_template_filename,
                               track_cluster_filename,
                               stamp_command,
                               stamp_motif_database_filename,
                               stamp_score_distribution_filename,
                               rankings_metadata_zip_filename,
                               batch_report_template_filename)

        self.rankings_folder = rankings_folder
        self.db_id2filename = db_id2filename
        self.db_id2name = db_id2name
        self.results_folder = results_folder
        self.filename2description = dict(
            (os.path.join(rankings_folder, filename), self.get_description_for_db_id(db_id))
            for db_id, filename in self.db_id2filename.iteritems())
        self.webserver_url = webserver_url
        self.webserver_upload_folder = webserver_upload_folder
        self.webserver_email_address = webserver_email_address
        self.mail_server_host = mail_server_host
        self.mail_server_port = mail_server_port
        self.two_bit_to_fa_command = two_bit_to_fa_command
        self.genome_2bit_filename = genome_2bit_filename
        self.cluster_buster_command = cluster_buster_command
        self.motif_collection = motif_collection
        self.motif_scan_template_filename = motif_scan_template_filename
        self.db_cache_size = db_cache_size

    def is_database_id_defined(self, db_id):
        return db_id in self.db_id2filename

    def get_full_path_for_db_id(self, db_id):
        if not db_id:
            return None

        assert self.is_database_id_defined(db_id), "Database with ID {0:s} is unknown.".format(db_id)

        return os.path.join(self.rankings_folder, self.db_id2filename[db_id])

    def get_filename_for_db_id(self, db_id):
        return self.db_id2filename[db_id]

    def get_description_for_db_id(self, db_id):
        return self.db_id2name[db_id]

    def get_description_for_db_filename(self, db_filename):
        return self.filename2description[db_filename]
