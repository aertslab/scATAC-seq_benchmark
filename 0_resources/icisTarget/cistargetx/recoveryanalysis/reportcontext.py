import collections
import os
import shutil
import tempfile
import zipfile
from collections import defaultdict

from cistargetx.conversion.featurelist import FeatureList


TEMPLATE_CLI_REGIONS_MAIN = os.path.join(os.path.dirname(__file__), 'cli-regions-main.vm')
TEMPLATE_CLI_GENES_MAIN = os.path.join(os.path.dirname(__file__), 'cli-genes-main.vm')
TEMPLATE_CLI_REGIONS_CANDIDATE_TARGETS = os.path.join(os.path.dirname(__file__), 'cli-regions-candidate-targets.vm')

TEMPLATE_CLI_REGIONS_TOP_RANKED_TARGETS = os.path.join(os.path.dirname(__file__), 'cli-regions-top-ranked-targets.vm')
TEMPLATE_CLI_GENES_CANDIDATE_TARGETS = os.path.join(os.path.dirname(__file__), 'cli-genes-candidate-targets.vm')
TEMPLATE_CLI_GENES_TOP_RANKED_TARGETS = os.path.join(os.path.dirname(__file__), 'cli-genes-top-ranked-targets.vm')

TEMPLATE_CLI_BATCH_REPORT = os.path.join(os.path.dirname(__file__), 'cli-batch-report.vm')


def load_gene_ids(filename):
    with open(filename, 'r') as input_fh:
        return set(line.rstrip() for line in input_fh.readlines())


def load_dictionary(filename):
    lut = dict()

    with open(filename, 'r') as input_fh:
        for line in input_fh:
            if line.startswith('#'):
                continue

            columns = line.rstrip('\n').split('\t')

            for elem in columns[1:]:
                lut.setdefault(columns[0], set()).add(elem)

    return lut


FeatureMetaData = collections.namedtuple('FeatureMetaData',
                                         'source_id source_version description logo_filename cb_filename gene_ids')


class ReportContext:
    @staticmethod
    def _iterate_motif_annotation_table(filename):
        with open(filename, 'r') as input_fh:
            for line in input_fh:
                if line.startswith('#') or not line.strip():
                    continue

                columns = line.rstrip('\n').split('\t')

                if len(columns) == 4:
                    # motif_id  motif_description  gene_name  description
                    yield columns[0], columns[2]
                if len(columns) == 13:
                    # motif_id  motif_name  motif_description  source_name  source_version  gene_name
                    # motif_similarity_qvalue  similar_motif_id  similar_motif_description
                    # orthologous_identity  orthologous_gene_name  orthologous_species  description
                    yield columns[0], columns[5]

    @staticmethod
    def _load_motif_annotation_table(filename):
        feature_id2gene_ids = defaultdict(set)

        for motif_id, gene_id in ReportContext._iterate_motif_annotation_table(filename):
            feature_id2gene_ids[motif_id].add(gene_id)

        return feature_id2gene_ids

    @staticmethod
    def _read_feature_meta_data_from_zip_archive(rankings_metadata_zip_filename, gene_id_annotation=None):
        feature_id2gene_ids = (ReportContext._load_motif_annotation_table(gene_id_annotation)
                               if gene_id_annotation
                               else dict())
        feature_id2metadata = dict()

        with zipfile.ZipFile(rankings_metadata_zip_filename, mode='r', allowZip64=True) as zip_input_fh:
            with zip_input_fh.open('motifs.tsv', 'r') as input_fh:
                for line in input_fh:
                    columns = line.rstrip('\n').split('\t')
                    feature_id = columns[3]
                    logo_filename = 'logos/{0:s}.png'.format(feature_id)
                    cb_filename = 'singletons/{0:s}.cb'.format(feature_id)

                    feature_id2metadata[feature_id] = FeatureMetaData(columns[0],
                                                                      columns[1],
                                                                      columns[4],
                                                                      logo_filename,
                                                                      cb_filename,
                                                                      feature_id2gene_ids.get(feature_id, set()))

        return feature_id2metadata

    @staticmethod
    def _read_feature_meta_data_from_folder(logo_folder, logo_extension, gene_id_annotation=None):
        feature_id2gene_ids = (ReportContext._load_motif_annotation_table(gene_id_annotation)
                               if gene_id_annotation
                               else dict())
        feature_id2metadata = dict()

        for feature_id, logo_filename in ((filename[:-(1 + len(logo_extension))], os.path.join(logo_folder, filename))
                                          for filename in os.listdir(logo_folder)
                                          if filename.endswith(logo_extension)):

            feature_id2metadata[feature_id] = FeatureMetaData('?',
                                                              '?',
                                                              feature_id,
                                                              logo_filename,
                                                              '',
                                                              feature_id2gene_ids.get(feature_id, set()))

        return feature_id2metadata

    @staticmethod
    def _iterate_track_annotation_table(filename):
        # source_id source_version track_id track_description
        with open(filename, 'r') as input_fh:
            for line in input_fh:
                if line.startswith('#') or not line.strip():
                    continue

                columns = line.rstrip('\n').split('\t')

                yield columns[0], columns[1], columns[2], columns[3]

    @staticmethod
    def _read_feature_track_annotations_from_file(track_annotations_filename, gene_id_annotation=None):
        feature_id2gene_ids = (ReportContext._load_motif_annotation_table(gene_id_annotation)
                               if gene_id_annotation
                               else dict())
        feature_id2metadata = dict()

        for source_id, source_version, feature_id, description in ReportContext._iterate_track_annotation_table(
                track_annotations_filename):

            feature_id2metadata[feature_id] = FeatureMetaData(source_id,
                                                              source_version,
                                                              description if description else feature_id,
                                                              '',
                                                              '',
                                                              feature_id2gene_ids.get(feature_id, set()))

        return feature_id2metadata

    @staticmethod
    def _read_gene_or_region_meta_data(id_description_filename,
                                       id_location_bed_filename=None,
                                       id_description_highlighting_filename=None):
        gene_or_region_id2description = (load_dictionary(os.path.expandvars(id_description_filename))
                                         if id_description_filename
                                         else dict())
        gene_or_region_id_locations = (FeatureList.from_bed_file(os.path.expandvars(id_location_bed_filename))
                                       if id_location_bed_filename
                                       else None)

        if id_description_highlighting_filename:
            elements_to_highlight = (load_gene_ids(os.path.expandvars(id_description_highlighting_filename))
                                     if id_description_highlighting_filename
                                     else set())

            return gene_or_region_id2description, gene_or_region_id_locations, elements_to_highlight
        else:
            return gene_or_region_id2description, gene_or_region_id_locations

    @staticmethod
    def create_from_zip_archive(id_type,
                                id_description_filename,
                                id_description_highlighting_filename,
                                id_location_bed_filename,
                                rankings_metadata_zip_filename,
                                track_cluster_filename=None,
                                stamp_command=None,
                                motif_annotations_filename=None,
                                track_annotations_filename=None):
        assert id_type in ['regions', 'genes']

        if id_type == 'regions':
            main_template_filename = TEMPLATE_CLI_REGIONS_MAIN
            candidate_targets_template_filename = TEMPLATE_CLI_REGIONS_CANDIDATE_TARGETS
            top_ranked_targets_template_filename = TEMPLATE_CLI_REGIONS_TOP_RANKED_TARGETS
        else:
            main_template_filename = TEMPLATE_CLI_GENES_MAIN
            candidate_targets_template_filename = TEMPLATE_CLI_GENES_CANDIDATE_TARGETS
            top_ranked_targets_template_filename = TEMPLATE_CLI_GENES_TOP_RANKED_TARGETS

        batch_report_template_filename = TEMPLATE_CLI_BATCH_REPORT

        if id_description_highlighting_filename:
            gene_or_region_id2description, gene_or_region_id_locations, elements_to_highlight = \
                ReportContext._read_gene_or_region_meta_data(id_description_filename,
                                                             id_location_bed_filename,
                                                             id_description_highlighting_filename)
        else:
            gene_or_region_id2description, gene_or_region_id_locations = \
                ReportContext._read_gene_or_region_meta_data(id_description_filename,
                                                             id_location_bed_filename)
            elements_to_highlight = None

        feature_id2metadata = ReportContext._read_feature_meta_data_from_zip_archive(rankings_metadata_zip_filename,
                                                                                     motif_annotations_filename)
        if track_annotations_filename:
            feature_id2metadata.update(
                ReportContext._read_feature_track_annotations_from_file(track_annotations_filename,
                                                                        motif_annotations_filename))

        return ReportContext(id_type,
                             gene_or_region_id2description,
                             elements_to_highlight,
                             gene_or_region_id_locations,
                             feature_id2metadata,
                             main_template_filename,
                             candidate_targets_template_filename,
                             top_ranked_targets_template_filename,
                             track_cluster_filename,
                             stamp_command,
                             None,
                             None,
                             rankings_metadata_zip_filename,
                             batch_report_template_filename)

    @staticmethod
    def create_from_files(id_type,
                          id_description_filename,
                          id_description_highlighting_filename,
                          id_location_bed_filename,
                          logo_folder, logo_extension='png',
                          track_cluster_filename=None,
                          stamp_command=None,
                          stamp_motif_database_filename=None,
                          stamp_score_distribution_filename=None,
                          motif_annotations_filename=None,
                          track_annotations_filename=None):
        assert id_type in ['regions', 'genes']

        if id_type == 'regions':
            main_template_filename = TEMPLATE_CLI_REGIONS_MAIN
            candidate_targets_template_filename = TEMPLATE_CLI_REGIONS_CANDIDATE_TARGETS
            top_ranked_targets_template_filename = TEMPLATE_CLI_REGIONS_TOP_RANKED_TARGETS
        else:
            main_template_filename = TEMPLATE_CLI_GENES_MAIN
            candidate_targets_template_filename = TEMPLATE_CLI_GENES_CANDIDATE_TARGETS
            top_ranked_targets_template_filename = TEMPLATE_CLI_GENES_TOP_RANKED_TARGETS

        batch_report_template_filename = TEMPLATE_CLI_BATCH_REPORT

        if id_description_highlighting_filename:
            gene_or_region_id2description, gene_or_region_id_locations, elements_to_highlight = \
                ReportContext._read_gene_or_region_meta_data(id_description_filename,
                                                             id_location_bed_filename,
                                                             id_description_highlighting_filename)
        else:
            gene_or_region_id2description, gene_or_region_id_locations = \
                ReportContext._read_gene_or_region_meta_data(id_description_filename,
                                                             id_location_bed_filename)
            elements_to_highlight = None

        feature_id2metadata = dict()

        if logo_folder and os.path.exists(logo_folder):
            feature_id2metadata = ReportContext._read_feature_meta_data_from_folder(logo_folder,
                                                                                    logo_extension,
                                                                                    motif_annotations_filename)

            if track_annotations_filename and os.path.isfile(track_annotations_filename):
                feature_id2metadata.update(
                    ReportContext._read_feature_track_annotations_from_file(track_annotations_filename,
                                                                            motif_annotations_filename))
        elif track_annotations_filename and os.path.isfile(track_annotations_filename):
            feature_id2metadata = ReportContext._read_feature_track_annotations_from_file(track_annotations_filename,
                                                                                          motif_annotations_filename)

        return ReportContext(id_type,
                             gene_or_region_id2description,
                             elements_to_highlight,
                             gene_or_region_id_locations,
                             feature_id2metadata,
                             main_template_filename,
                             candidate_targets_template_filename,
                             top_ranked_targets_template_filename,
                             track_cluster_filename,
                             stamp_command,
                             stamp_motif_database_filename,
                             stamp_score_distribution_filename,
                             None,
                             batch_report_template_filename)

    def __init__(self,
                 ylabel='regions',
                 gene_or_region_id2description=None,
                 elements_to_highlight=None,
                 gene_or_region_id_locations=None,
                 feature_id2metadata=dict(),
                 main_template_filename=TEMPLATE_CLI_REGIONS_MAIN,
                 candidate_targets_template_filename=TEMPLATE_CLI_REGIONS_CANDIDATE_TARGETS,
                 top_ranked_targets_template_filename=TEMPLATE_CLI_REGIONS_TOP_RANKED_TARGETS,
                 track_cluster_filename=None,
                 stamp_command=None,
                 stamp_motif_database_filename=None,
                 stamp_score_distribution_filename=None,
                 rankings_metadata_zip_filename=None,
                 batch_report_template_filename=None):
        self.ylabel = ylabel

        self.gene_or_region_id2description = gene_or_region_id2description
        self.elements_to_highlight = elements_to_highlight
        self.gene_or_region_id_locations = gene_or_region_id_locations
        self.feature_id2metadata = feature_id2metadata

        self.main_template_filename = main_template_filename
        self.candidate_targets_template_filename = candidate_targets_template_filename
        self.top_ranked_targets_template_filename = top_ranked_targets_template_filename
        self.batch_report_template_filename = batch_report_template_filename
        self.track_cluster_filename = track_cluster_filename

        self.stamp_command = os.path.expandvars(stamp_command) if stamp_command else None
        self.stamp_motif_database_filename = (os.path.expandvars(stamp_motif_database_filename)
                                              if stamp_motif_database_filename
                                              else None)
        self.stamp_score_distribution_filename = (os.path.expandvars(stamp_score_distribution_filename)
                                                  if stamp_score_distribution_filename
                                                  else None)

        self.rankings_metadata_zip_filename = rankings_metadata_zip_filename
        self.zip_input_fh = None

    @property
    def _zip_input_fh(self):
        if self.rankings_metadata_zip_filename and not self.zip_input_fh:
            self.zip_input_fh = zipfile.ZipFile(self.rankings_metadata_zip_filename, mode='r', allowZip64=True)

        return self.zip_input_fh

    def _has_cb_filename_for_feature(self, feature_id):
        if feature_id in self.feature_id2metadata:
            return True if self.feature_id2metadata[feature_id].cb_filename else False
        else:
            return False

    def _get_cb_content_for_feature(self, feature_id):
        if feature_id in self.feature_id2metadata:
            if self.rankings_metadata_zip_filename:
                with self._zip_input_fh.open(self.feature_id2metadata[feature_id].cb_filename, 'r') as input_fh:
                    return input_fh.read()
            else:
                return ''
        else:
            return ''

    def write_clusterbuster_motifs(self, feature_ids, cb_motifs_output_filename, include_transfac_pro=False):
        with open(cb_motifs_output_filename, 'w') as fh_w:
            for feature_id in feature_ids:
                if not include_transfac_pro and (
                            feature_id.startswith('transfac_pro') or feature_id.startswith('transfac__pro')):
                    continue
                if self._has_cb_filename_for_feature(feature_id):
                    fh_w.write(self._get_cb_content_for_feature(feature_id))

    def _has_logo_for_feature(self, feature_id):
        if feature_id in self.feature_id2metadata:
            return True if self.feature_id2metadata[feature_id].logo_filename else False
        else:
            return False

    def _get_logo_filename_for_feature(self, feature_id):
        if feature_id in self.feature_id2metadata:
            return self.feature_id2metadata[feature_id].logo_filename
        else:
            return None

    def copy_logo_for_feature_to(self, feature_id, destination):
        if self._has_logo_for_feature(feature_id):
            original_logo_filename = self._get_logo_filename_for_feature(feature_id)
            basename, extension = os.path.splitext(os.path.basename(original_logo_filename))
            logo_filename = basename + ".logo" + extension

            if self.rankings_metadata_zip_filename:
                with open(os.path.join(destination, logo_filename), 'wb') as output_fh:
                    with self._zip_input_fh.open(original_logo_filename, 'r') as input_fh:
                        output_fh.write(input_fh.read())
            else:
                shutil.copy(original_logo_filename,
                            os.path.join(destination, logo_filename))

            return logo_filename
        else:
            return None

    def get_description_for_feature(self, feature_id):
        if feature_id in self.feature_id2metadata:
            return self.feature_id2metadata[feature_id].description
        else:
            return feature_id

    def get_gene_annotations_for_feature(self, feature_id):
        if feature_id in self.feature_id2metadata:
            return self.feature_id2metadata[feature_id].gene_ids
        else:
            return set()

    def has_stamp(self):
        return (self.stamp_command
                and (self.rankings_metadata_zip_filename
                     or (self.stamp_motif_database_filename and self.stamp_score_distribution_filename)))

    def has_cluster(self):
        return self.track_cluster_filename

    def get_track_cluster_filename(self):
        return self.track_cluster_filename

    def get_stamp_score_distribution_filename(self):
        if self.rankings_metadata_zip_filename:
            self.stamp_score_distribution_filename = tempfile.mktemp("stamp.scores")

            with open(self.stamp_score_distribution_filename, 'w') as output_fh:
                with self._zip_input_fh.open('motifs.scores', 'r') as input_fh:
                    output_fh.write(input_fh.read())

        return self.stamp_score_distribution_filename

    def get_stamp_motif_database_filename(self):
        if self.rankings_metadata_zip_filename:
            self.stamp_motif_database_filename = tempfile.mktemp("stamp.transfac")

            with open(self.stamp_motif_database_filename, 'w') as output_fh:
                with self._zip_input_fh.open('motifs.transfac', 'r') as input_fh:
                    output_fh.write(input_fh.read())

        return self.stamp_motif_database_filename

    def remove_stamp_intermediate_files(self):
        if self.rankings_metadata_zip_filename:
            os.remove(self.stamp_motif_database_filename)
            os.remove(self.stamp_score_distribution_filename)
            self.stamp_motif_database_filename = None
            self.stamp_score_distribution_filename = None

    def get_ranked_gene_or_region_ids_for_template(self, rank_target_gene_or_region_ids, elements_to_highlight=None):
        if not elements_to_highlight:
            elements_to_highlight = self.elements_to_highlight

        def highlight_if_necessary(element):
            if not elements_to_highlight or element not in elements_to_highlight:
                return element
            else:
                return '<b>{0:s}</b>'.format(element)

        results = []

        for rank, gene_or_region_id in rank_target_gene_or_region_ids:
            results.append({
                'rank': rank,
                'id': gene_or_region_id,
                'description': ' '.join(map(highlight_if_necessary,
                                            self.gene_or_region_id2description.get(gene_or_region_id, '')))
            })

        results.sort(key=lambda e: e['rank'])

        return results

    def get_ranked_gene_or_region_ids_for_tsv(self, rank_target_gene_or_region_ids):
        results = []

        for rank, gene_or_region_id in rank_target_gene_or_region_ids:
            results.append({
                'rank': rank,
                'id': gene_or_region_id,
                'description': ' '.join(self.gene_or_region_id2description.get(gene_or_region_id, ''))
            })

        results.sort(key=lambda e: e['rank'])

        return results

    def get_locations_for_template(self, rank_target_gene_or_region_ids):
        if not self.gene_or_region_id_locations:
            return None

        results = []

        for rank, gene_or_region_id in rank_target_gene_or_region_ids:
            for feature in self.gene_or_region_id_locations.filter_by_name([gene_or_region_id]):
                results.append({
                    'chromosome': feature.chromosome,
                    'start': feature.interval[0],
                    'end': feature.interval[1],
                    'id': gene_or_region_id
                })

        return results
