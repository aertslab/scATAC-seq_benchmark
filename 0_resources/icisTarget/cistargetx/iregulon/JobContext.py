import os
import tempfile
import zipfile

from cistargetx.recoveryanalysis.reportcontext import ReportContext
from cistargetx.conversion.featurelist import FeatureList


def memoize(fnc):
    cache = {}

    def _helper(*args):
        if args not in cache:
            cache[args] = fnc(*args)

        return cache[args]

    return _helper


class JobContext(object):
    @staticmethod
    def create_from_cfg_file(cfg):
        motif2tf_version = cfg.get('compserver', 'motif2tf_version')
        stamp_command = os.path.expandvars(cfg.get('STAMP', 'command'))

        rankings_metadata = os.path.expandvars(cfg.get('compserver', 'rankings_metadata'))
        feature_id2metadata = ReportContext._read_feature_meta_data_from_zip_archive(rankings_metadata)

        # Loading information related to supported nomenclatures.
        nomenclature2code = dict()
        code2nomenclature_name = dict()

        for nomenclature, db_code in cfg.items('nomenclature2db_code'):
            nomenclature2code[nomenclature] = int(db_code)
            code2nomenclature_name[int(db_code)] = nomenclature

        code2nomenclature_description = dict()

        for nomenclature, name in cfg.items('nomenclature2name'):
            code2nomenclature_description[nomenclature2code[nomenclature]] = name

        # Loading information related to databases.
        rankings_folder = os.path.expandvars(cfg.get('compserver', 'rankings_folder'))
        db_name2db_filename = dict((db_name, os.path.join(rankings_folder, filename))
                                   for db_name, filename in cfg.items('databases'))
        nomenclature2db_name = dict((nomenclature, db_name)
                                    for nomenclature, db_name in cfg.items('defaultdb4nomenclature'))

        # Loading information related to gene to region conversion.
        default_overlap_fraction = cfg.getfloat('conversion_parameters', 'overlap_fraction')
        delineation_name2bed_filename = dict((name, os.path.expandvars(filename))
                                             for name, filename in cfg.items('delineation2bed'))
        nomenclature2delineation_name = dict((nomenclature_name, delineation_name)
                                             for nomenclature_name, delineation_name in cfg.items('defaultdelineation4nomenclature'))
        db_name2bed_filename = dict((name, os.path.expandvars(filename))
                                    for name, filename in cfg.items('regions4database'))

        # Loading gene annotations.
        nomenclature2annotation_filename = dict((name, os.path.expandvars(filename))
                                                for name, filename in cfg.items('genes4nomenclature'))

        return JobContext(feature_id2metadata,
                          stamp_command,
                          rankings_metadata,
                          db_name2db_filename,
                          nomenclature2db_name,
                          code2nomenclature_description,
                          default_overlap_fraction,
                          delineation_name2bed_filename,
                          nomenclature2delineation_name,
                          nomenclature2annotation_filename,
                          db_name2bed_filename,
                          code2nomenclature_name,
                          motif2tf_version)

    def __init__(self,
                 feature_id2metadata,
                 stamp_command,
                 rankings_metadata_zip_filename,
                 db_name2db_filename,
                 nomenclature2db_name,
                 code2nomenclature_description,
                 default_overlap_fraction,
                 delineation_name2bed_filename,
                 nomenclature2delineation_name,
                 nomenclature2annotation_filename,
                 db_name2bed_filename,
                 code2nomenclature_name,
                 motif2tf_version):

        self.feature_id2metadata = feature_id2metadata
        self.stamp_command = stamp_command
        self.rankings_metadata_zip_filename = rankings_metadata_zip_filename
        self.zip_input_fh = None

        self.db_name2db_filename = db_name2db_filename
        self.nomenclature2db_name = nomenclature2db_name

        self.code2nomenclature_name = code2nomenclature_name
        self.code2nomenclature_description = code2nomenclature_description

        self.default_overlap_fraction = default_overlap_fraction
        self.nomenclature2delineation_name = nomenclature2delineation_name
        self.delineation_name2bed_filename = delineation_name2bed_filename
        self.nomenclature2annotation_filename = nomenclature2annotation_filename

        self.db_name2bed_filename = db_name2bed_filename

        self.motif2tf_version = motif2tf_version

    def get_motif2tf_version(self):
        return self.motif2tf_version

    def get_default_database_name_for_nomenclature_name(self, nomenclature_name):
        return self.nomenclature2db_name[nomenclature_name]

    def get_filename_for_database_name(self, database_name):
        return self.db_name2db_filename[database_name]

    def get_name_for_nomenclature_code(self, db_code):
        assert db_code in self.code2nomenclature_description, "Unknown database code {0:d}".format(db_code)

        return self.code2nomenclature_description[db_code]

    @property
    def default_fraction_of_overlap(self):
        return self.default_overlap_fraction

    def get_delineation_for_nomenclature_name(self, nomenclature_name):
        return self.nomenclature2delineation_name[nomenclature_name]

    def is_region_based(self, database_name):
        return database_name in self.db_name2bed_filename

    @memoize
    def get_regions_for_delineation_name(self, delineation_name):
        return FeatureList.from_bed_file(self.delineation_name2bed_filename[delineation_name],
                                         transform=lambda gi: gi.split('#')[0])

    @memoize
    def get_regions_for_database_name(self, database_name):
        return FeatureList.from_bed_file(self.db_name2bed_filename[database_name])

    @memoize
    def get_gene_annotations(self, nomenclature_name):
        return FeatureList.from_bed_file(self.nomenclature2annotation_filename[nomenclature_name])

    def get_transcript_for_gene_ids(self, nomenclature_name, gene_ids):
        return self.get_gene_annotations(nomenclature_name).filter_by_name(gene_ids)

    def get_description_for_feature(self, feature_id):
        if feature_id in self.feature_id2metadata:
            return self.feature_id2metadata[feature_id].description
        else:
            return feature_id

    ####################################################################################################################
    # STAMP clustering.
    ####################################################################################################################

    @property
    def _zip_input_fh(self):
        if self.rankings_metadata_zip_filename and not self.zip_input_fh:
            self.zip_input_fh = zipfile.ZipFile(self.rankings_metadata_zip_filename, mode='r', allowZip64=True)

        return self.zip_input_fh

    def has_stamp(self):
        return self.stamp_command and (self.rankings_metadata_zip_filename
                                       or (self.stamp_motif_database_filename
                                           and self.stamp_score_distribution_filename))

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
