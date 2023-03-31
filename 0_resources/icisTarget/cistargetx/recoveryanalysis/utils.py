import functools
import itertools
import operator
import os
import re

import airspeed
from recoverycurves import RecoveryCurves


#
# Dictionary wrappers for Context instances.
#

class EmptyDictionary(object):
    def __init__(self):
        pass

    def __getitem__(self, item):
        return ""


class ContextWrapper(object):
    def __init__(self, context, attribute_name='description', with_feature_id_as_default=True):
        self.context = context
        self.attribute_name = attribute_name
        self.with_feature_id_as_default = with_feature_id_as_default

    def __getitem__(self, item):
        if item in self.context.feature_id2metadata:
            metadata = self.context.feature_id2metadata[item]

            return getattr(metadata, self.attribute_name)
        elif self.with_feature_id_as_default:
            return item
        else:
            return ''


#
# Velocity template auxiliary methods.
#

def write_template(template_filename, output_filename, namespace):
    with open(output_filename, 'w') as output_fh:
        write_template2file_handle(template_filename, output_fh, namespace)


def write_template2file_handle(template_filename, output_fh, namespace):
    template_base_folder = os.path.dirname(template_filename)

    with open(template_filename, 'r') as input_fh:
        template = airspeed.Template(input_fh.read())

    output_fh.write(str(template.merge(namespace, loader=airspeed.CachingFileLoader(template_base_folder))))


#
# Generation of overviews.
#

def _iterate_features(named_curves_iterable, remove_random_features=True):
    for name, curves in named_curves_iterable:
        iterators = map(functools.partial(RecoveryCurves.iterate_enriched_features,
                                          remove_random_features=remove_random_features),
                        curves)
        features = sorted(itertools.chain(*iterators), key=operator.itemgetter(2), reverse=True)

        for idx, (feature_id, auc, nes, curve) in enumerate(features):
            candidate_target_ids = map(operator.itemgetter(1), curve.get_candidate_target_gene_ids(feature_id))
            top_ranked_ids = map(operator.itemgetter(1), curve.get_top_ranked_gene_ids(feature_id))
            top_ranked_ranks = map(operator.itemgetter(0), curve.get_top_ranked_gene_ids(feature_id))
            database_name = curve.database.get_group_name(feature_id)

            yield (name,
                   idx + 1,
                   feature_id,
                   database_name,
                   auc,
                   nes,
                   candidate_target_ids,
                   top_ranked_ids,
                   top_ranked_ranks)


def _write_statistics_as_tsv(named_curves_iterable,
                             output_fh,
                             feature_id2description,
                             feature_id2gene_annotations,
                             remove_random_features=True):
    separator = ";"

    for (name, rank, feature_id, database_name, auc, nes, candidate_target_ids,
         top_ranked_ids, top_ranked_ranks) in _iterate_features(named_curves_iterable, remove_random_features):

        print >> output_fh, "{0:s}\t{1:d}\t{2:s}\t{3:s}\t{4:s}\t{5:s}\t{6:.10g}\t{7:.10g}\t{8:s}\t{9:s}\t{10:s}".format(
            name,
            rank,
            feature_id,
            feature_id2description[feature_id],
            ",".join(feature_id2gene_annotations[feature_id]),
            database_name,
            auc,
            nes,
            separator.join(candidate_target_ids),
            separator.join(top_ranked_ids),
            separator.join(map(str, top_ranked_ranks)))


def write_overview_as_tsv(named_curves_iterable,
                          output_fh,
                          feature_id2description=EmptyDictionary(),
                          feature_id2gene_annotations=EmptyDictionary()):
    print >> output_fh, "#GeneSignatureID\tRank\tFeatureID\tFeatureDescription\tFeatureAnnotations\tFeatureDatabase\tAUC\tNES\tCandidateTargetIDs\tTopTargetIDs\tTopTargetRanks"

    _write_statistics_as_tsv(named_curves_iterable, output_fh, feature_id2description, feature_id2gene_annotations)


def write_statistics_as_tsv(name,
                            curves,
                            output_fh,
                            feature_id2description=EmptyDictionary(),
                            feature_id2gene_annotations=EmptyDictionary()):
    print >> output_fh, "#GeneSignatureID\tRank\tFeatureID\tFeatureDescription\tFeatureAnnotations\tFeatureDatabase\tAUC\tNES\tCandidateTargetIDs\tTopTargetIDs\tTopTargetRanks"

    _write_statistics_as_tsv({name: curves}.iteritems(), output_fh, feature_id2description, feature_id2gene_annotations)


def write_full_statistics_as_tsv(name,
                                 curves,
                                 output_fh,
                                 feature_id2description=EmptyDictionary(),
                                 feature_id2gene_annotations=EmptyDictionary()):
    print >> output_fh, "#GeneSignatureID\tRank\tFeatureID\tFeatureDescription\tFeatureAnnotations\tFeatureDatabase\tAUC\tNES\tCandidateTargetIDs\tTopTargetIDs\tTopTargetRanks"

    _write_statistics_as_tsv({name: curves}.iteritems(),
                             output_fh,
                             feature_id2description,
                             feature_id2gene_annotations,
                             False)

    _write_non_enriched_features_as_tsv({name: curves}.iteritems(),
                                        output_fh,
                                        feature_id2description,
                                        feature_id2gene_annotations)


def _iterate_non_enriched_features(named_curves_iterable):
    for name, curves in named_curves_iterable:
        iterators = map(RecoveryCurves.iterate_all_features, curves)
        features = sorted(itertools.chain(*iterators), key=operator.itemgetter(2), reverse=True)

        for idx, (feature_id, auc, nes, curve) in enumerate(features):
            if nes >= curve.nes_threshold:
                continue

            database_name = curve.database.get_group_name(feature_id)

            yield name, idx + 1, feature_id, database_name, auc, nes


def _write_non_enriched_features_as_tsv(named_curves_iterable,
                                        output_fh,
                                        feature_id2description,
                                        feature_id2gene_annotations):
    for (name, rank, feature_id, database_name,
         auc, nes) in _iterate_non_enriched_features(named_curves_iterable):

        print >> output_fh, "{0:s}\t{1:d}\t{2:s}\t{3:s}\t{4:s}\t{5:s}\t{6:.10g}\t{7:.10g}\t\t\t".format(
            name,
            rank,
            feature_id,
            feature_id2description[feature_id],
            ",".join(feature_id2gene_annotations[feature_id]),
            database_name,
            auc,
            nes)


def write_overview_as_html(named_curves_iterable,
                           output_fh,
                           template_filename,
                           template_variables=dict(),
                           feature_id2description=EmptyDictionary(),
                           feature_id2gene_annotations=EmptyDictionary()):
    separator = "\t"

    class Row:
        def __init__(self, name, rank, feature_id, database_name, auc, nes, candidate_target_ids, top_ranked_ids,
                     top_ranked_ranks):
            self.name = name
            self.rank = rank
            self.feature_id = feature_id
            self.feature_description = feature_id2description[feature_id]
            self.feature_gene_annotations = ",".join(feature_id2gene_annotations[feature_id])
            self.database_name = re.sub('>', '&gt;', re.sub('<', '&lt;', database_name))
            self.auc = auc
            self.nes = format(nes, '.5f')
            self.candidate_targets = separator.join(candidate_target_ids)
            self.top_targets = separator.join(top_ranked_ids)
            self.target_ranks = separator.join(map(str, top_ranked_ranks))

    data = [Row(*row) for row in _iterate_features(named_curves_iterable)]
    namespace = dict(template_variables.iteritems())
    namespace['data'] = data

    write_template2file_handle(template_filename, output_fh, namespace)
