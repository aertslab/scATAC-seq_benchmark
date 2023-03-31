import itertools
import numpy
import operator


RANDOM_FEATURE_PREFIX = "random"
NAN = float('nan')


class RecoveryCurves:
    def __init__(self, database, auc_threshold, rank_threshold, nes_threshold):
        self.database = database
        self.auc_threshold = auc_threshold
        self.rank_threshold = rank_threshold
        self.nes_threshold = nes_threshold

        assert (self.rank_threshold > 0
                and self.rank_threshold < self.database.total_gene_count), \
            "Rank threshold must be an integer between 1 and {0:d}".format(self.database.total_gene_count)
        assert (self.auc_threshold > 0.0
                and self.auc_threshold <= 1.0), "AUC threshold must be a fraction between 0.0 and 1.0"
        assert self.rank_cutoff <= rank_threshold, \
            "An AUC threshold of {0:f} corresponds to {1:d} top ranked genes/regions in the database. " \
            "Please increase the rank threshold or decrease the AUC threshold.".format(
                self.auc_threshold,
                self.rank_cutoff
            )

        # Pre-allocation of result arrays.
        # Do not use unsigned integers.
        self.rccs = numpy.full(shape=(rank_threshold, database.feature_count),
                               fill_value=0,
                               dtype=numpy.int_)
        self.aucs = numpy.full(shape=database.feature_count,
                               fill_value=0,
                               dtype=numpy.float64)
        self.normalized_enrichment_scores = numpy.empty(shape=database.feature_count,
                                                        dtype=numpy.float64)
        # Not really pre-allocation.
        self.avgrcc = numpy.empty(shape=(rank_threshold, 0))
        self.avg2stdrcc = numpy.empty(shape=(rank_threshold, 0))

        rank_cutoff = self.rank_cutoff
        maxauc = float(rank_cutoff * database.gene_count)

        for idx in numpy.arange(database.feature_count):
            # TODO: weighted recovery curve can be speeded up via the following numpy procedure:
            #         numpy.histogram(a, bins=10, range=None, normed=False, weights=None, density=None)
            #         numpy.histogram(database.rankings[:, idx], bins=rank_cutoff, range=(0, rank_cutoff), weights=weights)
            curranking = numpy.append(database.rankings[:, idx], database.total_gene_count)
            self.rccs[:, idx] = numpy.cumsum(numpy.bincount(curranking)[:rank_threshold])
            self.aucs[idx] = numpy.sum(self.rccs[:rank_cutoff, idx]) / maxauc

        self.avgrcc = numpy.average(self.rccs, axis=1)
        self.stdrcc = numpy.std(self.rccs, axis=1)
        self.avg2stdrcc = self.avgrcc + 2.0 * self.stdrcc
        self.normalized_enrichment_scores = (self.aucs - self.aucs.mean()) / self.aucs.std()
        # logging.debug("AUCs for '{0:s}': ".format(database.name) + " ".join(map(str, self.aucs)))

    @property
    def gene_ids(self):
        return self.database.gene_ids

    @property
    def rank_cutoff(self):
        return int(round(self.auc_threshold * self.database.total_gene_count))

    @property
    def enriched_feature_count(self):
        return (self.normalized_enrichment_scores >= self.nes_threshold).sum()

    @property
    def enriched_feature_count_without_randoms(self):
        rsum = 0

        for feature_id, auc, nes in self.enriched_features:
            if feature_id.startswith(RANDOM_FEATURE_PREFIX):
                break

            rsum += 1

        return rsum

    @property
    def enriched_features(self):
        return sorted(((feature_id, auc, nes)
                       for feature_id, auc, nes in zip(self.database.feature_ids, self.aucs, self.normalized_enrichment_scores)
                       if nes > self.nes_threshold),
                      reverse=True,
                      key=operator.itemgetter(2))

    @property
    def all_features(self):
        return sorted(zip(self.database.feature_ids, self.aucs, self.normalized_enrichment_scores),
                      reverse=True,
                      key=operator.itemgetter(2))

    def get_recovery_curve(self, feature_id):
        idx = list(self.database.feature_ids).index(feature_id)

        return self.rccs[:, idx]

    def get_critical_point(self, feature_id):
        """ Returns (x,y). """
        rcc = self.get_recovery_curve(feature_id)
        x_values = numpy.arange(1, self.rank_threshold + 1)
        y_values = rcc - self.avg2stdrcc
        y_max = y_values.max()
        x_max = int(x_values[y_values == y_max][0])

        return x_max, rcc[x_max - 1]

    def get_ppv_tpr(self, feature_id, true_targets):
        """ (PPV, TPR, TP, FN, FP) """
        predicted_targets_set = set(map(operator.itemgetter(1), self.get_candidate_target_gene_ids(feature_id)))
        true_targets_set = set(true_targets)
        true_target_count = len(true_targets_set & predicted_targets_set)

        return ((NAN
                 if len(predicted_targets_set) == 0
                 else float(true_target_count) / float(len(predicted_targets_set))),
                (NAN
                 if len(true_targets_set) == 0
                 else float(true_target_count) / float(len(true_targets_set))),
                true_target_count,
                len(true_targets_set) - true_target_count,
                len(predicted_targets_set) - true_target_count)

    def iterate_enriched_features(self, true_targets=None, remove_random_features=True):
        for feature_id, auc, nes in self.enriched_features:
            if remove_random_features and feature_id.startswith(RANDOM_FEATURE_PREFIX):
                break

            columns = [feature_id, auc, nes, self]

            if true_targets:
                ppv, tpr, tp, fn, fp = self.get_ppv_tpr(feature_id, true_targets)
                columns += [ppv, tpr, tp, fn, fp]

            yield columns

    @property
    def has_enriched_random_features(self):
        return any(feature_id.startswith(RANDOM_FEATURE_PREFIX) for feature_id, auc, nes in self.enriched_features)

    @property
    def adjusted_enrichement_threshold(self):
        try:
            return max(nes
                       for feature_id, auc, nes in self.enriched_features
                       if feature_id.startswith(RANDOM_FEATURE_PREFIX))
        except ValueError:
            # ValueError: max() arg is an empty sequence.
            return self.nes_threshold

    @property
    def has_no_enriched_features(self):
        return self.enriched_feature_count_without_randoms == 0

    def iterate_all_features(self):
        for feature_id, auc, nes in self.all_features:
            yield feature_id, auc, nes, self

    def fraction_of_mapped_gene_ids(self, gene_ids):
        gene_ids_as_set = set(gene_ids)
        total_count = len(gene_ids_as_set)

        if total_count == 0:
            return 0.0
        else:
            return float(len(set(self.gene_ids) & gene_ids_as_set)) / total_count

    @property
    def has_top_regions(self):
        return self.database.has_top_regions

    def get_top_regions(self, rank_threshold, feature_id):
        return self.database.get_top_regions(rank_threshold, feature_id)

    def _get_gene_ids(self, feature_id, rank_threshold):
        gene_ids = self.database.gene_ids
        ranking = self.database.get_ranking(feature_id)
        sorted_idx = numpy.argsort(ranking)
        ranking = ranking[sorted_idx]
        gene_ids = gene_ids[sorted_idx]
        filtered_idx = ranking < rank_threshold

        return zip(ranking[filtered_idx] + 1, gene_ids[filtered_idx])

    def get_candidate_target_gene_ids(self, feature_id):
        return self._get_gene_ids(feature_id, self.get_critical_point(feature_id)[0])

    def get_top_ranked_gene_ids(self, feature_id):
        return self._get_gene_ids(feature_id, self.rank_threshold)

    def get_all_candidate_target_gene_ids(self, feature_id):
        rank_threshold = self.get_critical_point(feature_id)[0]
        gene_id_set = set(self.database.gene_ids)

        return [(rank, gene_id, gene_id in gene_id_set)
                for rank, gene_id in self.get_top_regions(rank_threshold, feature_id)]

    def get_all_top_ranked_gene_ids(self, feature_id):
        rank_threshold = self.rank_threshold
        gene_id_set = set(self.database.gene_ids)

        return [(rank, gene_id, gene_id in gene_id_set)
                for rank, gene_id in self.get_top_regions(rank_threshold, feature_id)]
