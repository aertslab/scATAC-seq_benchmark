import logging
import operator
from collections import defaultdict
from textwrap import wrap

import MySQLdb

from cistargetx.common.rankingsdatabase import RankingsDatabase
from cistargetx.conversion.featurelist import Feature
from cistargetx.recoveryanalysis.motifclustering import cluster
from cistargetx.recoveryanalysis.recoverycurves import RecoveryCurves


START_JOB_STATEMENT = r"""
  UPDATE jobQueue
  SET jobStatusCode = 2, jobStartTime = now()
  WHERE ID = %s AND jobStatusCode = 1;
"""

STOP_JOB_STATEMENT = r"""
  UPDATE jobQueue
  SET jobStatusCode = 3, jobFinishTime = now()
  WHERE ID = %s;
"""

ERROR_JOB_STATEMENT = r"""
  UPDATE jobQueue
  SET jobStatusCode = -1, jobFinishTime = now(), errorMessage = %s
  WHERE ID = %s;
"""

INSERT_ENRICHED_FEATURE_STATEMENT = r"""
  INSERT INTO enrichedFeature(jobID, rankingsdb, rank, name, description, AUCValue, NEScore, clusterNumber,
                              candidateTargetIDs, candidateTargetRanks)
  VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
"""

INSERT_CANDIDATE_ENHANCER_STATEMENT = r"""
  INSERT INTO candidateEnhancer(featureID, rank, targetGeneName, chromosome, startInclusive, endExclusive, regionName)
  VALUES (%s, %s, %s, %s, %s, %s, %s);
"""

INSERT_FEATURE_ANNOTATION_STATEMENT = r"""
  INSERT INTO featureAnnotation(featureID, transcriptionFactorName, orthologousIdentity, orthologousGeneName,
                                orthologousSpecies, motifSimilarityFDR, similarMotifName, similarMotifDescription)
  VALUES (%s, %s, %s, %s, %s, %s, %s, %s);
"""

QUERY_MOTIF_ANNOTATION_STATEMENT = r"""
  SELECT geneName, orthologousIdentity, orthologousGeneName, orthologousSpecies,
         motifSimilarityQValue, similarMotifName, similarMotifDescription
  FROM motifFeature_{0:s}
  WHERE motifName = %s
    AND sourceName = %s
    AND (similarMotifName IS NULL OR motifSimilarityQValue <= %s)
    AND (orthologousGeneName IS NULL OR orthologousIdentity >= %s)
    AND nomenclatureCode = %s;
"""

FRACTION_OF_MAPPED_GENE_IDS = 0.60


class JobRunningInOtherInstanceError(Exception):
    pass


class UnsupportedRankingsDatabaseError(Exception):
    pass


class Job(object):
    def __init__(self,
                 connection,
                 configuration,
                 context,
                 job_id,
                 name,
                 nomenclature_code,
                 ranking_database_codes,
                 conversion_delineation,
                 conversion_upstream_region_in_bp,
                 conversion_downstream_region_in_bp,
                 conversion_fraction_of_overlap,
                 rank_threshold,
                 auc_rank_threshold_as_percentage,
                 nes_threshold,
                 min_orthologous_identity,
                 max_motif_similarity_fdr,
                 gene_ids,
                 ip,
                 user_agent):

        self.connection = connection
        self.context = context
        self.configuration = configuration

        self.job_id = job_id
        self.name = unicode(name, errors='replace').encode('ascii', 'replace')

        self.nomenclature_code = nomenclature_code
        self.nomenclature_name = self.context.code2nomenclature_name[nomenclature_code]

        self.ranking_db_codes = [None if ranking_db_code == '' else ranking_db_code
                                 for ranking_db_code in ranking_database_codes.split(';')]

        self.conversion_delineation = conversion_delineation if conversion_delineation else None
        self.conversion_upstream_region_in_bp = (int(conversion_upstream_region_in_bp)
                                                 if conversion_upstream_region_in_bp
                                                 else None)
        self.conversion_downstream_region_in_bp = (int(conversion_downstream_region_in_bp)
                                                   if conversion_downstream_region_in_bp
                                                   else None)
        self.conversion_fraction_of_overlap = conversion_fraction_of_overlap if conversion_fraction_of_overlap else None

        self.rank_threshold = rank_threshold
        self.auc_threshold = auc_rank_threshold_as_percentage
        self.nes_threshold = nes_threshold

        self.min_orthologous_identity = min_orthologous_identity
        self.max_motif_similarity_fdr = max_motif_similarity_fdr

        self.gene_ids = set(map(lambda s: s.strip(), gene_ids.split(';')))

        self.ip = ip
        self.user_agent = user_agent

        self.motif2tf_version = self.context.get_motif2tf_version()


    ####################################################################################################################
    # Basic job operations ...
    ####################################################################################################################

    def reset(self):
        self.connection = MySQLdb.connect(host=self.configuration.servername,
                                          port=self.configuration.port,
                                          user=self.configuration.username,
                                          passwd=self.configuration.password,
                                          db=self.configuration.database)

        # Automatically reconnect to the MySQL server if the connection goes away.
        self.connection.ping(True)

    def start(self):
        # Only handle the job if all rankings databases which are requested are specified in the config file.
        for ranking_db in self.ranking_db_codes:
            if ranking_db not in self.context.db_name2db_filename:
                raise UnsupportedRankingsDatabaseError

        cursor = self.connection.cursor()
        cursor.execute(START_JOB_STATEMENT, (self.job_id,))
        self.connection.commit()
        job_already_running = cursor.rowcount == 0
        cursor.close()

        if job_already_running:
            # Skip this job if it was already started in the meanwhile by another worker thread.
            raise JobRunningInOtherInstanceError

        # Check sanity of input parameters.
        assert self.name.strip(), "A non-empty job name must be supplied."
        assert self.auc_threshold > 0.0 and self.auc_threshold <= 1.0, "AUC threshold must be a fraction between 0.0 and 1.0 and cannot be zero."
        assert self.min_orthologous_identity >= 0.0 and self.min_orthologous_identity <= 1.0, "Minimum orthologous identity should be fraction between 0.0 and 1.0."
        assert self.max_motif_similarity_fdr >= 0.0 and self.max_motif_similarity_fdr <= 1.0, "Maximum motif similarity FDR should be fraction between 0.0 and 1.0."
        assert len(self.gene_ids) > 1, "Number of supplied genes must be greater than 1."

    def finish(self):
        cursor = self.connection.cursor()
        cursor.execute(STOP_JOB_STATEMENT, (self.job_id,))
        self.connection.commit()
        cursor.close()

    def error(self, message):
        cursor = self.connection.cursor()
        cursor.execute(ERROR_JOB_STATEMENT, (message, self.job_id))
        self.connection.commit()
        cursor.close()

    ####################################################################################################################
    # Parameters.
    ####################################################################################################################

    @property
    def database_names(self):
        return [self.context.get_default_database_name_for_nomenclature_name(self.nomenclature_name)
                if ranking_db_code is None
                else ranking_db_code
                for ranking_db_code in self.ranking_db_codes]

    def database_filename(self, database_name):
        return self.context.get_filename_for_database_name(database_name)

    def region_based_database(self, database_name):
        return self.context.is_region_based(database_name)

    @property
    def fraction_of_overlap(self):
        if self.conversion_fraction_of_overlap:
            return self.conversion_fraction_of_overlap
        else:
            return self.context.default_fraction_of_overlap

    def conversion_delineation_name(self, database_name):
        assert self.region_based_database(database_name)

        return (self.conversion_delineation
                if self.conversion_delineation
                else self.context.get_delineation_for_nomenclature_name(self.nomenclature_name))

    def conversion_interval(self, database_name):
        assert self.region_based_database(database_name)

        if self.conversion_upstream_region_in_bp and self.conversion_downstream_region_in_bp:
            return self.conversion_upstream_region_in_bp, self.conversion_downstream_region_in_bp
        else:
            return ()

    def conversion_intervals_iterator(self, database_name):
        assert self.conversion_interval(database_name)

        # Retrieve all transcripts for the genes.
        transcripts = self.context.get_transcript_for_gene_ids(self.nomenclature_name, self.gene_ids)

        # Create putative regulatory region features.
        upstream, downstream = self.conversion_interval(database_name)

        for transcript in transcripts:
            name = transcript.name
            strand = transcript.strand
            chromosome = transcript.chromosome

            if strand == "+":
                tss = transcript.interval[0]

                yield Feature(chromosome, max(0, tss - upstream), tss + downstream, name=name, strand=strand)
            elif strand == "-":
                tss = transcript.interval[1] - 1

                yield Feature(chromosome, max(0, tss - downstream), tss + upstream, name=name, strand=strand)
            else:
                raise ValueError('Invalid strand "{0:s}" in gene annotation file'.format(strand))

    def gene_ids2regions(self, gene_ids, database_name):
        scored_regions = self.context.get_regions_for_database_name(database_name)
        regions = set()
        region_id2_gene_ids = defaultdict(set)

        if self.conversion_interval(database_name):
            for feature in self.conversion_intervals_iterator(database_name):
                for region in scored_regions.find_overlap_with_feature(feature, self.fraction_of_overlap):
                    regions.add(region.name)
                    region_id2_gene_ids[region.name].add(feature.name)
        else:
            putative_regulatory_regions = self.context.get_regions_for_delineation_name(
                self.conversion_delineation_name(database_name))

            for gene_id in gene_ids:
                for putative_regulatory_region in putative_regulatory_regions.filter_by_name([gene_id]):
                    for scored_region in scored_regions.find_overlap_with_feature(putative_regulatory_region,
                                                                                  self.fraction_of_overlap):
                        regions.add(scored_region.name)
                        region_id2_gene_ids[scored_region.name].add(gene_id)

        return regions, region_id2_gene_ids

    def fraction_of_mapped_gene_ids(self, database_name):
        assert self.region_based_database(database_name)

        putative_regulatory_regions = self.context.get_regions_for_delineation_name(
            self.conversion_delineation_name(database_name))
        all_gene_ids = set(putative_regulatory_regions.feature_ids)

        return float(len(all_gene_ids & self.gene_ids)) / len(self.gene_ids)

    def load_database(self, database_name):
        if self.region_based_database(database_name):
            regions_ids, region_id2_gene_ids = self.gene_ids2regions(self.gene_ids, database_name)

            return (RankingsDatabase.load(self.database_filename(database_name),
                                          database_name,
                                          gene_ids=regions_ids),
                    region_id2_gene_ids)
        else:
            return (RankingsDatabase.load(self.database_filename(database_name),
                                          database_name,
                                          gene_ids=self.gene_ids),
                    dict(zip(self.gene_ids, map(lambda e: [e], self.gene_ids))))

    ####################################################################################################################
    # Database operations.
    ####################################################################################################################

    def _iterate_tfs(self, feature_id):
        cursor = self.connection.cursor()

        columns_double_underscore_splitted = feature_id.split('__', 1)
        columns_dash_splitted = feature_id.split('-', 1)

        # Check if we have a feature ID in the form of "motifcollection__motifid" or "motifcollection-motifid".
        if len(columns_double_underscore_splitted) == 2:
            source_name, source_id = columns_double_underscore_splitted

            if len(columns_dash_splitted) == 2:
                if len(columns_double_underscore_splitted[0]) > len(columns_dash_splitted[0]):
                    source_name, source_id = columns_dash_splitted
        elif len(columns_dash_splitted) == 2:
            source_name, source_id = columns_dash_splitted
        else:
            return

        cursor.execute(QUERY_MOTIF_ANNOTATION_STATEMENT.format(self.motif2tf_version),
                       (source_id, source_name, self.max_motif_similarity_fdr, self.min_orthologous_identity,
                        self.nomenclature_code))

        for (gene_name, orthologous_identity, orthologous_gene_name, orthologous_species, motif_similarity_q_value,
             similar_motif_name, similar_motif_description) in cursor:

            yield (gene_name,
                   orthologous_identity if orthologous_gene_name else None,
                   orthologous_gene_name if orthologous_gene_name else None,
                   orthologous_species if orthologous_gene_name else None,
                   motif_similarity_q_value if similar_motif_name else None,
                   similar_motif_name if similar_motif_name else None,
                   similar_motif_description if similar_motif_name else None)

        cursor.close()

    def _expand_candidate_target_gene_ids(self, candidate_target_gene_ids_plus_ranks, db_id2user_ids):
        for rank, db_id in candidate_target_gene_ids_plus_ranks:
            for user_id in db_id2user_ids[db_id]:
                yield rank, user_id

    def _write_enriched_feature_row(self,
                                    job_id,
                                    database_name,
                                    feature_rank,
                                    feature_id,
                                    feature_description,
                                    auc,
                                    nes,
                                    cluster_number,
                                    candidate_target_gene_ids_plus_ranks):

        candidate_target_gene_ids = ";".join(map(operator.itemgetter(1), candidate_target_gene_ids_plus_ranks))
        candidate_target_gene_ranks = ";".join(map(str,
                                                   map(operator.itemgetter(0), candidate_target_gene_ids_plus_ranks)))

        cursor = self.connection.cursor()
        cursor.execute(INSERT_ENRICHED_FEATURE_STATEMENT,
                       (job_id,
                        database_name,
                        feature_rank,
                        feature_id,
                        feature_description,
                        auc,
                        nes,
                        cluster_number,
                        candidate_target_gene_ids,
                        candidate_target_gene_ranks))

        enriched_feature_row_id = self.connection.insert_id()
        self.connection.commit()
        cursor.close()

        return enriched_feature_row_id

    def _write_feature_annotation_row(self,
                                      enriched_feature_row_id,
                                      tf_name,
                                      orthologous_identity,
                                      orthologous_gene_name,
                                      orthologous_species,
                                      motif_similarity_fdr,
                                      similar_motif_name,
                                      similar_motif_description):

        cursor = self.connection.cursor()
        cursor.execute(INSERT_FEATURE_ANNOTATION_STATEMENT,
                       (enriched_feature_row_id,
                        tf_name,
                        orthologous_identity,
                        orthologous_gene_name,
                        orthologous_species,
                        motif_similarity_fdr,
                        similar_motif_name,
                        similar_motif_description))

        self.connection.commit()
        cursor.close()

    def _write_candidate_enhancer_row(self,
                                      enriched_feature_row_id,
                                      rank,
                                      target_gene_name,
                                      chromosome,
                                      start,
                                      end,
                                      region_name):

        cursor = self.connection.cursor()
        cursor.execute(INSERT_CANDIDATE_ENHANCER_STATEMENT,
                       (enriched_feature_row_id,
                        rank,
                        target_gene_name,
                        chromosome,
                        start,
                        end,
                        region_name))

        self.connection.commit()
        cursor.close()

    ####################################################################################################################
    # Actual analysis operation.
    ####################################################################################################################

    def _find_cluster_number(self, clusters, feature_id):
        for idx, feature_ids in enumerate(clusters):
            if feature_id in feature_ids:
                return idx

        return -1

    def do_analysis(self):
        for database_name in self.database_names:
            logging.info("start analysis for rankings database '{0:s}'.".format(database_name))

            db, db_id2user_ids = self.load_database(database_name)

            if self.region_based_database(database_name):
                mapped_gene_ids_fraction = self.fraction_of_mapped_gene_ids(database_name)

                if mapped_gene_ids_fraction < FRACTION_OF_MAPPED_GENE_IDS:
                    putative_regulatory_regions = self.context.get_regions_for_delineation_name(
                        self.conversion_delineation_name(database_name))
                    all_gene_ids = set(putative_regulatory_regions.feature_ids)
                    lost_gene_ids = self.gene_ids - all_gene_ids

                    message = "Less than {0:d}% (namely {1:d}%) of supplied gene IDs map to ranked genes.\nPlease check that you use the correct species and gene nomenclature.\n\nThe following genes were lost:\n".format(
                        int(FRACTION_OF_MAPPED_GENE_IDS * 100.0),
                        int(mapped_gene_ids_fraction * 100.0))
                    message_lost_gene_ids = "; ".join(lost_gene_ids)

                    raise AssertionError(message + "\n".join(wrap(message_lost_gene_ids)))
            else:
                mapped_gene_ids_fraction = float(db.gene_count) / len(self.gene_ids)

                if mapped_gene_ids_fraction < FRACTION_OF_MAPPED_GENE_IDS:
                    lost_gene_ids = self.gene_ids - set(db.gene_ids)

                    message = "Less than {0:d}% (namely {1:d}%) of supplied gene IDs map to ranked genes.\nPlease check that you use the correct species and gene nomenclature.\n\nThe following genes were lost:\n".format(
                        int(FRACTION_OF_MAPPED_GENE_IDS * 100.0),
                        int(mapped_gene_ids_fraction * 100.0))
                    message_lost_gene_ids = "; ".join(lost_gene_ids)

                    raise AssertionError(message + "\n".join(wrap(message_lost_gene_ids)))

            # Create recovery curves.
            rccs = RecoveryCurves(db, self.auc_threshold, self.rank_threshold, self.nes_threshold)

            # Do clustering analysis.
            enriched_feature_ids = map(operator.itemgetter(0), rccs.enriched_features)

            if self.context.has_stamp():
                logging.info("clustering analysis")

                clusters, garbage = cluster(self.context.get_stamp_motif_database_filename(),
                                            self.context.get_stamp_score_distribution_filename(),
                                            enriched_feature_ids,
                                            cmd=self.context.stamp_command)
                self.context.remove_stamp_intermediate_files()
            else:
                clusters = [enriched_feature_ids]

            # Dump data to database.
            enriched_feature_count = rccs.enriched_feature_count_without_randoms

            for idx, (feature_id, auc, nes) in enumerate(sorted(rccs.enriched_features,
                                                                key=operator.itemgetter(2),
                                                                reverse=True)):
                logging.info("analysis of significant ranking {0:s} (NES = {1:.10g}, {2:d}/{3:d})".format(
                    feature_id,
                    nes,
                    idx + 1,
                    enriched_feature_count)
                )

                enriched_feature_row_id = self._write_enriched_feature_row(
                    self.job_id,
                    database_name,
                    idx + 1,
                    feature_id,
                    self.context.get_description_for_feature(feature_id),
                    auc,
                    nes,
                    self._find_cluster_number(clusters, feature_id),
                    list(self._expand_candidate_target_gene_ids(rccs.get_candidate_target_gene_ids(feature_id),
                                                                db_id2user_ids))
                )

                if self.region_based_database(database_name):
                    scored_regions = self.context.get_regions_for_database_name(database_name)

                    for rank, region_id in rccs.get_candidate_target_gene_ids(feature_id):
                        regions = scored_regions.filter_by_name([region_id])

                        if not regions:
                            continue
                        for gene_id in db_id2user_ids[region_id]:
                            self._write_candidate_enhancer_row(enriched_feature_row_id,
                                                               rank,
                                                               gene_id,
                                                               regions[0].chromosome,
                                                               regions[0].interval[0],
                                                               regions[0].interval[1],
                                                               region_id)

                for (tf_name, orthologous_identity, orthologous_gene_name, orthologous_species, motif_similarity_fdr,
                     similar_motif_name, similar_motif_description) in self._iterate_tfs(feature_id):

                    self._write_feature_annotation_row(enriched_feature_row_id,
                                                       tf_name,
                                                       orthologous_identity,
                                                       orthologous_gene_name,
                                                       orthologous_species,
                                                       motif_similarity_fdr,
                                                       similar_motif_name,
                                                       similar_motif_description)
