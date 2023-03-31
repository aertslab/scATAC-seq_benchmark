import copy
import logging
import numpy
import operator
import os
import re
import scipy
import sqlite3

from lrucache import LRUCache
from orderstatistics import combine_via_order_statistics
from utils import derive_dtype


GENE_ID_COUNT_QUERY = r"SELECT COUNT(*) FROM rankings;"
GENE_IDS_QUERY = r"SELECT geneID FROM rankings WHERE geneID IN ({0:s}) ORDER BY geneID;"
ALL_GENE_IDS_QUERY = r"SELECT geneID FROM rankings ORDER BY geneID;"
RANKINGS_QUERY = r"SELECT geneID, ranking FROM rankings WHERE geneID IN ({0:s}) ORDER BY geneID;"
ALL_RANKINGS_QUERY = r"SELECT geneID, ranking FROM rankings ORDER BY geneID;"
FEATURE_IDS_QUERY = r"SELECT motifName FROM motifs ORDER BY idx;"

CREATE_TABLE_STATEMENTS = r"""
DROP TABLE IF EXISTS rankings;
DROP TABLE IF EXISTS motifs;
CREATE TABLE rankings (geneID VARCHAR(255), ranking BLOB);
CREATE TABLE motifs (motifName VARCHAR(255), idx INTEGER);
"""

CREATE_INDEX_STATEMENT = r"CREATE UNIQUE INDEX id ON rankings (geneID)"
INSERT_FEATURE_STATEMENT = r"INSERT INTO motifs (idx, motifName) VALUES (?, ?);"
INSERT_RANKING_STATEMENT = r"INSERT INTO rankings (geneID, ranking) VALUES (?, ?);"


class RankingsDatabase:
    # TODO: explicitly chosen for a closed arithmetic interface to rankings databases. Converting
    # TODO: this to a mutator interface might be beneficient for the performance (no copying of data).
    @staticmethod
    def _fetch_feature_ids(filename):
        """ Fetch features from SQLite3 database with supplied filename. """
        if not os.path.exists(filename):
            raise ValueError("Database {0:s} doesn't exist.".format(filename))
        with sqlite3.connect(filename) as db:
            cursor = db.cursor()
            all_feature_ids = map(operator.itemgetter(0), cursor.execute(FEATURE_IDS_QUERY).fetchall())
            cursor.close()
        return all_feature_ids

    @staticmethod
    def _quoted_csv(values):
        # Escape single quotes (') by using (''), because sometimes ID's contain a single quote.
        def quote(value):
            return "'" + value.replace("'", "''") + "'"

        return ','.join(map(quote, values))

    @staticmethod
    def _create_empty_db(name):
        return RankingsDatabase(name,
                                numpy.empty(shape=0, dtype='|S255'),
                                numpy.empty(shape=0, dtype='|S255'),
                                numpy.empty(shape=(0, 0), dtype=numpy.float64),
                                numpy.zeros(shape=0, dtype=numpy.int_),
                                {0: name},
                                0)

    @staticmethod
    def load_with_feature_patterns(filename, name, feature_patterns, gene_ids=None):
        masks = [re.compile(pattern) for pattern in feature_patterns]
        all_feature_ids = RankingsDatabase._fetch_feature_ids(filename)
        feature_ids = filter(lambda feature_id: any(m.match(feature_id) for m in masks),
                             all_feature_ids)
        if feature_ids:
            logging.info("For {0:s} only {1:d} of the {2:d} features remain after application of masks.".format(
                filename, len(feature_ids), len(all_feature_ids)))
            return RankingsDatabase.load(filename, name, feature_ids, gene_ids)
        else:
            logging.info("For {0:s} no features remain after application of masks.".format(filename))
            return RankingsDatabase._create_empty_db(name)

    @staticmethod
    def load(filename, name, feature_ids=None, gene_ids=None):
        logging.info("Start: loading database {0:s} ( {1:s} / {2:s}) ) in memory.".format(
            filename,
            "No feature ID" if feature_ids is None else "{0:d} feature IDs".format(len(feature_ids)),
            "No gene/region ID" if gene_ids is None else "{0:d} gene/region IDs".format(len(gene_ids)))
        )

        """ Fetch rankings from SQLite3 database with supplied filename. """
        if not os.path.exists(filename):
            raise ValueError("Database {0:s} doesn't exist.".format(filename))

        with sqlite3.connect(filename) as db:
            cursor = db.cursor()

            all_feature_ids = map(operator.itemgetter(0), cursor.execute(FEATURE_IDS_QUERY).fetchall())

            if feature_ids:
                selected_feature_idxs_ids = [(all_feature_ids.index(feature_id), feature_id)
                                             for feature_id in feature_ids
                                             if feature_id in all_feature_ids]
                if selected_feature_idxs_ids:
                    feature_idxs, feature_ids = zip(*selected_feature_idxs_ids)
                    feature_idxs = list(feature_idxs)
                else:
                    feature_idxs, feature_ids = None, []
            else:
                feature_idxs, feature_ids = None, all_feature_ids

            feature_id_count = len(feature_ids)

            total_gene_id_count = int(cursor.execute(GENE_ID_COUNT_QUERY).fetchone()[0])
            dtype = derive_dtype(total_gene_id_count)

            if gene_ids:
                cursor.execute(GENE_IDS_QUERY.format(RankingsDatabase._quoted_csv(gene_ids)))
                gene_ids = set(gene_id for gene_id, in cursor)
                gene_id_count = len(gene_ids)
            else:
                gene_ids = []
                gene_id_count = total_gene_id_count

            # Fill rankings array with zeros, but do not use numpy.zeros as it does not preallocate the whole array.
            rankings = numpy.full(shape=(gene_id_count, feature_id_count), fill_value=0, dtype=dtype)
            if not feature_id_count or not gene_id_count:
                cursor.close()
                return RankingsDatabase._create_empty_db(name)
            else:
                try:
                    if gene_ids:
                        cursor.execute(RANKINGS_QUERY.format(RankingsDatabase._quoted_csv(gene_ids)))
                    else:
                        cursor.execute(ALL_RANKINGS_QUERY)
                except sqlite3.OperationalError as e:
                    if str(e) == 'disk I/O error':
                        raise IOError("An disk I/O error occurred during execution of a large SQLite query, try TMPDIR=/path/with/enough/free/space <program> <args>")
                    else:
                        raise sqlite3.OperationalError(str(e))

                row_idx = 0
                gene_ids = []

                if feature_idxs:
                    for ID, ranking in cursor:
                        rankings[row_idx, :] = numpy.frombuffer(ranking, dtype=dtype)[feature_idxs]
                        gene_ids.append(ID)
                        row_idx += 1
                else:
                    for ID, ranking in cursor:
                        rankings[row_idx, :] = numpy.frombuffer(ranking, dtype=dtype)
                        gene_ids.append(ID)
                        row_idx += 1

                cursor.close()

                logging.info("End: loading database {0:s} in memory.".format(filename))

                return RankingsDatabase(name,
                                        numpy.array(feature_ids, dtype='|S255'),
                                        numpy.array(gene_ids, dtype='|S255'),
                                        rankings,
                                        numpy.zeros(shape=feature_id_count, dtype=numpy.int_),
                                        {0: name},
                                        total_gene_id_count,
                                        filename)

    @staticmethod
    def create(name, feature_ids, gene_ids, rankings):
        return RankingsDatabase(name,
                                numpy.array(feature_ids, dtype='|S255'),
                                numpy.array(gene_ids, dtype='|S255'),
                                rankings,
                                numpy.zeros(shape=len(feature_ids), dtype=numpy.int_),
                                {0: name},
                                len(gene_ids))

    def __init__(self, name, feature_ids, gene_ids, rankings, groups, group2name, total_gene_count, filename=None):
        self.name = name
        self.feature_ids = feature_ids
        self.gene_ids = gene_ids
        self.rankings = rankings
        self.groups = groups
        self.group2name = group2name
        self.total_gene_count = total_gene_count
        self.filename = filename
        self.preloaded_top_ranked_genes = False

        assert self.name.strip() != ""
        assert self.gene_count == self.rankings.shape[0]
        assert self.feature_count == self.rankings.shape[1]
        assert self.groups.size == self.feature_count
        assert self.total_gene_count >= self.gene_count

    def __str__(self):
        return "RankingsDatabase(name = '{0:s}', #gene IDs = {1:d}, #features = {2:d})".format(
            self.name, self.gene_count, self.feature_count)


    @property
    def dtype(self):
        return self.rankings.dtype

    @property
    def gene_count(self):
        return len(self.gene_ids)

    @property
    def feature_count(self):
        return len(self.feature_ids)

    @property
    def is_empty(self):
        return self.feature_count == 0 or self.gene_count == 0

    def _find_index(self, feature_id):
        try:
            return list(self.feature_ids).index(feature_id)
        except ValueError:
            return -1

    def get_group_name(self, feature_id):
        feature_idx = self._find_index(feature_id)

        if feature_idx >= 0:
            return self.group2name[self.groups[feature_idx]]
        else:
            return None

    def get_ranking(self, feature_id):
        feature_idx = self._find_index(feature_id)

        if feature_idx >= 0:
            return self.rankings[:, feature_idx].flat
        else:
            return None

    def get_single_ranking(self, gene_id, feature_id):
        row_idx = numpy.where(self.gene_ids == gene_id)[0][0]
        col_idx = self._find_index(feature_id)

        if row_idx >= 0 and col_idx >= 0:
            return self.rankings[row_idx, col_idx]
        else:
            return None

    def write(self, filename):
        with sqlite3.connect(filename) as db:
            db.text_factory = str
            cursor = db.cursor()
            cursor.executescript(CREATE_TABLE_STATEMENTS)

            for feature_idx, feature_id in enumerate(self.feature_ids):
                cursor.execute(INSERT_FEATURE_STATEMENT, (feature_idx, feature_id))

            for row_idx, gene_id in enumerate(self.gene_ids):
                # Caveat: buffer() is absent from python 3.0.
                # To avoid "ValueError: could not convert BLOB to buffer",
                # flatten the array which makes a contiguous copy of the array.
                cursor.execute(INSERT_RANKING_STATEMENT, (gene_id, buffer(self.rankings[row_idx, :].flatten(order='C'))))

            cursor.execute(CREATE_INDEX_STATEMENT)
            cursor.close()

    def add_with_optimized_dtype(self, other):
        assert self.gene_count == other.gene_count
        assert numpy.all(numpy.in1d(self.gene_ids, other.gene_ids))
        assert self.total_gene_count == other.total_gene_count

        db_name = "(" + self.name + "+" + other.name + ")"

        feature_ids = numpy.empty(shape=self.feature_count + other.feature_count, dtype='|S255')
        feature_ids[:self.feature_count] = self.feature_ids
        feature_ids[self.feature_count:] = other.feature_ids

        # TODO: weird bug when using adddatabases.py: "could not convert BLOB to buffer".
        #rankings = numpy.hstack( (self.rankings, other.rankings) )
        feature_count = len(feature_ids)
        # Fill rankings array with zeros, but do not use numpy.zeros as it does not preallocate the whole array.
        rankings = numpy.full(shape=(self.gene_count, feature_count), fill_value=0, dtype=derive_dtype(self.gene_count))
        rankings[:, :self.feature_count] = self.rankings
        rankings[:, self.feature_count:] = other.rankings

        groups = numpy.empty(shape=feature_ids.shape, dtype=numpy.int_)
        groups[:self.feature_count] = self.groups
        group_inc = self.groups.max() + 1
        groups[self.feature_count:] = (other.groups + group_inc)
        group2name = copy.deepcopy(self.group2name)

        for group, name in other.group2name.iteritems():
            group2name[group + group_inc] = name

        return RankingsDatabase(db_name,
                                feature_ids,
                                self.gene_ids,
                                rankings,
                                groups,
                                group2name,
                                self.total_gene_count)

    def retain_features_with_optimized_dtype(self, feature_ids):
        if not feature_ids or len(feature_ids) == 0:
            return RankingsDatabase._create_empty_db(self.name)

        keep_feature_idxs = numpy.in1d(self.feature_ids, numpy.array(list(feature_ids), '|S255'))
        all_feature_idxs = numpy.arange(self.feature_count)
        feature_idxs = list(all_feature_idxs[keep_feature_idxs])
        assert len(feature_ids) == len(feature_idxs), "Reduce feature ids: number of IDs is not equal to the number of idxs."
        feature_count = len(feature_ids)
        # Fill rankings array with zeros, but do not use numpy.zeros as it does not preallocate the whole array.
        rankings = numpy.full(shape=(self.gene_count, feature_count), fill_value=0, dtype=derive_dtype(self.gene_count))
        rankings[:, ] = self.rankings[:, feature_idxs]
        groups = self.groups[feature_idxs]
        group2name = dict((group, name) for group, name in self.group2name.iteritems() if group in groups)

        return RankingsDatabase(self.name,
                                self.feature_ids[feature_idxs],
                                self.gene_ids, rankings,
                                groups,
                                group2name,
                                self.total_gene_count,
                                self.filename)

    def __add__(self, other):
        assert self.gene_count == other.gene_count
        assert numpy.all(numpy.in1d(self.gene_ids, other.gene_ids))
        assert self.total_gene_count == other.total_gene_count

        db_name = "(" + self.name + "+" + other.name + ")"

        feature_ids = numpy.empty(shape=self.feature_count + other.feature_count, dtype='|S255')
        feature_ids[:self.feature_count] = self.feature_ids
        feature_ids[self.feature_count:] = other.feature_ids

        rankings = numpy.hstack((self.rankings, other.rankings))

        groups = numpy.empty(shape=feature_ids.shape, dtype=numpy.int_)
        groups[:self.feature_count] = self.groups
        group_inc = self.groups.max() + 1
        groups[self.feature_count:] = (other.groups + group_inc)
        group2name = copy.deepcopy(self.group2name)

        for group, name in other.group2name.iteritems():
            group2name[group + group_inc] = name

        return RankingsDatabase(db_name,
                                feature_ids,
                                self.gene_ids,
                                rankings,
                                groups, group2name,
                                self.total_gene_count)

    def reduce_feature_ids(self, feature_ids):
        if not feature_ids or len(feature_ids) == 0:
            return RankingsDatabase._create_empty_db(self.name)

        keep_feature_idxs = numpy.in1d(self.feature_ids, numpy.array(list(feature_ids), '|S255'))
        all_feature_idxs = numpy.arange(self.feature_count)
        feature_idxs = list(all_feature_idxs[keep_feature_idxs])
        assert len(feature_ids) == len(feature_idxs), "Reduce feature ids: number of IDs is not equal to the number of idxs."
        groups = self.groups[feature_idxs]
        group2name = dict((group, name) for group, name in self.group2name.iteritems() if group in groups)

        return RankingsDatabase(self.name,
                                self.feature_ids[feature_idxs],
                                self.gene_ids,
                                self.rankings[:, feature_idxs],
                                groups,
                                group2name,
                                self.total_gene_count,
                                self.filename)

    def reduce_gene_ids(self, gene_ids):
        keep_gene_idxs = numpy.in1d(self.gene_ids, numpy.array(list(gene_ids), '|S255'))
        all_gene_idxs = numpy.arange(self.gene_count)
        gene_idxs = all_gene_idxs[keep_gene_idxs]
        assert len(gene_ids) == len(gene_idxs), "Reduce gene ids: number of IDs is not equal to the number of idxs."

        return RankingsDatabase(self.name,
                                self.feature_ids,
                                self.gene_ids[gene_idxs],
                                self.rankings[gene_idxs, :],
                                self.groups,
                                self.group2name,
                                self.total_gene_count,
                                self.filename)

    def reduce_gene_ids_with_optimized_dtype(self, gene_ids):
        gene_count = len(gene_ids)

        if not gene_ids or gene_count == 0:
            return RankingsDatabase._create_empty_db(self.name)

        keep_gene_idxs = numpy.in1d(self.gene_ids, numpy.array(list(gene_ids), '|S255'))
        all_gene_idxs = numpy.arange(self.gene_count)
        gene_idxs = all_gene_idxs[keep_gene_idxs]
        assert gene_count == len(gene_idxs), "Reduce gene ids: number of IDs is not equal to the number of idxs."

        # Select the requested gene ids from the rankings matrix and rebuild the ranking and convert to the right dtype.
        rankings = numpy.argsort(numpy.argsort(self.rankings[gene_idxs, :], axis=0),
                                 axis=0
                                 ).astype(dtype=derive_dtype(self.gene_count))

        return RankingsDatabase(self.name,
                                self.feature_ids,
                                self.gene_ids[gene_idxs],
                                rankings,
                                self.groups,
                                self.group2name,
                                gene_count,
                                self.filename)

    def __mul__(self, other):
        """
        Gene IDs should be ordered alphabetically.
        And only works for databases with all gene IDs loaded into memory.
        """
        assert self.gene_count == other.gene_count
        assert numpy.all(numpy.in1d(self.gene_ids, other.gene_ids))
        assert self.total_gene_count == other.total_gene_count
        assert self.gene_count == self.total_gene_count

        db_name = "(" + self.name + "*" + other.name + ")"

        feature_count = self.feature_count + other.feature_count
        # TODO: scipy.misc appears to be unavailable in version 0.10 ...
        combination_count = scipy.misc.comb(feature_count, 2, exact=1)
        # Fill rankings array with zeros, but do not use numpy.zeros as it does not preallocate the whole array.
        rankings = numpy.full(shape=(self.gene_count, combination_count), fill_value=0, dtype=self.dtype)
        feature_ids = []
        idx = 0
        rankings_pair = numpy.full(shape=(self.gene_count, 2), fill_value=0, dtype=numpy.float)

        for feature_idx1 in range(0, self.feature_count):
            for feature_idx2 in range(0, other.feature_count):
                feature_name = self.feature_ids[feature_idx1] + "*" + other.feature_ids[feature_idx2]
                feature_ids.append(feature_name)
                logging.info("combining {0:s} via orderstatistics".format(feature_name))

                # The relative rankings have to be created from 1-based absolute rankings,
                # to avoid zeros in the relative rankings (annihilator in multiplication).
                rankings_pair[:, 0] = (self.rankings[:, feature_idx1] + 1.0) / float(self.gene_count)
                rankings_pair[:, 1] = (other.rankings[:, feature_idx2] + 1.0) / float(self.gene_count)

                # This only works because gene_ids are ordered alphabetically.
                combined_gene_ids = self.gene_ids[numpy.argsort(combine_via_order_statistics(rankings_pair))]
                rankings[:, idx] = numpy.argsort(combined_gene_ids)
                idx += 1

        return RankingsDatabase(db_name,
                                numpy.array(feature_ids, dtype='|S255'),
                                self.gene_ids,
                                rankings,
                                numpy.zeros(shape=combination_count, dtype=numpy.int_),
                                {0: db_name},
                                self.total_gene_count)

    def __pow__(self, power, modulo=None):
        if power != 2:
            raise ValueError('Only power of 2 supported.')

        db_name = self.name + "**2"
        combination_count = scipy.misc.comb(self.feature_count, 2, exact=1)
        # Fill rankings array with zeros, but do not use numpy.zeros as it does not preallocate the whole array.
        rankings = numpy.full(shape=(self.gene_count, combination_count), fill_value=0, dtype=self.dtype)
        feature_ids = []
        idx = 0
        rankings_pair = numpy.empty(shape=(self.gene_count, 2), dtype=numpy.float)

        for feature_idx1 in range(0, self.feature_count):
            for feature_idx2 in range(feature_idx1 + 1, self.feature_count):
                feature_name = self.feature_ids[feature_idx1] + "*" + self.feature_ids[feature_idx2]
                feature_ids.append(feature_name)
                logging.info("combining {0:s} via orderstatistics".format(feature_name))

                # The relative rankings have to be created from 1-based absolute rankings,
                # to avoid zeros in the relative rankings (annihilator in multiplication).
                rankings_pair[:, 0] = (self.rankings[:, feature_idx1] + 1.0) / float(self.gene_count)
                rankings_pair[:, 1] = (self.rankings[:, feature_idx2] + 1.0) / float(self.gene_count)

                # This only works because gene_ids are ordered alphabetically.
                combined_gene_ids = self.gene_ids[numpy.argsort(combine_via_order_statistics(rankings_pair))]
                rankings[:, idx] = numpy.argsort(combined_gene_ids)
                idx += 1

        return RankingsDatabase(db_name,
                                numpy.array(feature_ids, dtype='|S255'),
                                self.gene_ids, rankings,
                                numpy.zeros(shape=combination_count, dtype=numpy.int_),
                                {0: db_name},
                                self.total_gene_count)

    def collapse(self):
        logging.info("collapsing {0:s} via orderstatistics".format(self.name))

        db_name = "<" + self.name + ">"
        feature_ids = ['*'.join(self.feature_ids)]

        # The relative rankings have to be created from 1-based absolute rankings,
        #  to avoid zeros in the relative rankings (annihilator in multiplication).
        relative_rankings = numpy.float64(self.rankings)
        relative_rankings += 1.0
        relative_rankings /= float(self.gene_count)

        # This only works because gene_ids are ordered alphabetically.
        combined_gene_ids = self.gene_ids[numpy.argsort(combine_via_order_statistics(relative_rankings))]
        rankings = numpy.argsort(combined_gene_ids)[:, numpy.newaxis]

        return RankingsDatabase(db_name,
                                numpy.array(feature_ids, dtype='|S255'),
                                self.gene_ids,
                                rankings,
                                numpy.zeros(shape=1, dtype=numpy.int_),
                                {0: db_name},
                                self.total_gene_count)

    @property
    def has_top_regions(self):
        return self.filename is not None

    def load_top_regions(self, rank_threshold, feature_ids=None):
        if self.has_loaded_top_regions:
            return

        assert rank_threshold > 0
        assert self.has_top_regions, "Top ranked regions can not be queried because filename for database is not specified."
        logging.info("loading top {0:d} ranked genes of database '{1:s}' into memory".format(rank_threshold, self.name))

        db = RankingsDatabase.load(self.filename, self.name)
        self.feature_id2top_rank_genes = dict()

        for feature_id in (feature_ids if feature_ids is not None else db.feature_ids):
            rankings = db.get_ranking(feature_id)
            ranking_idxs = rankings < rank_threshold
            self.feature_id2top_rank_genes[feature_id] = zip(rankings[ranking_idxs] + 1, db.gene_ids[ranking_idxs])

        self.preloaded_top_ranked_genes = True

    @property
    def has_loaded_top_regions(self):
        return self.preloaded_top_ranked_genes

    def clear_top_regions(self):
        if not self.has_loaded_top_regions: return
        del self.feature_id2top_rank_genes
        self.preloaded_top_ranked_genes = False

    def get_top_regions(self, rank_threshold, feature_id):
        assert self.has_loaded_top_regions
        return filter(lambda rank_gene_id: rank_gene_id[0] <= rank_threshold,
                      self.feature_id2top_rank_genes[feature_id])


def _load_database(filename, name, gene_ids=None, feature_masks=[]):
    if feature_masks:
        return RankingsDatabase.load_with_feature_patterns(filename, name, feature_masks, gene_ids=gene_ids)
    else:
        return RankingsDatabase.load(filename, name, gene_ids=gene_ids)


DATABASE_CACHE = None


def _create_combination_database(filenames, gene_ids, features_to_combine=[], control_db=None,
                                 name_translation=os.path.basename, combination_method='pairwise',
                                 db_cache_size=0):
    dbs = []

    for filename in filenames:
        if db_cache_size > 0:
            global DATABASE_CACHE

            if not DATABASE_CACHE:
                DATABASE_CACHE = LRUCache(db_cache_size)

            if filename not in DATABASE_CACHE:
                logging.info("'{0:s}' is not in cache. Loading database into memory.".format(filename))
                DATABASE_CACHE[filename] = RankingsDatabase.load(filename, name_translation(filename))
            else:
                logging.info("'{0:s}' already in cache. Using cached version.".format(filename))

            filtered_feature_ids = filter(lambda feature_id:
                                          feature_id in DATABASE_CACHE[filename].feature_ids, features_to_combine)

            db = DATABASE_CACHE[filename].reduce_feature_ids(filtered_feature_ids)
        else:
            db = RankingsDatabase.load(filename, name_translation(filename), feature_ids=features_to_combine)
        if not db.is_empty:
            dbs.append(db)
    if len(dbs) > 0:
        if combination_method == 'pairwise':
            db = reduce(operator.add, dbs) ** 2
        elif combination_method == 'collapse':
            db = reduce(operator.add, dbs).collapse()
        else:
            raise ValueError

        db = db.reduce_gene_ids(gene_ids)

        if control_db:
            db += RankingsDatabase.load(control_db, name_translation(control_db), gene_ids=gene_ids)

        return db
    else:
        return None


def load_databases(filenames, gene_ids=None, feature_masks=[], features_to_combine=[],
                   enrichment_within_db=False, control_db=None,
                   name_translation=os.path.basename,
                   combination_method='pairwise',
                   db_cache_size=0):
    assert combination_method in ['pairwise', 'collapse']

    assert gene_ids.issubset([]) is False, "No genes are given."

    filename2db = dict()

    for filename in filenames:
        db = _load_database(filename, name_translation(filename), gene_ids, feature_masks)
        if not db.is_empty:
            assert db.gene_count > 1, ("Please, check carefully that the correct species and input type were selected "
                                       "as only one of the provided gene or i-cisTarget region IDs could be found in "
                                       "the rankings database.")
            filename2db[filename] = db

    if features_to_combine:
        tmp_control_db = control_db if enrichment_within_db and control_db not in filename2db else None
        combination_db = _create_combination_database(filenames, gene_ids, features_to_combine, tmp_control_db,
                                                      name_translation, combination_method, db_cache_size)
    else:
        combination_db = None

    if not filename2db:
        raise AssertionError(
            "Please, check carefully that the correct species and input type were selected as none of the provided "
            "gene or i-cisTarget region IDs could be found in the rankings database.")
    elif combination_db:
        if enrichment_within_db and control_db in filename2db:
            dbs = [filename2db[control_db] + combination_db]
            dbs.extend(db for filename, db in filename2db.iteritems() if filename != control_db)
            return dbs
        else:
            dbs = list(filename2db.values())
            dbs.append(combination_db)
            
            if enrichment_within_db:
                return dbs
            else:
                return [reduce(operator.add, dbs)]
    elif enrichment_within_db:
        return list(filename2db.values())
    else:
        return [reduce(operator.add, list(filename2db.values()))]
