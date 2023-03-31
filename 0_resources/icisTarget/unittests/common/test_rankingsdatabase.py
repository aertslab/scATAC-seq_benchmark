import os
import unittest
from cistargetx.common.rankingsdatabase import RankingsDatabase

TEST_DATABASE = os.path.expandvars('${DATADIR}/cistargetx/rankings/hg19-limited-upstream10kb-5utr-intron1-10species.db')

class RankingsDatabaseTestCase(unittest.TestCase):
    def test_load(self):
        db = RankingsDatabase.load(TEST_DATABASE, os.path.basename(TEST_DATABASE))
        self.assertEqual(db.name, os.path.basename(TEST_DATABASE))
        self.assertEqual(db.feature_count, 3731)
        self.assertEqual(db.gene_count, 21798)
        self.assertEqual(db.get_group_name('transfac10.1-M01651-V-P53_03'), os.path.basename(TEST_DATABASE))

    def test_load_with_gene_ids(self):
        genes = ['TP53', 'E2F1']
        db = RankingsDatabase.load(TEST_DATABASE, os.path.basename(TEST_DATABASE), gene_ids=genes)
        self.assertEqual(db.feature_count, 3731)
        self.assertEqual(db.gene_count, 2)
        self.assertTrue(db.gene_ids[0] in genes)
        self.assertTrue(db.gene_ids[1] in genes)

    def test_load_with_feature_ids(self):
        features = ['transfac10.1-M01651-V-P53_03', 'transfac10.1-M01637-F-RAP1_02']
        db = RankingsDatabase.load(TEST_DATABASE,
                                   os.path.basename(TEST_DATABASE),
                                   feature_ids=features)
        self.assertEqual(db.feature_count, 2)
        self.assertEqual(db.gene_count, 21798)
        self.assertTrue(db.feature_ids[0] in features)
        self.assertTrue(db.feature_ids[1] in features)

    def test_reduce_feature_ids(self):
        features = ['transfac10.1-M01651-V-P53_03', 'transfac10.1-M01637-F-RAP1_02']
        db = RankingsDatabase.load(TEST_DATABASE, os.path.basename(TEST_DATABASE)).reduce_feature_ids(features)
        self.assertEqual(db.feature_count, 2)
        self.assertEqual(db.gene_count, 21798)
        self.assertTrue(db.feature_ids[0] in features)
        self.assertTrue(db.feature_ids[1] in features)

    def test_reduce_gene_ids(self):
        genes = ['TP53', 'E2F1']
        db = RankingsDatabase.load(TEST_DATABASE, os.path.basename(TEST_DATABASE)).reduce_gene_ids(genes)
        self.assertEqual(db.feature_count, 3731)
        self.assertEqual(db.gene_count, 2)
        self.assertEqual(db.total_gene_count, 21798)
        self.assertTrue(db.gene_ids[0] in genes)
        self.assertTrue(db.gene_ids[1] in genes)

    def test_add(self):
        features = ['transfac10.1-M01651-V-P53_03', 'transfac10.1-M01637-F-RAP1_02']
        db1 = RankingsDatabase.load(TEST_DATABASE,
                                   os.path.basename(TEST_DATABASE),
                                   feature_ids=[features[0]])
        db2 = RankingsDatabase.load(TEST_DATABASE,
                                   os.path.basename(TEST_DATABASE),
                                   feature_ids=[features[1]])
        db = db1 + db2
        self.assertEqual(db.feature_count, 2)
        self.assertEqual(db.gene_count, 21798)
        self.assertTrue(db.feature_ids[0] in features)
        self.assertTrue(db.feature_ids[1] in features)

    def test_mul(self):
        features = ['transfac10.1-M01651-V-P53_03', 'transfac10.1-M01637-F-RAP1_02']
        db1 = RankingsDatabase.load(TEST_DATABASE,
                                   os.path.basename(TEST_DATABASE),
                                   feature_ids=[features[0]])
        db2 = RankingsDatabase.load(TEST_DATABASE,
                                   os.path.basename(TEST_DATABASE),
                                   feature_ids=[features[1]])
        db = db1 * db2
        self.assertEqual(db.feature_count, 1)
        self.assertEqual(db.gene_count, 21798)
        self.assertTrue(db.feature_ids[0], features[0] + '*' + features[1])
        self.assertEqual(db.get_group_name(features[0] + '*' + features[1]), "(" + db1.name + '*' + db2.name + ")")
        self.assertEqual(len(db.group2name), 1)

    def test_load_top_regions(self):
        db = RankingsDatabase.load(TEST_DATABASE,
                                   os.path.basename(TEST_DATABASE),
                                   feature_ids=['transfac10.1-M01651-V-P53_03'])
        db.load_top_regions(5000)

if __name__ == '__main__':
    unittest.main()
