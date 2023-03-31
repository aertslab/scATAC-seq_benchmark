import unittest, os
from cistargetx.common.rankingsdatabase import RankingsDatabase
from cistargetx.recoveryanalysis.recoverycurves import RecoveryCurves

TEST_DATABASE = os.path.expandvars('${DATADIR}/cistargetx/rankings/hg19-limited-upstream10kb-5utr-intron1-10species.db')
TEST_QUERY_FILENAME = os.path.join(os.path.dirname(__file__), 'p53_Kannan.hgnc')


def load_ids(filename):
    with open(filename, 'r') as input_fh:
        for line in input_fh:
            yield line.rstrip()

class RecoveryCurvesTestCase(unittest.TestCase):
    def setUp(self):
        self.db = RankingsDatabase.load(TEST_DATABASE, os.path.basename(TEST_DATABASE), gene_ids=list(load_ids(TEST_QUERY_FILENAME)))
        self.curves = RecoveryCurves(self.db, rank_threshold=1000, auc_threshold=0.03, enrichment_threshold=2.5)

    def test_get_enriched_feature_count(self):
        self.assertEqual(self.curves.enriched_feature_count, 106)

if __name__ == '__main__':
    unittest.main()
