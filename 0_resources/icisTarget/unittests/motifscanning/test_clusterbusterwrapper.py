import unittest, os
from cistargetx.motifscanning.clusterbusterwrapper import scan
from cistargetx.conversion.featurelist import FeatureList


CBUST_COMMAND = '/usr/local/bin/cbust'
FASTA_FILENAME = os.path.join(os.path.dirname(__file__), 'twi.6-8.peaks.fasta')
BED_FILENAME = os.path.join(os.path.dirname(__file__), 'twi.6-8.peaks.bed')
CB_FILENAME = os.path.join(os.path.dirname(__file__), 'flyfactorsurvey-twi_FlyReg_FBgn0003900.cb')


class ClusterBusterWrapperTestCase(unittest.TestCase):
    def setUp(self):
        self.feature_list = FeatureList.from_bed_file(BED_FILENAME)

    def test_scan(self):
        features = scan(CB_FILENAME, FASTA_FILENAME, CBUST_COMMAND, self.feature_list)
        for feature in features:
            print str(feature) + " | len = " + str(len(feature))


if __name__ == '__main__':
    unittest.main()
