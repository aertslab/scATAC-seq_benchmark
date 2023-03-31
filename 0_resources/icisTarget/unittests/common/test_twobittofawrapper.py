import unittest, os
from cistargetx.conversion.featurelist import FeatureList
from cistargetx.common.twobittofawrapper import convert_to_fasta

TEST_BED_FILENAME = os.path.join(os.path.dirname(__file__), 'test.bed')
TWO_BIT_FILENAME = os.path.expandvars('${DATADIR}/genomes/dm3.2bit')
COMMAND = '/usr/local/bin/twoBitToFa'


class MyTestCase(unittest.TestCase):
    def setUp(self):
        self.feature_list = FeatureList.from_bed_file(TEST_BED_FILENAME)

    def test_convert_to_fasta(self):
        id2seq = convert_to_fasta(self.feature_list, TWO_BIT_FILENAME, COMMAND)
        self.assertEqual(10, len(id2seq))
        self.assertTrue('twi_2-4_5.69897000434764_1.39' in id2seq)
        self.assertTrue('twi_2-4_5.37365963261988_1.24' in id2seq)
        self.assertTrue('twi_2-4_5.35654732351228_1.24' in id2seq)
        self.assertTrue('twi_2-4_6.0177287669379_1.95' in id2seq)
        self.assertTrue('twi_2-4_5.96657624451842_1.55' in id2seq)
        self.assertTrue('twi_2-4_5.96657624451842_2.22' in id2seq)
        self.assertTrue('twi_2-4_5.97469413473639_1.51' in id2seq)
        self.assertTrue('twi_2-4_4.89790947448997_1.23' in id2seq)
        self.assertTrue('twi_2-4_5.66354026614676_1.49' in id2seq)
        self.assertTrue('twi_2-4_5.84771165560785_1.56' in id2seq)


if __name__ == '__main__':
    unittest.main()
