import unittest, tempfile, os
from cistargetx.motifscanning.motifscanner import MotifCollection

MOTIF_FOLDER = '${DATADIR}/motifs/cbust-collection-nov2010/singletons/'


class MotifScannerTestCase(unittest.TestCase):
    def test_motif_collection_init(self):
        collection = MotifCollection.load_from_folder(os.path.expandvars(MOTIF_FOLDER))
        self.assertEquals(3731, len(collection))

    def test_motif_collection_reduce(self):
        collection = MotifCollection.load_from_folder(os.path.expandvars(MOTIF_FOLDER)).reduce(['jaspar4-PH0054.1-Hoxa7_1', 'jaspar4-PH0056.1-Hoxa9'])
        self.assertEquals(2, len(collection))

    def test_motif_collection_write(self):
        collection = MotifCollection.load_from_folder(os.path.expandvars(MOTIF_FOLDER)).reduce(['jaspar4-PH0054.1-Hoxa7_1', 'jaspar4-PH0056.1-Hoxa9'])
        temp_filename = tempfile.mktemp()

        collection.write(temp_filename)
        with open(temp_filename, 'r') as input_fh:
            print input_fh.read()

        os.remove(temp_filename)


if __name__ == '__main__':
    unittest.main()
