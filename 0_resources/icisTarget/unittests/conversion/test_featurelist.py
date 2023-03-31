import unittest
from cistargetx.conversion.featurelist import Feature, FeatureList

class FeatureTestCase(unittest.TestCase):
    def setUp(self):
        self.f = Feature('chr1', 20, 30, 'ID1')

    def test_no_name(self):
        f = Feature('chr1', 20, 30)
        self.assertEqual(f.name, Feature.DUMMY_NAME)

    def test_len(self):
        self.assertEqual(len(self.f), 10)

    def test_contains(self):
        self.assertTrue(Feature('chr1', 20, 30, 'ID2') in self.f)
        self.assertTrue(Feature('chr1', 19, 30, 'ID3') not in self.f)
        self.assertTrue(Feature('chr1', 21, 23, 'ID4') in self.f)
        self.assertTrue(Feature('chr2', 21, 23, 'ID4') not in self.f)

    def test_from_string(self):
        f = Feature.from_string("chr1 20 30 ID1")
        self.assertEqual(f.chromosome, 'chr1')
        self.assertEqual(f.interval, (20, 30))
        self.assertEqual(f.name, 'ID1')
        f = Feature.from_string("chr1\t20\t30\tID1")
        self.assertEqual(f.chromosome, 'chr1')
        self.assertEqual(f.interval, (20, 30))
        self.assertEqual(f.name, 'ID1')
        f = Feature.from_string("chr1 \t 20 \t 30 \t ID1  ")
        self.assertEqual(f.chromosome, 'chr1')
        self.assertEqual(f.interval, (20, 30))
        self.assertEqual(f.name, 'ID1')
        f = Feature.from_string("chr1 \t 20 \t 30 \t ")
        self.assertEqual(f.chromosome, 'chr1')
        self.assertEqual(f.interval, (20, 30))
        self.assertEqual(f.name, Feature.DUMMY_NAME)

    def test_has_overlap_with(self):
        self.assertTrue(Feature('chr1', 20, 30, 'ID2').has_overlap_with(self.f))
        self.assertFalse(Feature('chr1', 10, 20, 'ID2').has_overlap_with(self.f))
        self.assertTrue(Feature('chr1', 10, 21, 'ID2').has_overlap_with(self.f))
        self.assertFalse(Feature('chr1', 10, 19, 'ID2').has_overlap_with(self.f))
        self.assertFalse(Feature('chr2', 20, 30, 'ID2').has_overlap_with(self.f))
        self.assertFalse(Feature('chr1', 30, 32, 'ID2').has_overlap_with(self.f))

    def test_get_overlap_in_bp_with(self):
        self.assertEqual(Feature('chr1', 10, 20, 'ID2').get_overlap_in_bp_with(self.f), 0)
        self.assertEqual(Feature('chr1', 10, 21, 'ID2').get_overlap_in_bp_with(self.f), 1)
        self.assertEqual(Feature('chr1', 20, 30, 'ID2').get_overlap_in_bp_with(self.f), len(self.f))
        self.assertEqual(Feature('chr2', 10, 21, 'ID2').get_overlap_in_bp_with(self.f), 0)


class FeatureListTestCase(unittest.TestCase):
    def setUp(self):
        self.f = FeatureList([Feature('chr1', 20, 30, 'ID1'),
                              Feature('chr1', 35, 40, 'ID2'),
                              Feature('chr1', 45, 50, 'ID3')])

    def test_len(self):
        self.assertEqual(len(self.f), 3)

    def test_from_data(self):
        l = FeatureList.from_string("chr1 20 30 ID1\nchr1 35 40 ID2")
        self.assertEqual(len(l), 2)
        f1 = l[0]
        self.assertEqual(f1.chromosome, 'chr1')
        self.assertEqual(f1.interval, (20, 30))
        self.assertEqual(f1.name, 'ID1')
        f2 = l[1]
        self.assertEqual(f2.chromosome, 'chr1')
        self.assertEqual(f2.interval, (35, 40))
        self.assertEqual(f2.name, 'ID2')

    def test_filter_by_name(self):
        l = self.f.filter_by_name(['ID2'])
        self.assertEqual(len(l), 1)
        f1 = l[0]
        self.assertEqual(f1.chromosome, 'chr1')
        self.assertEqual(f1.interval, (35, 40))
        self.assertEqual(f1.name, 'ID2')
        l = self.f.filter_by_name(['ID2', 'ID1'])
        self.assertEqual(len(l), 2)

    def test_find_overlap_with_feature(self):
        l = self.f.find_overlap_with_feature(Feature('chr1', 29, 36, 'ID4'))
        self.assertEqual(len(l), 2)
        f1 = l[0]
        self.assertEqual(f1.chromosome, 'chr1')
        self.assertEqual(f1.interval, (20, 30))
        self.assertEqual(f1.name, 'ID1')
        f2 = l[1]
        self.assertEqual(f2.chromosome, 'chr1')
        self.assertEqual(f2.interval, (35, 40))
        self.assertEqual(f2.name, 'ID2')
        l = self.f.find_overlap_with_feature(Feature('chr1', 29, 36, 'ID3'), fraction=1.0)
        self.assertEqual(len(l), 0)
        l = self.f.find_overlap_with_feature(Feature('chr1', 35, 40, 'ID5'), fraction=1.0)
        self.assertEqual(len(l), 1)
        self.assertEqual(l[0].name, 'ID2')

    def test_getitem(self):
        f = self.f['ID1']
        self.assertEqual(f.name, 'ID1')
        f = self.f[0]
        self.assertEqual(f.name, 'ID1')

    def test_iterator(self):
        i = self.f.__iter__()
        f = i.next()
        self.assertEqual(f.name, 'ID1')
        f = i.next()
        self.assertEqual(f.name, 'ID2')
        f = i.next()
        self.assertEqual(f.name, 'ID3')
        try:
            i.next()
            self.assertRaises()
        except StopIteration:
            pass


if __name__ == '__main__':
    unittest.main()
