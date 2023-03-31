import unittest, os, numpy
from cistargetx.common.orderstatistics import combine_via_order_statistics


TEST_FILENAMES = [
    os.path.join(os.path.dirname(__file__), 'test1.rr'),
    os.path.join(os.path.dirname(__file__), 'test2.rr'),
    os.path.join(os.path.dirname(__file__), 'test3.rr'),
    os.path.join(os.path.dirname(__file__), 'test4.rr'),
    os.path.join(os.path.dirname(__file__), 'test5.rr')]
RESULT_FILENAME = os.path.join(os.path.dirname(__file__), 'result.r')
MACHINE_EPS = 1E-15


def _load_rr_file(filename):
    record_structure = {'names': ('id', 'rr'), 'formats':('S255', 'float64')}
    return numpy.loadtxt(filename, dtype=record_structure, comments=';')


def _append_column(rrs, filename, idx):
    table = _load_rr_file(filename)
    rrs[:, idx] = table['rr'][numpy.argsort(table['id'])]


def load_relative_rankings(filenames):
    table = _load_rr_file(filenames[0])
    rrs = numpy.empty(shape=(len(table['id']), len(filenames)) , dtype=numpy.float64)
    sorted_idx = numpy.argsort(table['id'])
    rrs[:, 0] = table['rr'][sorted_idx]
    gene_ids = table['id'][sorted_idx]
    for idx in range(1, len(filenames)):
        _append_column(rrs, filenames[idx], idx)
    return gene_ids, rrs


def millis():
    import time as time_ #make sure we don't override time
    return int(round(time_.time() * 1000))


class OrderStatisticsTestCase(unittest.TestCase):
    def setUp(self):
        self.gene_ids, self.rrs = load_relative_rankings(TEST_FILENAMES)
        ids, self.result = load_relative_rankings([RESULT_FILENAME])


    def test_combine_via_order_statistics(self):
        millis_before = millis()
        r = combine_via_order_statistics(self.rrs)
        millis_after = millis()
        print "Execution time = {0:d}ms".format(millis_after-millis_before)
        self.assertTrue(numpy.where(self.result.flat - r >= MACHINE_EPS)[0].size == 0)


if __name__ == '__main__':
    unittest.main()
