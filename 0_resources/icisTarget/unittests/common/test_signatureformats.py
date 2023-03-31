import unittest, os
from cistargetx.common.signatureformats import GeneSignatureFileFormat, GMT_COMMA, LEGACY, GeneSignature


TEST_GMT_FILENAME = os.path.join(os.path.dirname(__file__), 'benchmark.hgnc.gmt')
TEST_LEGACY_FILENAME = os.path.join(os.path.dirname(__file__), 'CREB_Zhang.hgnc')


class FormatCase(unittest.TestCase):
    def test_guess_format(self):
        with open(TEST_GMT_FILENAME, 'r') as input_fh:
            gene_signature_file_format = GeneSignatureFileFormat.guess_file_format(input_fh.readlines())
            self.assertIsNotNone(gene_signature_file_format)
            self.assertEqual(GMT_COMMA, gene_signature_file_format)
        with open(TEST_LEGACY_FILENAME, 'r') as input_fh:
            gene_signature_file_format = GeneSignatureFileFormat.guess_file_format(input_fh.readlines())
            self.assertIsNotNone(gene_signature_file_format)
            self.assertEqual(LEGACY, gene_signature_file_format)


class GeneSignatureCase(unittest.TestCase):
    def test_load_from_file(self):
        sigs = GeneSignature.load_from_file(TEST_GMT_FILENAME)
        self.assertEqual(len(sigs), 20)
        sigs = GeneSignature.load_from_file(TEST_LEGACY_FILENAME)
        self.assertEqual(len(sigs), 1)


if __name__ == '__main__':
    unittest.main()
