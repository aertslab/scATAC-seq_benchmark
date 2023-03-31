import os
import re

# GMT File format:
#   http://www.broadinstitute.org/cancer/software/genepattern/gp_guides/file-formats/sections/gmt
#
# Line format:
#   (gene set name) (tab) (description) (tab) (gene 1) (tab) (gene 2) (tab) ... (gene n)
#
#   - column 1:   gene set name. Duplicate names are not allowed.
#   - column 2:   gene set description.
#   - column 3-n: genes in the gene set.


class GeneSignatureFileFormat(object):
    @staticmethod
    def create_gmt_format(separator):
        def _parse(lines):
            for line in lines:
                if line.startswith("#") or not line.strip():
                    continue

                columns = re.split('\t', line.rstrip())

                if separator == '\t':
                    yield columns[0], columns[1], columns[2:]
                else:
                    yield columns[0], columns[1], columns[2].split(separator)

        return GeneSignatureFileFormat("GMT", separator, _parse)

    @staticmethod
    def guess_file_format(lines):
        data_lines = filter(lambda l: not l.startswith("#") and l.strip(), lines)

        def count_columns(line):
            return line.count('\t')

        column_counts = map(count_columns, data_lines)

        def create_predicate_equal_to(ref_value):
            return lambda value: value == ref_value

        if all(map(create_predicate_equal_to(0), column_counts)):
            return LEGACY
        elif all(map(create_predicate_equal_to(2), column_counts)):
            signature_column = map(lambda l: re.split('\t', l.rstrip())[2], data_lines)

            def create_predicate_contains(c):
                return lambda l: c in l

            if any(map(create_predicate_contains(';'), signature_column)):
                return GMT_SEMICOLON
            # Caveat: a comma might be accidentally be present in a gene name (FlyBase).
            elif any(map(create_predicate_contains(','), signature_column)):
                return GMT_COMMA
            # The case were each signature contains one gene only.
            else:
                return GMT_TAB
        # Counting TABs and not columns.
        elif all(map(lambda c: c >= 2, column_counts)):
            return GMT_TAB
        else:
            return None

    def __init__(self, name, separator, parse_method):
        assert name, "Error: Name for gene signature file format object must be specified."
        self.name = name
        self.separator = separator
        self.parse = parse_method

    def __hash__(self):
        return hash((self.name, self.separator))

    def __eq__(self, other):
        return (self.name, self.separator) == (other.name, other.separator)


LEGACY = GeneSignatureFileFormat('legacy', '', lambda lines: [('?', '?', map(lambda l: l.rstrip(), lines))])
GMT_TAB = GeneSignatureFileFormat.create_gmt_format('\t')
GMT_SEMICOLON = GeneSignatureFileFormat.create_gmt_format(';')
GMT_COMMA = GeneSignatureFileFormat.create_gmt_format(',')


def _strip_blanks(string):
    string = re.sub('[ \t]+$', '', string)
    return re.sub('^[ \t]+', '', string)


def _is_not_empty(string):
    return len(string) > 0


class GeneSignature(object):
    @staticmethod
    def load_from_string(data, name, gene_signature_file_format=None):
        return GeneSignature._load_from_lines(re.split('[\n\r]+', data), name, gene_signature_file_format)

    @staticmethod
    def load_from_stream(stream, name, gene_signature_file_format=None):
        return GeneSignature._load_from_lines(stream.readlines(), name, gene_signature_file_format)

    @staticmethod
    def _load_from_lines(data, name, gene_signature_file_format=None):
        lines = filter(_is_not_empty, map(_strip_blanks, data))

        if not gene_signature_file_format:
            gene_signature_file_format = GeneSignatureFileFormat.guess_file_format(lines)

        assert gene_signature_file_format, "Error: Gene signature file format could not be detected."

        if gene_signature_file_format == LEGACY:
            for dummy_name, dummy_description, gene_ids in gene_signature_file_format.parse(lines):
                return [GeneSignature(name, name, gene_ids)]
        else:
            signatures = []
            prev_name = set()
            for name, description, gene_ids in gene_signature_file_format.parse(lines):
                if name in prev_name:
                    assert False, "Error: ID '{0:s}' for gene signature is not unique.".format(name)
                if len(gene_ids) == 0:
                    assert False, "Error: ID '{0:s}' has no genes in its gene set.".format(name)
                signatures.append(GeneSignature(name, description, gene_ids))
                prev_name.add(name)
            return signatures

    @staticmethod
    def load_from_file(filename, gene_signature_file_format=None,
                       transform=lambda f: os.path.splitext(os.path.basename(f))[0]):
        with open(filename, 'r') as input_fh:
            return GeneSignature.load_from_stream(input_fh, transform(filename), gene_signature_file_format)

    def __init__(self, name, description, identifiers):
        self.name = name
        self.description = description
        self.identifiers = set(identifiers)

    def __len__(self):
        return len(self.identifiers)
