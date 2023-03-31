import logging
from collections import namedtuple, defaultdict

from cistargetx.common.signatureformats import GeneSignature


CollectionMetaData = namedtuple('CollectionMetaData', 'name description version sourceFile')


class MetadataEnrichedGeneSignature(GeneSignature):
    @staticmethod
    def load_collection_metadata(lines):
        name, description, version = [None] * 3
        for line in lines:
            if not line.startswith('#') or "=" not in line: continue
            columns = line.rstrip()[1:].split("=")
            key = columns[0].strip()
            value = "=".join(columns[1:]).strip()
            if key == 'name':
                name = value
            elif key == 'description':
                description = value
            elif key == 'version':
                version = value
        return CollectionMetaData(name, description, version, "\n".join(lines))


    @staticmethod
    def load_signatures(lines):
        result = list()
        for idx, line in enumerate(lines):
            if line.startswith('#') or not line.strip():
                continue
            columns = line.rstrip().split('\t')
            if len(columns) != 3:
                logging.error('Invalid gene signature on line {0:d}'.format(idx + 1))
            name = columns[0]
            description, nomenclature, pubmed_id = [None] * 3
            attributes = defaultdict(set)
            for attribute in columns[1].split(';'):
                if not attribute or not attribute.strip():
                    continue
                # TODO: Parsing should be done more intelligent: ignore ; and = inside double quotes ...
                try:
                    key, value = attribute.split("=")
                except ValueError:
                    logging.error("Cannot parse attribute '{0:s}' on line {1:d}.".format(attribute, idx + 1))
                key = key.strip()
                text_value = value.strip()[1:-1]
                if key == 'nomenclature':
                    nomenclature = text_value
                elif key == 'description':
                    description = text_value
                elif key == 'pubmedid':
                    pubmed_id = int(text_value)
                else:
                    try:
                        value = float(text_value)
                    except ValueError:
                        value = text_value
                    attributes[key].add(value)
            identifiers = filter(lambda name: name.strip(), columns[2].split(';'))
            result.append(
                MetadataEnrichedGeneSignature(name, description, identifiers, nomenclature, pubmed_id, attributes))
        return result

    def __init__(self, name, description, identifiers, nomenclature, pubmedID=None, attributes=dict()):
        super(MetadataEnrichedGeneSignature, self).__init__(name, description, identifiers)
        self.nomenclature = nomenclature
        self.pubmedID = pubmedID
        self.attributes = attributes
