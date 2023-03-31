import os
import re
import sys
from ConfigParser import RawConfigParser
from base64 import b64decode
from collections import namedtuple

from MySQLdb import connect


INSERT_METADATA_STATEMENT = r"""
  INSERT INTO metatargetomeMetaData(sourceName, tfName, speciesNomenclatureCode)
  VALUES(%s, %s, %s);
"""

INSERT_DATA_STATEMENT = r"""
  INSERT INTO metatargetome(metaID, geneName, occurenceCount)
  VALUES(%s, %s, %s);
"""


def create_db_connection(cfg):
    connection = connect(host=cfg.get('dbserver', 'servername'),
                         port=cfg.getint('dbserver', 'port'),
                         user=cfg.get('dbserver', 'username'),
                         passwd=b64decode(cfg.get('dbserver', 'password')),
                         db=cfg.get('dbserver', 'database'))

    # Automatically reconnect to the MySQL server if the connection goes away.
    connection.ping(True)

    return connection


def load_metatargetomes(filename):
    Metatargetome = namedtuple('Metatargetome', 'tf targetome')
    metatargetomes = []

    def compare(t1, t2):
        r = cmp(t2[1], t1[1])

        if r:
            return r

        return cmp(t1[0], t2[0])

    with open(filename, 'r') as input_fh:
        cur_tf, targetome = None, []

        for line in input_fh:
            if not line.strip():
                continue

            if line.startswith('#tf'):
                if cur_tf:
                    targetome_as_tuple = tuple(
                        sorted([(target_gene, int(occurence_count)) for target_gene, occurence_count in targetome],
                               cmp=compare)
                    )

                    metatargetomes.append(Metatargetome(cur_tf, targetome_as_tuple))
                    cur_tf, targetome = None, []

                cur_tf = line.rstrip().split("=")[1]
            elif not line.startswith('#'):
                columns = re.split('[ \t]+', line.rstrip())
                gene_id, rank = columns[0], int(columns[1])
                targetome.append((gene_id, rank))

    if cur_tf:
        targetome_as_tuple = tuple(
            sorted([(target_gene, int(occurence_count)) for target_gene, occurence_count in targetome],
                   cmp=compare)
        )

        metatargetomes.append(Metatargetome(cur_tf, targetome_as_tuple))

    return metatargetomes


def insert_metatargetome(connection, metatargetome, sourcename, species_nomenclature_code):
    cursor = connection.cursor()

    cursor.execute(INSERT_METADATA_STATEMENT, (sourcename, metatargetome.tf, species_nomenclature_code))
    meta_id = connection.insert_id()

    for target_gene, occurence_count in metatargetome.targetome:
        cursor.execute(INSERT_DATA_STATEMENT, (meta_id, target_gene, occurence_count))

    cursor.close()


def main():
    ini_filename, upload_filename, species_nomenclature = sys.argv[1:]
    cfg = RawConfigParser()
    cfg.read(ini_filename)

    metatargetomes = load_metatargetomes(upload_filename)
    species_nomenclature_code = int(species_nomenclature)
    sourcename = os.path.basename(upload_filename).split('.')[0]

    connection = create_db_connection(cfg)

    for metatargetome in metatargetomes:
        insert_metatargetome(connection, metatargetome, sourcename, species_nomenclature_code)
        connection.commit()

    connection.close()


if __name__ == "__main__":
    main()
