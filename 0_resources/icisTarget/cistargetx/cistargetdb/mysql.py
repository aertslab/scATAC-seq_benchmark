import base64
import itertools
import logging

import MySQLdb
import _mysql_exceptions


INSERT_COLLECTION_METADATA_STATEMENT = """INSERT INTO geneSignatureCollection(name, description, version, sourceFile)
                                          VALUES (%s, %s, %s, %s)"""
GET_COLLECTION_ID_STATEMENT = """SELECT ID FROM geneSignatureCollection WHERE name = %s AND version = %s"""
REMOVE_COLLECTION_STATEMENT = """DELETE FROM geneSignatureCollection WHERE name = %s AND version = %s;"""

QUERY_GENE_SIGNATURE_METADATA_STATEMENT = """SELECT ID FROM geneSignature WHERE name = %s"""
INSERT_GENE_SIGNATURE_METADATA_STATEMENT = """INSERT INTO geneSignature(name, description, collectionID, nomenclatureCode, pubMedIdentifier)
                                          VALUES (%s, %s, %s, %s, %s)"""
INSERT_GENE_SIGNATURE_GENE_STATEMENT = """INSERT INTO geneSignatureGene(signatureID, geneName) VALUES (%s, %s)"""

QUERY_GENE_SIGNATURE_STATEMENT = """SELECT gs.ID FROM geneSignature AS gs, geneSignatureCollection as col
                                     WHERE gs.pubMedIdentifier = %s AND gs.collectionID = col.ID AND (col.name <> %s OR col.version <> %s)"""
REMOVE_GENE_SIGNATURES_STATEMENT = """DELETE FROM geneSignature AS gs WHERE gs.pubMedIdentifier = %s
                                      AND NOT EXISTS (SELECT * FROM geneSignatureCollection AS col WHERE col.ID = gs.collectionID AND col.name = %s AND col.version = %s);"""

QUERY_GENE_SIGNATURE_NUMERIC_ATTRIBUTE_STATEMENT = """SELECT * FROM geneSignatureAttribute WHERE signatureID = %s AND name = %s AND numericValue = %s"""
QUERY_GENE_SIGNATURE_TEXTUAL_ATTRIBUTE_STATEMENT = """SELECT * FROM geneSignatureAttribute WHERE signatureID = %s AND name = %s AND textualValue = %s"""
INSERT_GENE_SIGNATURE_NUMERIC_ATTRIBUTE_STATEMENT = """INSERT INTO geneSignatureAttribute(signatureID, name, numericValue) VALUES (%s, %s, %s)"""
INSERT_GENE_SIGNATURE_TEXTUAL_ATTRIBUTE_STATEMENT = """INSERT INTO geneSignatureAttribute(signatureID, name, textualValue) VALUES (%s, %s, %s)"""

INSERT_ANALYSIS_STATEMENT = """INSERT INTO analysis(signatureID, calculatedOn, configurationFile, databaseName, rankThreshold, AUCThreshold, NESThreshold, fractionOfMappedGeneIDs)
                               VALUES (%s, now(), %s, %s, %s, %s, %s, %s)"""

INSERT_ENRICHED_FEATURE_STATEMENT = """INSERT INTO enrichedFeature(analysisID, featureName, rank, AUCValue, NEScore)
                                       VALUES (%s, %s, %s, %s, %s)"""

INSERT_TARGET_GENE_STATEMENT = """INSERT INTO targetGene(enrichedFeatureID, geneName, rank, isCandidateTarget)
                                  VALUES (%s, %s, %s, %s)"""

ANALYSIS_COUNT_STATEMENT = """SELECT COUNT(*) FROM analysis WHERE signatureID = %s"""


def create_db_connection(cfg, new=False):
    if new:
        connection = MySQLdb.connect(host=cfg.get('dbserver', 'servername'),
                                     port=cfg.getint('dbserver', 'port'),
                                     user=cfg.get('dbserver', 'username'),
                                     passwd=base64.b64decode(cfg.get('dbserver', 'password')))
    else:
        connection = MySQLdb.connect(host=cfg.get('dbserver', 'servername'),
                                     port=cfg.getint('dbserver', 'port'),
                                     user=cfg.get('dbserver', 'username'),
                                     passwd=base64.b64decode(cfg.get('dbserver', 'password')),
                                     db=cfg.get('dbserver', 'database'))

    # Automatically reconnect to the MySQL server if the connection goes away.
    connection.ping(True)

    return connection


class CollectionAlreadyPresentException(Exception):
    def __init__(self, collection_id):
        super(CollectionAlreadyPresentException, self).__init__()
        self.collection_id = collection_id


def insert_collection_metadata(connection, metadata):
    cursor = connection.cursor()
    try:
        cursor.execute(INSERT_COLLECTION_METADATA_STATEMENT,
                       (metadata.name, metadata.description, metadata.version, metadata.sourceFile))
        connection.commit()
    except _mysql_exceptions.IntegrityError:
        cursor.execute(GET_COLLECTION_ID_STATEMENT, (metadata.name, metadata.version))
        collection_id = cursor.fetchone()[0]
        raise CollectionAlreadyPresentException(collection_id)
    collection_id = connection.insert_id()
    cursor.close()
    return collection_id


def remove_collection(connection, metadata):
    cursor = connection.cursor()
    cursor.execute(REMOVE_COLLECTION_STATEMENT, (metadata.name, metadata.version))
    connection.commit()
    logging.warning('{0:d} collections where deleted from database.'.format(cursor.rowcount))
    cursor.close()


def insert_gene_signature(connection, collection_id, metadata, signature, nomenclature2code, duplicate_strategy):
    cursor = connection.cursor()

    pubmed_id = signature.pubmedID
    if duplicate_strategy == 'rm_new':
        cursor.execute(QUERY_GENE_SIGNATURE_STATEMENT, (pubmed_id, metadata.name, metadata.version))
        if cursor.fetchone():
            return None, False
    elif duplicate_strategy == 'rm_old':
        cursor.execute(REMOVE_GENE_SIGNATURES_STATEMENT, (pubmed_id, metadata.name, metadata.version))
        connection.commit()
        logging.warning('{0:d} signatures where deleted from database because they share PubMedID with {0:s}.'.format(
            cursor.rowcount, signature.name))

    cursor.execute(QUERY_GENE_SIGNATURE_METADATA_STATEMENT, (signature.name,))
    row = cursor.fetchone()
    if row:
        signature_id = row[0]
        new = False
    else:
        cursor.execute(INSERT_GENE_SIGNATURE_METADATA_STATEMENT, (signature.name,
                                                                  signature.description,
                                                                  collection_id,
                                                                  nomenclature2code[signature.nomenclature],
                                                                  signature.pubmedID))
        signature_id = connection.insert_id()
        cursor.executemany(INSERT_GENE_SIGNATURE_GENE_STATEMENT,
                           zip(itertools.repeat(signature_id), signature.identifiers))
        connection.commit()
        new = True

    for name, values in signature.attributes.iteritems():
        for value in values:
            numeric_value = isinstance(value, float)
            if numeric_value:
                cursor.execute(QUERY_GENE_SIGNATURE_NUMERIC_ATTRIBUTE_STATEMENT, (signature_id, name, value))
            else:
                cursor.execute(QUERY_GENE_SIGNATURE_TEXTUAL_ATTRIBUTE_STATEMENT, (signature_id, name, value))
            row = cursor.fetchone()
            if not row:
                if numeric_value:
                    cursor.execute(INSERT_GENE_SIGNATURE_NUMERIC_ATTRIBUTE_STATEMENT, (signature_id, name, value))
                else:
                    cursor.execute(INSERT_GENE_SIGNATURE_TEXTUAL_ATTRIBUTE_STATEMENT, (signature_id, name, value))

    connection.commit()

    cursor.execute(ANALYSIS_COUNT_STATEMENT, (signature_id,))
    n_analyses = cursor.fetchone()[0]
    analysis_present = n_analyses > 0

    cursor.close()

    return signature_id, new, analysis_present


def insert_analysis(connection, signature_id, signature, parameters, curves):
    cursor = connection.cursor()
    for curve in curves:
        cursor.execute(INSERT_ANALYSIS_STATEMENT, (signature_id,
                                                   parameters.inifile,
                                                   curve.database.name,
                                                   parameters.rankthreshold,
                                                   parameters.aucthreshold,
                                                   parameters.nes_threshold,
                                                   curve.fraction_of_mapped_gene_ids(signature.identifiers)))
        analysis_id = connection.insert_id()
        for idx, (feature_id, auc, nes, feature) in enumerate(curve.iterate_enriched_features()):
            cursor.execute(INSERT_ENRICHED_FEATURE_STATEMENT, (analysis_id, feature_id, idx + 1, auc, nes))
            enriched_feature_id = connection.insert_id()
            candidate_target_threshold = curve.get_critical_point(feature_id)[0]
            for rank, geneName, known in feature.get_top_ranked_gene_ids(feature_id):
                is_candidate_target = rank <= candidate_target_threshold
                cursor.execute(INSERT_TARGET_GENE_STATEMENT,
                               (enriched_feature_id, geneName, rank, int(is_candidate_target)))

    connection.commit()

    cursor.close()
