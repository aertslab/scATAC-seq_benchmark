import timeit
import base64
import sys
import os
import getopt
from itertools import imap
from ConfigParser import RawConfigParser

import MySQLdb

from bunch import Bunch
from cistargetx.iregulon.daemon import load_nomenclature2code


DATABASE_CONFIGURATION_FILE = os.path.join(os.path.dirname(__file__), 'seq-srv-01.ini')

FIND_TARGETOME_STATEMENT = """CALL findTargetome(%s, %s, %s, %s, %s, %s)"""
FIND_MOTIFOME_STATEMENT = """CALL findTargetomeForMotif(%s, %s, %s, %s)"""
FIND_REGULOME_STATEMENT = """CALL findRegulome(%s, %s, %s, %s, %s, %s)"""
FIND_COREGULOME_STATEMENT = """CALL findCoregulome(%s, %s, %s, %s, %s, %s)"""

SIGNATURES_QUERY = """
   SELECT s.name, s.description, ef.AUCValue, ef.NEScore
   FROM enrichedFeature AS ef, targetGene AS g, analysis AS a, geneSignature AS s, geneSignatureCollection AS c
   WHERE ef.featureName = %s AND ef.NEScore >= %s
   AND g.enrichedFeatureID = ef.ID AND g.geneName = %s
   AND ef.analysisID = a.ID AND a.signatureID = s.ID AND s.nomenclatureCode = %s
   AND s.collectionID = c.ID AND c.name = %s
   ORDER BY ef.NEScore DESC;"""

QUERY_TYPE2STATEMENT = {
    'TARGETOME': FIND_TARGETOME_STATEMENT,
    'REGULOME': FIND_REGULOME_STATEMENT,
    'COREGULOME': FIND_COREGULOME_STATEMENT,
    'MOTIFOME': FIND_MOTIFOME_STATEMENT,
    "SIGNATURES": SIGNATURES_QUERY
}

FIND_REGULOME_DYNAMIC_QUERY = """
  SELECT fa.geneName AS transcriptionFactor, tg.geneName,
         count(distinct ef.analysisID) AS occurenceCount, avg(distinct ef.NEScore) AS averageNEScore
  FROM targetGene AS tg, enrichedFeature AS ef, analysis AS a,
       geneSignature AS gs, geneSignatureCollection AS gsc, {0:s}, featureAnnotation AS fa
  WHERE tg.geneName = %s AND tg.enrichedFeatureID = ef.ID AND ef.NEScore >= %s
  AND ef.analysisID = a.ID AND a.signatureID = gs.ID AND gsc.name LIKE %s AND gs.nomenclatureCode = %s
  AND {1:s}
  AND ef.featureName = fa.featureName AND fa.nomenclatureCode = %s
  AND fa.orthologousIdentity >= %s AND fa.motifSimilarityQValue <= %s
  GROUP BY fa.geneName, tg.geneName
  ORDER BY count(distinct ef.analysisID) DESC;"""

FIND_TARGETOME_DYNAMIC_QUERY = """
   SELECT fa.geneName AS transcriptionFactor, tg.geneName AS targetGene,
          count(distinct ef.analysisID) AS occurenceCount, avg(distinct ef.NEScore) AS averageNEScore
   FROM featureAnnotation AS fa, enrichedFeature AS ef, targetGene AS tg, analysis AS a, geneSignature AS gs, geneSignatureCollection AS gsc, {0:s}
   WHERE fa.geneName = %s AND fa.nomenclatureCode = %s
   AND fa.orthologousIdentity >= %s AND fa.motifSimilarityQValue <= %s
   AND fa.featureName = ef.featureName AND ef.NEScore >= %s
   AND tg.isCandidateTarget AND tg.enrichedFeatureID = ef.ID AND ef.analysisID = a.ID AND a.signatureID = gs.ID
   AND gsc.name LIKE %s AND gs.nomenclatureCode = %s
   AND {1:s}
   GROUP BY fa.geneName, tg.geneName
   ORDER BY count(distinct ef.analysisID) DESC;"""

FIND_COREGULOME_DYNAMIC_QUERY = """
   SELECT fa1.geneName AS transcriptionFactor, fa2.geneName AS coFactor, count(distinct a.ID) AS occurenceCount, avg(distinct ef2.NEScore) AS averageNEScore
   FROM featureAnnotation AS fa1, enrichedFeature AS ef1, analysis AS a, geneSignature AS gs, geneSignatureCollection AS gsc, {0:s},
        enrichedFeature AS ef2,featureAnnotation AS fa2
   WHERE fa1.geneName = %s AND fa1.nomenclatureCode = %s
   AND fa1.orthologousIdentity >= %s AND fa1.motifSimilarityQValue <= %s
   AND fa1.featureName = ef1.featureName AND ef1.NEScore >= %s
   AND ef1.analysisID = a.ID AND a.signatureID = gs.ID
   AND gsc.name LIKE %s AND gs.nomenclatureCode = %s
   AND {1:s}
   AND a.ID = ef2.analysisID AND ef2.NEScore >= %s AND ef2.featureName = fa2.featureName AND fa2.nomenclatureCode = %s
   AND fa2.orthologousIdentity >= %s AND fa2.motifSimilarityQValue <= %s AND fa2.geneName != %s
   GROUP BY fa1.geneName, fa2.geneName
   ORDER BY count(distinct ef1.analysisID) DESC;"""

FIND_MOTIFOME_DYNAMIC_QUERY = """
   SELECT ef.featureName as motifID, tg.geneName AS targetGene,
          count(distinct ef.analysisID) AS occurenceCount, avg(distinct ef.NEScore) AS averageNEScore
   FROM enrichedFeature AS ef, targetGene AS tg, analysis AS a, geneSignature AS gs, geneSignatureCollection AS gsc, {0:s},
   WHERE ef.featureName = %s AND ef.NEScore >= %s
   AND tg.isCandidateTarget AND tg.enrichedFeatureID = ef.ID AND ef.analysisID = a.ID AND a.signatureID = gs.ID
   AND gsc.name LIKE %s AND gs.nomenclatureCode = %s
   AND {1:s}
   GROUP BY ef.featureName, tg.geneName
   ORDER BY count(distinct ef.analysisID) DESC;"""


def display_usage(output_fh=sys.stderr, command_name=os.path.basename(sys.argv[0])):
    print >> output_fh, 'Usage: python {0:s} <HGNC|MGI> <TARGETOME|MOTIFOME|REGULOME|COREGULOME> <tfName> [<attribute1>="<value>";<attribute2>="<value>";...]'.format(
        command_name)
    print >> output_fh, 'Usage: python {0:s} <HGNC|MGI> <SIGNATURES> <motifID> <geneID> <sourceName>'.format(
        command_name)
    print >> output_fh, "Options: 1. --cfg=<database_configuration_file>            : the database configuration file (default MySQL database cistargetdb on seq-srv-01)"
    print >> output_fh, "         2. --nes=<minNES> or -n <minNES>                  : minimum NES for motifs to be enriched (default 2.5)"
    print >> output_fh, "         3. --motifsim=<maxMotifSimilarityFDR> or -m <fdr> : maximum motif similarity FDR for TF annotation of enriched motifs (default 0.05)"
    print >> output_fh, "         4. --orthoid=<minOrthologousIdentity> or -o <id>  : minimum orthologous identity fraction for TF annotation of enriched motifs (default 0.0)"
    print >> output_fh, "Possible attributes:"
    print >> output_fh, "         1. sourceName: e.g. sourceName=\"GeneSigDB\""
    print >> output_fh, "         2. meshTerm: e.g. meshTerm=\"Neoplasm\""
    print >> output_fh, "         3. disease: e.g. disease=\"Cancer\""
    print >> output_fh, "E.g.: {0:s} HGNC TARGETOME TP53 'sourceName=\"GeneSigDB\";disease=\"cancer\"'".format(
        command_name)
    print >> output_fh, "E.g.: {0:s} HGNC SIGNATURES ".format(command_name)


def parse_arguments(input_arguments):
    opts, args = getopt.getopt(input_arguments, 'm:n:o:', ['cfg=', 'motifsim=', 'nes=', 'orthoid='])
    if len(args) != 3 and len(args) != 4 and len(args) != 5:
        print >> sys.stderr, "Wrong number of input arguments."
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(1)

    gene_nomenclature = args[0].upper()
    if gene_nomenclature not in ('MGI', 'HGNC'):
        print >> sys.stderr, "'{0:s}' is an invalid gene nomenclature.".format(gene_nomenclature)
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(1)

    query_type = args[1].upper()
    if query_type not in QUERY_TYPE2STATEMENT:
        print >> sys.stderr, "'{0:s}' is an invalid query type.".format(query_type)
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(1)

    nes_threshold = 2.5
    min_orthologous_identity = 0.0
    max_motif_similarity_fdr = 0.05
    cfg_filename = DATABASE_CONFIGURATION_FILE
    for o, a in opts:
        if o in ('-n', '--nes'):
            nes_threshold = float(a)
        elif o == '--cfg':
            cfg_filename = a
        elif o in ('-m', '--motifsim'):
            max_motif_similarity_fdr = float(a)
        elif o in ('-o', '--orthoid'):
            min_orthologous_identity = float(a)
        else:
            raise RuntimeError('Invalid argument supplied.')

    if max_motif_similarity_fdr < 0 or max_motif_similarity_fdr > 1:
        print >> sys.stderr, "The maximum motif similarity FDR must be a real number between [0,1].".format(query_type)
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(1)

    if min_orthologous_identity < 0 or min_orthologous_identity > 1:
        print >> sys.stderr, "The minimum orthologous identity must be a real number between [0,1].".format(query_type)
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(1)

    if not os.path.isfile(cfg_filename):
        print >> sys.stderr, "'{0:s}' file doesn't exist.".format(cfg_filename)
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(1)
    cfg = RawConfigParser()
    cfg.read(cfg_filename)

    if query_type == 'SIGNATURES':
        if len(args) != 5:
            print >> sys.stderr, "Wrong number of input arguments."
            print >> sys.stderr
            display_usage(sys.stderr)
            sys.exit(1)

        motif_id, target_gene_id, source_name = args[2:]

        return cfg, Bunch(gene_nomenclature=gene_nomenclature,
                          motif_id=motif_id,
                          target_gene_id=target_gene_id,
                          query_type=query_type,
                          nes_threshold=nes_threshold,
                          min_orthologous_identity=min_orthologous_identity,
                          max_motif_similarity_fdr=max_motif_similarity_fdr,
                          source_name=source_name,
                          attributes=dict())
    else:
        gene_id = args[2]
        if not gene_id.strip():
            print >> sys.stderr, "Gene identifier cannot be empty.".format(query_type)
            print >> sys.stderr
            display_usage(sys.stderr)
            sys.exit(1)

        if len(args) == 4:
            elements = args[3].split(';')

            def remove_quotes(element):
                key, value = element
                value = value.strip()
                bgn_idx = 1 if value[0] == '"' else 0
                end_idx = -1 if value[-1] == '"' else sys.maxint
                return key.strip(), value[bgn_idx:end_idx]

            attributes = dict(imap(remove_quotes, (element.split('=') for element in elements)))
        elif len(args) == 5:
            print >> sys.stderr, "Gene identifier cannot be empty.".format(query_type)
            print >> sys.stderr
            display_usage(sys.stderr)
            sys.exit(1)
        else:
            attributes = dict()
        source_name = attributes.get('sourceName', '%')
        if 'sourceName' in attributes: del attributes['sourceName']

        return cfg, Bunch(gene_nomenclature=gene_nomenclature, gene_id=gene_id,
                          query_type=query_type,
                          nes_threshold=nes_threshold,
                          min_orthologous_identity=min_orthologous_identity,
                          max_motif_similarity_fdr=max_motif_similarity_fdr,
                          source_name=source_name,
                          attributes=attributes)


def create_db_connection(cfg):
    connection = MySQLdb.connect(host=cfg.get('dbserver', 'servername'),
                                 port=cfg.getint('dbserver', 'port'),
                                 user=cfg.get('dbserver', 'username'),
                                 passwd=base64.b64decode(cfg.get('dbserver', 'password')),
                                 db=cfg.get('dbserver', 'database'))

    # Automatically reconnect to the MySQL server if the connection goes away.
    connection.ping(True)

    return connection


def assemble_query(p, connection, nomenclature2code):
    nomenclature_code = nomenclature2code[p.gene_nomenclature]
    if p.query_type == 'SIGNATURES':
        return QUERY_TYPE2STATEMENT[p.query_type] % connection.literal(
            (p.motif_id, p.nes_threshold, p.target_gene_id, nomenclature_code,
             p.source_name))
    elif not p.attributes:
        if p.query_type == 'MOTIFOME':
            return QUERY_TYPE2STATEMENT[p.query_type] % connection.literal((p.gene_id,
                                                                            nomenclature_code,
                                                                            p.source_name,
                                                                            p.nes_threshold))
        else:
            return QUERY_TYPE2STATEMENT[p.query_type] % connection.literal((p.gene_id,
                                                                            nomenclature_code,
                                                                            p.source_name,
                                                                            p.min_orthologous_identity,
                                                                            p.max_motif_similarity_fdr,
                                                                            p.nes_threshold))
    else:
        if p.query_type == 'REGULOME':
            query = FIND_REGULOME_DYNAMIC_QUERY % connection.literal((p.gene_id,
                                                                      p.nes_threshold,
                                                                      p.source_name,
                                                                      nomenclature_code,
                                                                      nomenclature_code,
                                                                      p.min_orthologous_identity,
                                                                      p.max_motif_similarity_fdr))
        elif p.query_type == 'MOTIFOME':
            query = FIND_MOTIFOME_DYNAMIC_QUERY % connection.literal((p.gene_id,
                                                                      p.nes_threshold,
                                                                      p.source_name,
                                                                      nomenclature_code))
        elif p.query_type == 'TARGETOME':
            query = FIND_TARGETOME_DYNAMIC_QUERY % connection.literal((p.gene_id,
                                                                       nomenclature_code,
                                                                       p.min_orthologous_identity,
                                                                       p.max_motif_similarity_fdr,
                                                                       p.nes_threshold,
                                                                       p.source_name,
                                                                       nomenclature_code))
        else:
            query = FIND_COREGULOME_DYNAMIC_QUERY % connection.literal((p.gene_id,
                                                                        nomenclature_code,
                                                                        p.min_orthologous_identity,
                                                                        p.max_motif_similarity_fdr,
                                                                        p.nes_threshold,
                                                                        p.source_name,
                                                                        nomenclature_code,
                                                                        p.nes_threshold,
                                                                        nomenclature_code,
                                                                        p.min_orthologous_identity,
                                                                        p.max_motif_similarity_fdr,
                                                                        p.gene_id))

        tables_statement = ",".join("geneSignatureAttribute AS gsa{0:d}".format(idx)
                                    for idx in range(len(p.attributes))
                                    )

        attribute_statement = "AND".join(
            "gsa{0:d}.signatureID = gs.ID AND gsa{0:d}.name = '{1:s}' AND gsa{0:d}.textualValue = '{2:}'".format(idx,
                                                                                                                 key,
                                                                                                                 value)
            for idx, (key, value) in enumerate(p.attributes.iteritems())
        )

        return query.format(tables_statement, attribute_statement)


def main():
    start_time = timeit.default_timer()

    cfg, p = parse_arguments(sys.argv[1:])

    connection = create_db_connection(cfg)
    nomenclature2code = load_nomenclature2code(connection)

    cursor = connection.cursor()
    cursor.execute(assemble_query(p, connection, nomenclature2code))
    if p.query_type == 'SIGNATURES':
        print "#signatureName\tsignatureDescription\tAUCValue\tNEScore"
        for sig_name, sig_description, auc_value, nes in cursor:
            print "{0:s}\t{1:s}\t{2:f}\t{3:f}".format(sig_name, sig_description, auc_value, nes)
    else:
        if p.query_type == 'REGULOME':
            print "#geneName\ttfName\toccurenceCount\taverageNEScore"
        elif p.query_type == 'COREGULOME':
            print "#tfName\tcoFactorName\toccurenceCount\taverageNEScore"
        else:
            print "#tfName\ttargetGene\toccurenceCount\taverageNEScore"
        for tf_name, target_gene, occurence_count, average_nes in cursor:
            print "{0:s}\t{1:s}\t{2:d}\t{3:f}".format(tf_name, target_gene, occurence_count, average_nes)

    cursor.close()
    connection.close()

    print >> sys.stderr, "Elapsed time = {0:f}s".format(timeit.default_timer() - start_time)


if __name__ == "__main__":
    main()
