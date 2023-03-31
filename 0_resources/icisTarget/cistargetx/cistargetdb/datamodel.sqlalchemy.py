from sqlalchemy import Boolean, Column, DateTime, Float, ForeignKey, Index, Integer, MetaData, PrimaryKeyConstraint, \
    String, Table, Text, create_engine


# server_metadata = MetaData()
server = create_engine('mysql://cisTargetDB:bdtegratsic@10.35.46.100:3306/cistargetdb')
# TODO
server = create_engine('mysql://bram:bram@10.35.46.100:3306/tf_motif_mapping')
meta = MetaData()
meta.reflect(bind=server)
meta.tables.keys()

conn = server.connect()
conn.execute(meta.tables['motifAnnotation'].select())
conn.close()

metadata = MetaData()

gene_signature_collection_table = Table(
    'geneSignatureCollection',
    metadata,
    Column('ID', Integer, primary_key=True, autoincrement=True),
    Column('name', String(255), nullable=False),
    Column('description', Text),
    Column('version', String(32), nullable=False),
    Column('sourceFile', Text, nullable=False)
)
Index('geneSignatureCollectionIDIndex', gene_signature_collection_table.c.ID, unique=True)
Index('nameVersionIndex', gene_signature_collection_table.c.name, gene_signature_collection_table.c.version,
      unique=True)

gene_signature_table = Table(
    'geneSignature',
    metadata,
    Column('ID', Integer, primary_key=True, autoincrement=True),
    Column('name', String(255), nullable=False),
    Column('description', Text),
    Column('collectionID', None, ForeignKey(gene_signature_collection_table.c.ID, ondelete="CASCADE"), nullable=False),
    Column('nomenclatureCode', Integer, nullable=False),
    Column('pubMedIdentifier', Integer)
)
Index('geneSignatureIDIndex', gene_signature_table.c.ID, unique=True)
Index('collectionIDIndex', gene_signature_table.c.collectionID)

gene_signature_gene_table = Table(
    'geneSignatureGene',
    metadata,
    Column('signatureID', None, ForeignKey(gene_signature_table.c.ID, ondelete="CASCADE"), nullable=False),
    Column('geneName', String(255), nullable=False),
    PrimaryKeyConstraint('signatureID', 'geneName')
)
Index('geneSignatureGeneSignatureIDIndex', gene_signature_gene_table.c.signatureID)

gene_signature_attribute_table = Table(
    'geneSignatureAttribute',
    metadata,
    Column('signatureID', None, ForeignKey(gene_signature_table.c.ID, ondelete="CASCADE"), nullable=False),
    Column('name', String(255), nullable=False),
    Column('textualValue', String(255)),
    Column('numericalValue', Float),
    PrimaryKeyConstraint('signatureID', 'name', 'textualValue', 'numericalValue')
)
Index('geneSignatureAttributeSignatureIDIndex', gene_signature_attribute_table.c.signatureID)
Index('geneSignatureAttributeNameIndex', gene_signature_attribute_table.c.name)

analysis_table = Table(
    'analysis',
    metadata,
    Column('ID', Integer, primary_key=True, autoincrement=True),
    Column('signatureID', None, ForeignKey(gene_signature_table.c.ID), nullable=False),
    Column('calculatedOn', DateTime, nullable=False),
    Column('configurationFile', Text),
    Column('databaseName', String(255)),
    Column('rankThreshold', Float, nullable=False),
    Column('AUCThreshold', Float, nullable=False),
    Column('NESThreshold', Float, nullable=False),
    Column('fractionOfMappedGeneIDs', Float, nullable=False)
)
Index('analysisIDIndex', analysis_table.c.ID, unique=False)
Index('analysisSignatureIDIndex', analysis_table.c.signatureID)

enriched_feature_table = Table(
    'enrichedFeature',
    metadata,
    Column('ID', Integer, primary_key=True, autoincrement=True),
    Column('analysisID', None, ForeignKey(analysis_table.c.ID, ondelete="CASCADE"), nullable=False),
    Column('featureName', String(255), nullable=False),
    Column('rank', Integer, nullable=False),
    Column('AUCValue', Float, nullable=False),
    Column('NEScore', Float, nullable=False)
)
Index('enrichedFeatureIDIndex', enriched_feature_table.c.ID, unique=True)
Index('enrichedFeatureAnalysisIDIndex', enriched_feature_table.c.analysisID)
Index('featureNameIndex', enriched_feature_table.c.featureName)

target_gene_table = Table(
    'targetGene',
    metadata,
    Column('enrichedFeatureID', None, ForeignKey(enriched_feature_table.c.ID), nullable=False),
    Column('geneName', String(255), nullable=False),
    Column('rank', Integer, nullable=False),
    Column('isCandidateTarget', Boolean, nullable=False),
    PrimaryKeyConstraint('enrichedFeatureID', 'geneName')
)
Index('targetGeneEnrichedFeatureIDIndex', target_gene_table.c.enrichedFeatureID)
Index('targetGeneNameIndex', target_gene_table.c.geneName)

engine = create_engine('sqlite:///:memory:', echo=True)
metadata.create_all(engine)
