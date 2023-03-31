import os
from cistargetx.common.lrucache import cache
from featurelist import FeatureList


REGIONS_BED_FILENAME = {
    ('dm3', 'flybase_r5.37'): '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/cut_100_0.1_all_genome_noExons_noInsulators_min500_Dm3.bed',
    ('dm6', 'flybase_r6.02'): '/home/icistarget/data/regions/drosophila_melanogaster/dm6/flybase-r6.02/cut_50_0.1_all_genome_noExons_min500_dm6.bed',
    ('hg19', 'refseq_r45'): '/home/icistarget/data/regions/homo_sapiens/hg19/refseq-r45/hg19__refseq_r45__ClusteredUniformDHS_all_merge_cleaned2_features_rm-insul_rm-exons2_extend.regionid-location.bed',
    ('hg38', 'refseq_r109'): '/home/icistarget/data/regions/homo_sapiens/hg38/refseq-r109/region_based/hg38.ENCODE3_regions.blacklisted.regionid-location.bed',
    ('mm9', 'refseq_r45'): '/home/icistarget/data/regions/mus_musculus/mm9/refseq-r45/mm9__refseq_r45__regulatory_regions_DHS.regionid-location.bed',
    ('mm9', 'refseq_r70'): '/home/icistarget/data/regions/mus_musculus/mm9/refseq-r70/mm9__refseq_r70__regulatory_regions.regionid-location.bed',
    ('mm10', 'refseq_r108'): '/home/icistarget/data/regions/mus_musculus/mm10/refseq-r108/region_based/mm10.ENCODE3_regions.blacklisted.regionid-location.bed'
}


class RegionFeatures:
    REGION_FEATURES = dict()

    @staticmethod
    def initialise(assembly, gene_annotation_version):
        if (assembly, gene_annotation_version) not in RegionFeatures.REGION_FEATURES:
            RegionFeatures.REGION_FEATURES[(assembly, gene_annotation_version)] = RegionFeatures.REGION_FEATURES.get(
                    (assembly, gene_annotation_version),
                    FeatureList.from_bed_file(REGIONS_BED_FILENAME[(assembly, gene_annotation_version)])
            )


class GeneIDType:
    FBGN = 'fbgn' # FlyBase gene number
    SYMBOL = 'symbol' # Gene symbol. Genes that do not have a name yet have a CG (= Computed Gene number) associated with it.
    CG = 'cg' # CG gene number.
    HGNC_SYMBOL = 'hgnc_symbol' # HGNC gene symbol (human).
    MGI_SYMBOL = 'mgi_symbol' # MGI gene symbol (mouse).


class Delineation:
    UPSTREAM5KB_FULL_TX = 'up5kb+fullTx'
    UPSTREAM5KB_5UTR_INTRON1 = 'up5kb+5utr+intron1'
    UPSTREAM5KB_FULL_TX_DOWNSTREAM5KB = 'up5kb+fullTx+down5kb'
    UPSTREAM10KB_FULL_TX_DOWNSTREAM10KB = 'up10kb+fullTx+down10kb'
    UPSTREAM5KB_5UTR_INTRON1_UCSC = 'UCSC-up5kb+5utr+intron1'
    UPSTREAM5KB_ALL_INTRONS = 'up5kb+introns'
    UPSTREAM10KB_TSS_DOWNSTREAM10KB = 'up10kb+TSS+down10kb'


ID_TYPE_DELINEATION_2_BED_FILENAME = {
    ('dm3', 'flybase_r5.37', GeneIDType.FBGN, Delineation.UPSTREAM5KB_FULL_TX) : '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/flybase-dmel-r5.37-limited-upstream5000-full-transcript-fbgn.bed',
    ('dm3', 'flybase_r5.37', GeneIDType.SYMBOL, Delineation.UPSTREAM5KB_FULL_TX) : '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/flybase-dmel-r5.37-limited-upstream5000-full-transcript-symbol.bed',
    ('dm3', 'flybase_r5.37', GeneIDType.CG, Delineation.UPSTREAM5KB_FULL_TX) : '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/flybase-dmel-r5.37-limited-upstream5000-full-transcript-cg.bed',

    ('dm3', 'flybase_r5.37', GeneIDType.FBGN, Delineation.UPSTREAM5KB_5UTR_INTRON1) : '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/flybase-dmel-r5.37-limited-upstream5000-5utr-intron1-fbgn.bed',
    ('dm3', 'flybase_r5.37', GeneIDType.SYMBOL, Delineation.UPSTREAM5KB_5UTR_INTRON1) : '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/flybase-dmel-r5.37-limited-upstream5000-5utr-intron1-symbol.bed',
    ('dm3', 'flybase_r5.37', GeneIDType.CG, Delineation.UPSTREAM5KB_5UTR_INTRON1) : '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/flybase-dmel-r5.37-limited-upstream5000-5utr-intron1-cg.bed',

    ('dm3', 'flybase_r5.37', GeneIDType.FBGN, Delineation.UPSTREAM5KB_FULL_TX_DOWNSTREAM5KB) : '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/flybase-dmel-r5.37-limited-upstream5000-full-transcript-downstream5000-fbgn.bed',
    ('dm3', 'flybase_r5.37', GeneIDType.SYMBOL, Delineation.UPSTREAM5KB_FULL_TX_DOWNSTREAM5KB) : '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/flybase-dmel-r5.37-limited-upstream5000-full-transcript-downstream5000-symbol.bed',
    ('dm3', 'flybase_r5.37', GeneIDType.CG, Delineation.UPSTREAM5KB_FULL_TX_DOWNSTREAM5KB) : '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/flybase-dmel-r5.37-limited-upstream5000-full-transcript-downstream5000-cg.bed',

    ('dm3', 'flybase_r5.37', GeneIDType.FBGN, Delineation.UPSTREAM10KB_FULL_TX_DOWNSTREAM10KB) : '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/flybase-dmel-r5.37-limited-upstream10000-full-transcript-downstream10000-fbgn.bed',
    ('dm3', 'flybase_r5.37', GeneIDType.SYMBOL, Delineation.UPSTREAM10KB_FULL_TX_DOWNSTREAM10KB) : '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/flybase-dmel-r5.37-limited-upstream10000-full-transcript-downstream10000-symbol.bed',
    ('dm3', 'flybase_r5.37', GeneIDType.CG, Delineation.UPSTREAM10KB_FULL_TX_DOWNSTREAM10KB) : '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/flybase-dmel-r5.37-limited-upstream10000-full-transcript-downstream10000-cg.bed',

    ('dm3', 'flybase_r5.37', GeneIDType.CG, Delineation.UPSTREAM5KB_ALL_INTRONS) : '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/legacy-dmel-limited-upstream5000-all-introns-cg.bed',
    ('dm3', 'flybase_r5.37', GeneIDType.SYMBOL, Delineation.UPSTREAM5KB_ALL_INTRONS) : '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/legacy-dmel-limited-upstream5000-all-introns-symbol.bed',
    ('dm3', 'flybase_r5.37', GeneIDType.FBGN, Delineation.UPSTREAM5KB_ALL_INTRONS) : '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/legacy-dmel-limited-upstream5000-all-introns-fbgn.bed',

    ('dm3', 'flybase_r5.37', GeneIDType.CG, Delineation.UPSTREAM5KB_5UTR_INTRON1_UCSC) : '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/legacy-dmel-UCSC-limited-upstream5000-5utr-intron1-cg.bed',
    ('dm3', 'flybase_r5.37', GeneIDType.SYMBOL, Delineation.UPSTREAM5KB_5UTR_INTRON1_UCSC) : '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/legacy-dmel-UCSC-limited-upstream5000-5utr-intron1-symbol.bed',
    ('dm3', 'flybase_r5.37', GeneIDType.FBGN, Delineation.UPSTREAM5KB_5UTR_INTRON1_UCSC) : '/home/icistarget/data/regions/drosophila_melanogaster/dm3/flybase-r5.37/legacy-dmel-UCSC-limited-upstream5000-5utr-intron1-fbgn.bed',

    ('dm6', 'flybase_r6.02', GeneIDType.SYMBOL, Delineation.UPSTREAM5KB_FULL_TX) : '/home/icistarget/data/regions/drosophila_melanogaster/dm6/flybase-r6.02/flybase-dmel-r6.02_symbol-limited-upstream5000-full-transcript_chr.bed',

    ('hg19', 'refseq_r45', GeneIDType.HGNC_SYMBOL, Delineation.UPSTREAM10KB_TSS_DOWNSTREAM10KB) : '/home/icistarget/data/regions/homo_sapiens/hg19/refseq-r45/hg19__refseq_r45__hg19-regions-from-uniform-20kbAround-closest-genes.closestGenesToCtxCoordinates.bed',

    ('mm9', 'refseq_r45', GeneIDType.MGI_SYMBOL, Delineation.UPSTREAM10KB_TSS_DOWNSTREAM10KB) : '/home/icistarget/data/regions/mus_musculus/mm9/refseq-r45/mm9__refseq_r45__mm9_regulatory_regions_DHS_20kbAroundTSS.closestGenesToCtxCoordinates.bed',

    ('mm9', 'refseq_r70', GeneIDType.MGI_SYMBOL, Delineation.UPSTREAM10KB_TSS_DOWNSTREAM10KB) : '/home/icistarget/data/regions/mus_musculus/mm9/refseq-r70/mm9__refseq_r70_regions_to_TSS_10kb_Up_Down.bed'
}


def bed2regions(assembly, gene_annotation_version, data, fraction=None):
    """Convert BED file entries to a set of regions."""
    regions = set()

    RegionFeatures.initialise(assembly, gene_annotation_version)

    for feature in RegionFeatures.REGION_FEATURES[(assembly, gene_annotation_version)].find_overlap_with(
            FeatureList.from_string(data),
            fraction):
        regions.add(feature.name)
    return regions


def input_bed2input_bed_and_overlap_icistarget_regions(assembly, gene_annotation_version, data, fraction=None):
    """Get overlap of input BED file with i-cisTarget regions."""
    input_bed_overlap_with_icistarget_regions = dict()

    RegionFeatures.initialise(assembly, gene_annotation_version)

    for input_bed_feature, icistarget_regions_feature in zip(
            RegionFeatures.REGION_FEATURES[(assembly, gene_annotation_version)].find_overlap_with_and_return_other(
                    FeatureList.from_string(data),
                    fraction),
            RegionFeatures.REGION_FEATURES[(assembly, gene_annotation_version)].find_overlap_with(
                    FeatureList.from_string(data),
                    fraction)
    ):
        input_bed_feature_key = "\t".join([input_bed_feature.chromosome, str(input_bed_feature.start),
                                           str(input_bed_feature.end), input_bed_feature.name])
        icistarget_regions_feature = "\t".join([icistarget_regions_feature.chromosome,
                                                str(icistarget_regions_feature.start),
                                                str(icistarget_regions_feature.end),
                                                icistarget_regions_feature.name])

        input_bed_overlap_with_icistarget_regions.setdefault(
            input_bed_feature_key,
            []).append(icistarget_regions_feature)

    for input_bed_feature in FeatureList.from_string(data):
        input_bed_feature_key = "\t".join([input_bed_feature.chromosome, str(input_bed_feature.start),
                                           str(input_bed_feature.end), input_bed_feature.name])
        if input_bed_feature_key in input_bed_overlap_with_icistarget_regions:
            for icistarget_regions_feature in input_bed_overlap_with_icistarget_regions[input_bed_feature_key]:
                yield input_bed_feature_key + '\t' + icistarget_regions_feature
        else:
            yield input_bed_feature_key + '\t\t\t\t'


def _remove_suffix(name):
    if '#' in name:
        idx = name.index('#')
        return name[:idx]
    else:
        return name


@cache(12)
def _get_feature_list(assembly, gene_annotation_version, gene_id_type, delineation):
    key = (assembly, gene_annotation_version, gene_id_type, delineation)
    return FeatureList.from_bed_file(ID_TYPE_DELINEATION_2_BED_FILENAME[key], _remove_suffix)


def get_fraction_of_mapped_gene_ids(assembly, gene_annotation_version, gene_ids, gene_id_type=GeneIDType.FBGN, delineation=Delineation.UPSTREAM5KB_FULL_TX):
    all_gene_ids = set(_get_feature_list(assembly, gene_annotation_version, gene_id_type, delineation).feature_ids)
    return float(len(all_gene_ids & gene_ids))/len(gene_ids)


def gene_ids2regions(assembly, gene_annotation_version, gene_ids, gene_id_type=GeneIDType.FBGN, delineation=Delineation.UPSTREAM5KB_FULL_TX, fraction=None):
    """ Convert gene identifiers to a set of regions. """
    region_ids = set()

    RegionFeatures.initialise(assembly, gene_annotation_version)

    for feature in RegionFeatures.REGION_FEATURES[(assembly, gene_annotation_version)].find_overlap_with(
            _get_feature_list(assembly, gene_annotation_version, gene_id_type, delineation).filter_by_name(gene_ids),
            fraction):
        region_ids.add(feature.name)
    return region_ids


# TODO: Algorithm in pseudo-code:
# TODO: Prerequisites: TSSs for all IDs and all different gene nomeclatures (Use FeatureList with interval tree and lookup by name as fast operations)
# TODO: 1. For a given gene ID in signature => find TSS (or multiple TSSs if more than one transcript is defined for this gene)
# TODO: 2. For regions that are associated with this ID (cfr. previous algorithm):
# TODO:      2.A. If this region contains the TSS of this gene: include this region as part of the foreground set
# TODO:      2.B. Else if minimum distance of regions to TSS of this gene is below 2kb: include this region as part of the foreground set
# TODO:      2.C. Else if this candidate region does not contain any of the defined TSSs for the gene nomenclature: include this region as part of the foreground set
# TODO:      2.D. Exclude region from foreground set
def ids2regions_new(assembly, gene_annotation_version, gene_ids, gene_id_type=GeneIDType.FBGN, delineation=Delineation.UPSTREAM5KB_FULL_TX, fraction=None):
    pass
