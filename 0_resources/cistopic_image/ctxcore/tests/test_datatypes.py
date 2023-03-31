import pytest

from ctxcore.datatypes import (
    MotifOrTrackIDs,
    MotifsOrTracksType,
    RegionOrGeneIDs,
    RegionsOrGenesType,
    ScoresOrRankingsType,
)


def test_RegionsOrGenesType_from_str():
    """Check if a member of RegionsOrGenesType Enum can be made from a string."""
    assert RegionsOrGenesType.from_str("regions") == RegionsOrGenesType.REGIONS
    assert RegionsOrGenesType.from_str("REGIONS") == RegionsOrGenesType.REGIONS
    assert RegionsOrGenesType.from_str("genes") == RegionsOrGenesType.GENES
    assert RegionsOrGenesType.from_str("GENES") == RegionsOrGenesType.GENES

    with pytest.raises(
        ValueError,
        match=r'Unsupported RegionsOrGenesType "NON_EXISTING_REGIONS_OR_TRACKS_TYPE".',
    ):
        RegionsOrGenesType.from_str("NON_EXISTING_REGIONS_OR_TRACKS_TYPE")


def test_MotifsOrTracksType_from_str():
    """Check if a member of MotifsOrTracksType Enum can be made from a string."""
    assert MotifsOrTracksType.from_str("motifs") == MotifsOrTracksType.MOTIFS
    assert MotifsOrTracksType.from_str("MOTIFS") == MotifsOrTracksType.MOTIFS
    assert MotifsOrTracksType.from_str("tracks") == MotifsOrTracksType.TRACKS
    assert MotifsOrTracksType.from_str("TRACKS") == MotifsOrTracksType.TRACKS

    with pytest.raises(
        ValueError,
        match=r'Unsupported MotifsOrTracksType "NON_EXISTING_MOTIFS_OR_TRACKS_TYPE".',
    ):
        MotifsOrTracksType.from_str("NON_EXISTING_MOTIFS_OR_TRACKS_TYPE")


def test_ScoresOrRankingsType_from_str():
    """Check if a member of ScoresOrRankingsType Enum can be made from a string."""
    assert ScoresOrRankingsType.from_str("scores") == ScoresOrRankingsType.SCORES
    assert ScoresOrRankingsType.from_str("SCORES") == ScoresOrRankingsType.SCORES
    assert ScoresOrRankingsType.from_str("rankings") == ScoresOrRankingsType.RANKINGS
    assert ScoresOrRankingsType.from_str("RANKINGS") == ScoresOrRankingsType.RANKINGS

    with pytest.raises(
        ValueError,
        match=r'Unsupported ScoresOrRankingsType "NON_EXISTING_SCORES_OR_TRACK_TYPE".',
    ):
        ScoresOrRankingsType.from_str("NON_EXISTING_SCORES_OR_TRACK_TYPE")


def test_RegionOrGeneIDs_with_regions():
    """Check if a RegionOrGeneIDs object can be constructed from a list of region IDs."""
    region_or_gene_ids_instance = RegionOrGeneIDs(
        region_or_gene_ids=["reg2", "reg1", "reg6", "reg2"],
        regions_or_genes_type=RegionsOrGenesType.REGIONS,
    )
    assert region_or_gene_ids_instance.type == RegionsOrGenesType.REGIONS
    # Input region_or_gene_ids was a list ==> keep order.
    assert region_or_gene_ids_instance.ids == ("reg2", "reg1", "reg6")

    assert eval(repr(region_or_gene_ids_instance)) == region_or_gene_ids_instance
    assert len(region_or_gene_ids_instance) == 3


def test_RegionOrGeneIDs_with_genes():
    """Check if a RegionOrGeneIDs object can be constructed from a set of gene IDs."""
    region_or_gene_ids_instance = RegionOrGeneIDs(
        region_or_gene_ids={"gene2", "gene1", "gene6", "gene2"},
        regions_or_genes_type=RegionsOrGenesType.GENES,
    )
    assert region_or_gene_ids_instance.type == RegionsOrGenesType.GENES
    # Input region_or_gene_ids was a set ==> will be sorted.
    assert region_or_gene_ids_instance.ids == ("gene1", "gene2", "gene6")

    assert eval(repr(region_or_gene_ids_instance)) == region_or_gene_ids_instance
    assert len(region_or_gene_ids_instance) == 3


def test_RegionOrGeneIDs_with_regions_or_genes_type_str():
    """
    Check if a RegionOrGeneIDs object can be constructed from a tuple of gene IDs where regions_or_genes_type is given
    as a string.
    """
    region_or_gene_ids_instance = RegionOrGeneIDs(
        region_or_gene_ids=("gene2", "gene1", "gene6", "gene2"),
        regions_or_genes_type="gEnES",
    )
    assert region_or_gene_ids_instance.type == RegionsOrGenesType.GENES
    # Input region_or_gene_ids was a tuple ==> keep order.
    assert region_or_gene_ids_instance.ids == ("gene2", "gene1", "gene6")

    assert eval(repr(region_or_gene_ids_instance)) == region_or_gene_ids_instance
    assert len(region_or_gene_ids_instance) == 3


def test_RegionOrGeneIDs_subset_superset():
    """
    Check if region or gene IDs of a RegionOrGeneIDs object are a subset or a superset of another RegionOrGeneIDs
    object.
    """
    region_or_gene_ids_instance1 = RegionOrGeneIDs(
        region_or_gene_ids=["reg1", "reg2", "reg6"],
        regions_or_genes_type=RegionsOrGenesType.REGIONS,
    )
    region_or_gene_ids_instance2 = RegionOrGeneIDs(
        region_or_gene_ids=["reg1", "reg2", "reg4", "reg6"],
        regions_or_genes_type=RegionsOrGenesType.REGIONS,
    )

    assert region_or_gene_ids_instance1.issubset(region_or_gene_ids_instance2)
    assert region_or_gene_ids_instance2.issuperset(region_or_gene_ids_instance1)


def test_MotifsOrTracksIDs_with_motifs():
    """Check if a MotifOrTrackIDs object can be constructed from a list of motif IDs."""
    motif_or_track_ids_instance = MotifOrTrackIDs(
        motif_or_track_ids=["motif5", "motif10", "motif3", "motif10"],
        motifs_or_tracks_type=MotifsOrTracksType.MOTIFS,
    )
    assert motif_or_track_ids_instance.type == MotifsOrTracksType.MOTIFS
    # Input motif_or_track_ids was a list ==> keep order.
    assert motif_or_track_ids_instance.ids == ("motif5", "motif10", "motif3")

    assert eval(repr(motif_or_track_ids_instance)) == motif_or_track_ids_instance
    assert len(motif_or_track_ids_instance) == 3


def test_MotifsOrTracksIDs_with_tracks():
    """Check if a MotifOrTrackIDs object can be constructed from a set of track IDs."""
    motif_or_track_ids_instance = MotifOrTrackIDs(
        motif_or_track_ids={"track5", "track10", "track3", "track10"},
        motifs_or_tracks_type=MotifsOrTracksType.TRACKS,
    )
    assert motif_or_track_ids_instance.type == MotifsOrTracksType.TRACKS
    # Input motif_or_track_ids was a set ==> will be sorted.
    assert motif_or_track_ids_instance.ids == ("track10", "track3", "track5")

    assert eval(repr(motif_or_track_ids_instance)) == motif_or_track_ids_instance
    assert len(motif_or_track_ids_instance) == 3


def test_MotifsOrTracksIDs_with_motifs_or_tracks_type_str():
    """
    Check if a MotifOrTrackIDs object can be constructed from a tuple of track IDs,
    where motifs_or_tracks_type is given as a string.
    """
    motif_or_track_ids_instance = MotifOrTrackIDs(
        motif_or_track_ids=("track5", "track10", "track3", "track10"),
        motifs_or_tracks_type="tracks",
    )
    assert motif_or_track_ids_instance.type == MotifsOrTracksType.TRACKS
    # Input motif_or_track_ids was a tuple ==> keep order.
    assert motif_or_track_ids_instance.ids == ("track5", "track10", "track3")

    assert eval(repr(motif_or_track_ids_instance)) == motif_or_track_ids_instance
    assert len(motif_or_track_ids_instance) == 3
