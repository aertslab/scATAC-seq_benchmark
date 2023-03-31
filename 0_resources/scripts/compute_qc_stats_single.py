import collections as cl
import gc
import logging
import sys
from typing import Dict, List, Optional, Tuple, Union
import os
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyranges as pr
from scipy.stats import gaussian_kde, norm

from .cistopic_class import *
from .utils import (
    collapse_duplicates,
    multiplot_from_generator,
    read_fragments_from_file,
)



def make_argument_parser():
    """
    Creates an ArgumentParser to read the options for this script from
    sys.argv
    """
    parser = argparse.ArgumentParser(
        description="Create metadata",)
    parser.add_argument('--fragments', '-f', type=str, required=True,
                        help='Path to fragments file.')
    
    parser.add_argument('--label', '-l', type=str, required=True,
                        help='Label name')
    
    parser.add_argument('--regions_path', '-r', type=str, required=True,
                        help='Path to regions in bed format.')
    
    parser.add_argument('--valid_bc_path', '-b', type=str, required=True,
                        help='Path to valid bc pickle file.')
    
    parser.add_argument('--genome_path', '-g', type=str, required=True,
                        help='Path to genome_annotation.')
    
    parser.add_argument('--genome', '-g', type=str, required=True,
                        help='genome_version.')
    
    parser.add_argument('--n_frag', '-n', type=int, required=True,
                        help='TSS flank window.')
    
    parser.add_argument('--tss_flank_window', '-t', type=int, required=True,
                        help='TSS flank window.')
    
    parser.add_argument('--tss_window', '-tt', type=int, required=True,
                        help='tss_window')
    
    parser.add_argument('--tss_minimum_signal_window', '-ttt', type=int, required=True,
                        help='tss_minimum_signal_window')
    
    parser.add_argument('--tss_rolling_window', '-tttt', type=int, required=True,
                        help='tss_rolling_window')
    
    parser.add_argument('--min_norm', '-m', type=int, required=True,
                        help='Path to cisTopic object pickle file.')
    
    parser.add_argument('--partition', '-p', type=int, required=True,
                        help='Path to cisTopic object pickle file.')
    
    parser.add_argument('--check_for_duplicates', '-c', type=bool, required=True,
                        help='Path to cisTopic object pickle file.')
    
    parser.add_argument('--use_polars', '-c', type=bool, required=True,
                        help='Path to cisTopic object pickle file.')
    
    parser.add_argument('--metadata_out', '-c', type=str, required=True,
                        help='metadata out path.')
    
    parser.add_argument('--profile_data_out', '-c', type=str, required=True,
                        help='profile data out path.')
    
    return parser
    
def main():
    """
    The main executable function
    """

    parser = make_argument_parser()
    
    args = parser.parse_args()
    fragments = args.fragments
    label = args.fragments
    path_to_regions = args.regions_path
    valid_bc = args.regions_path
    tss_flank_window = args.tss_flank_window
    tss_window = args.tss_window
    tss_minimum_signal_window = args.tss_minimum_signal_window
    tss_rolling_window = args.tss_rolling_window
    min_norm = args.min_norm
    partition = args.partition
    min_frag = args.n_frag
    check_for_duplicates = args.check_for_duplicates
    use_polars = args.use_polars
    genome = args.genome_path
    metadata_out = args.metadata_out
    profile_data_out = args.profile_data_out

    stats = [
        "barcode_rank_plot",
        "duplicate_rate",
        "insert_size_distribution",
        "profile_tss",
        "frip",
    ]
    
    """
    Parameters
    ---
    fragments: str
            Path to fragments file.
    tss_annotation: pd.DataFrame or pr.PyRanges
            A data frame or pyRanges containing transcription start sites for each gene, with 'Chromosome', 'Start' and 'Strand' as columns (additional columns will be ignored).
    stats: list, optional
            A list with the statistics that have to be computed. Default: All ('barcode_rank_plot', 'duplicate_rate', 'insert_size_distribution', 'profile_tss', 'FRIP).
    label: str
            Sample label. Default: None.
    path_to_regions: str
            Path to regions file to use for FRIP.
    valid_bc: list, optional
            A list containing selected barcodes. This parameter is ignored if n_frag or n_bc are specified. Default: None.
    n_frag: int, optional
            Minimal number of fragments assigned to a barcode to be kept. Either n_frag or n_bc can be specified. Default: None.
    n_bc: int, optional
            Number of barcodes to select. Either n_frag or n_bc can be specified. Default: None.
    tss_window: int, optional
            Window around the TSS used to count fragments in the TSS when calculating the TSS enrichment per barcode. Default: 50 (+/- 50 bp).
    tss_flank_window: int, optional
            Flanking window around the TSS. Default: 1000 (+/- 1000 bp).
    tss_minimum_signal_window: int, optional
            Tail window use to normalize the TSS enrichment. Default: 100 (average signal in the 100bp in the extremes of the TSS window).
    tss_rolling_window: int, optional
            Rolling window used to smooth signal. Default: 10.
    min_norm: int, optional
            Minimum normalization score. If the average minimum signal value is below this value, this number is used to normalize the TSS signal. This approach penalizes cells with fewer reads.
    check_for_duplicates: bool, optional
            If no duplicate counts are provided per row in the fragments file, whether to collapse duplicates. Default: True.
    remove_duplicates: bool, optional
            Whether to remove duplicates. Default: True.
    use_polars: bool, optional
            Whether to use polars to read fragments files. Default: True.
    Return
    ---
    pd.DataFrame or list and list
            A list with the barcode statistics for all samples (or a combined data frame with a column 'Sample' indicating the sample of origin) and a list of dictionaries with the sample-level profiles for each sample.
    """
    pbm_genome_name_dict = {
        "hg38": "hsapiens_gene_ensembl",
        "hg37": "hsapiens_gene_ensembl",
        "mm10": "mmusculus_gene_ensembl",
        "dm6": "dmelanogaster_gene_ensembl",
    }

    pbm_host_dict = {
        "hg38": "http://www.ensembl.org",
        "hg37": "http://grch37.ensembl.org/",
        "mm10": "http://nov2020.archive.ensembl.org/",
        "dm6": "http://www.ensembl.org",
    }

    if os.path.exists(genome_path):
        print(f"Loading cached genome annotation...")
        annotation = pd.read_csv(genome_path, sep="\t", header=0, index_col=0)
    else:
        dataset = pbm.Dataset(name=pbm_genome_name_dict[genome], host=pbm_host_dict[genome])

        annotation = dataset.query(
            attributes=[
                "chromosome_name",
                "transcription_start_site",
                "strand",
                "external_gene_name",
                "transcript_biotype",
            ]
        )
        filter = annotation["Chromosome/scaffold name"].str.contains("CHR|GL|JH|MT")
        annotation = annotation[~filter]
        annotation["Chromosome/scaffold name"] = annotation[
            "Chromosome/scaffold name"
        ].str.replace(r"(\b\S)", r"chr\1")
        annotation.columns = ["Chromosome", "Start", "Strand", "Gene", "Transcript_type"]
        annotation = annotation[annotation.Transcript_type == "protein_coding"]
        annotation.to_csv(f"{genome}_annotation.tsv", sep="\t")

    tss_annotation = annotation

    # Create logger
    level = logging.INFO
    log_format = "%(asctime)s %(name)-12s %(levelname)-8s %(message)s"
    handlers = [logging.StreamHandler(stream=sys.stdout)]
    logging.basicConfig(level=level, format=log_format, handlers=handlers)
    log = logging.getLogger("cisTopic")
    
    # Compute stats
    metrics = {}
    metadata_bc_dict = {}
    profile_data_dict = {}
    # Prepare fragments
    if isinstance(fragments, str):
        log.info("Reading " + label)
        fragments_df = read_fragments_from_file(fragments, use_polars=use_polars).df
    else:
        fragments_df = fragments
    # Convert to category for memory efficiency
    fragments_df["Name"] = fragments_df["Name"].astype("category")
    # Check for duplicates
    if "Score" not in fragments_df or all(fragments_df["Score"] == "."):
        fragments_df = fragments_df[["Chromosome", "Start", "End", "Name"]]
        if check_for_duplicates:
            log.info("Collapsing duplicates")
            fragments_df = pd.concat(
                [
                    collapse_duplicates(fragments_df[fragments_df.Chromosome == x])
                    for x in fragments_df.Chromosome.cat.categories.values
                ]
            )
        else:
            fragments_df["Score"] = 1
    else:
        fragments_df = fragments_df[["Chromosome", "Start", "End", "Name", "Score"]]
    fragments_df["Score"] = fragments_df["Score"].astype("int32")
    # Prepare valid barcodes
    if valid_bc is not None:
        if n_bc is not None or n_frag is not None:
            valid_bc = None
    # Rank plot
    if "barcode_rank_plot" in stats:
        # Rank plot
        log.info("Computing barcode rank plot for " + label)
        metrics["barcode_rank_plot"] = barcode_rank_plot(
            fragments=fragments_df,
            valid_bc=valid_bc,
            n_frag=n_frag,
            n_bc=n_bc,
            remove_duplicates=remove_duplicates,
            plot=False,
            return_bc=True,
            return_plot_data=True,
        )
        if valid_bc is None:
            fragments_df = fragments_df[
                fragments_df.Name.isin(set(metrics["barcode_rank_plot"]["valid_bc"]))
            ]

    gc.collect()
    # Duplicate rate
    if "duplicate_rate" in stats:
        # Duplicate rate
        log.info("Computing duplicate rate plot for " + label)
        metrics["duplicate_rate"] = duplicate_rate(
            fragments=fragments_df, valid_bc=valid_bc, plot=False, return_plot_data=True
        )

    gc.collect()
    # Fragment size
    if "insert_size_distribution" in stats:
        # Fragment size
        log.info("Computing insert size distribution for " + label)
        metrics["insert_size_distribution"] = insert_size_distribution(
            fragments=fragments_df,
            valid_bc=valid_bc,
            remove_duplicates=remove_duplicates,
            plot=False,
            return_plot_data=True,
        )
    fragments_df = pr.PyRanges(fragments_df)
    gc.collect()
    # TSS
    if "profile_tss" in stats:
        # TSS
        log.info("Computing TSS profile for " + label)
        profile_tss_metrics = profile_tss(
            fragments=fragments_df,
            annotation=tss_annotation,
            valid_bc=valid_bc,
            plot=False,
            n_cpu=1,
            partition=partition,
            flank_window=tss_flank_window,
            tss_window=tss_window,
            minimum_signal_window=tss_minimum_signal_window,
            rolling_window=tss_rolling_window,
            min_norm=min_norm,
            return_TSS_enrichment_per_barcode=True,
            return_TSS_coverage_matrix_per_barcode=True,
            return_plot_data=True,
        )
        if profile_tss_metrics is not None:
            metrics["profile_tss"] = profile_tss_metrics
    gc.collect()
    # FRIP
    if "frip" in stats:
        # FRIP
        log.info("Computing FRIP profile for " + label)
        metrics["frip"] = frip(
            fragments=fragments_df,
            path_to_regions=path_to_regions,
            valid_bc=valid_bc,
            remove_duplicates=remove_duplicates,
            n_cpu=1,
            plot=False,
            return_plot_data=True,
        )
    del fragments_df
    gc.collect()
    metadata_bc, profile_data = metrics2data(metrics)

    if isinstance(metadata_bc, pd.DataFrame):
        metadata_bc = metadata_bc.fillna(0)
    
    with open(
        metadata_out, "wb"
    ) as f:
        pickle.dump(metadata_bc, f, protocol=4)

    with open(
        profile_data_out, "wb"
    ) as f:
        pickle.dump(profile_data, f, protocol=4)

if __name__ == "__main__":
    main()