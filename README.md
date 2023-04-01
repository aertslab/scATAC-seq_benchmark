# scATAC-seq Benchmark
These are all the Jupyter notebooks and scripts that were used to analyse data and generate figures for our paper "Systematic benchmarking of scATAC-seq protocols" (De Rop et al., 2023). With these scripts, you should be able to reproduce everything found in our manuscript.

## Directory structure
Here you can find the structure of the root directory, with descriptions of each subdirectory.
```
scATAC-seq_benchmark
├── 0_resources # all resources used (generic scripts, specific region sets, ...). Some files, such as SCREEN peak sets, or reference genomes, are too large for github, please contact us if you want to request these.
├── 1_data_repository # all the raw data (FASTQ) as well as fragments files. If you want to reproduce our analyses, you should source our raw data from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194028) 
├── fixedcells_1_vsn_preprocessing # processing steps of libds_1_vsn_preprocessing bams
├── fixedcells_2_cistopic # cisTopic pre-processing in SCREEN regions and calling consensus peaks
├── fixedcells_3_cistopic_consensus # cisTopic pre-processing in consensus peaks and calling master consensus peaks
├── fixedcells_4_merged # cisTopic pre-processing of merged dataset (169k cells), using master consensus peaks
├── fixedcells_5_cell_downsampling # random cell downsampling analyses
├── fixedcells_6_merged_equalcells # merged object containing equal number of cells from each technique, but not an equal number of cells for each cell type within each technique (not used in manuscript)
├── fixedcells_7_merged_equalcells_celltypefair # merged object containing equal number of cells from each technique and an equal number of cells for each cell type within each technique
├── fixedcells_8_individual_tech_cistopic_objects # individual objects for each technology containing equal number of cells from each technique and an equal number of cells for each cell type 
├── fixedcells_9_individual_malefemale_celltypefair # individual objects for each technology containing equal number of cells from each technique and an equal number of cells for each cell type and each donor
├── fixedcells_cellranger_arc # aligning downsampled data to reference genome using cellranger_arc, and performing some analyses on multiome-rna component
├── fixedcells_downsample_series # downsampling further to 35k, 30k, ... 5k reads per cell and analysing the resulting cisTopic objects
├── full_1_vsn_preprocessing # aligning full FASTQs to reference genome
├── full_2_cistopic # cisTopic pre-processing in SCREEN regions and calling consensus peaks
├── full_3_cistopic_consensus # cisTopic pre-processing in consensus peaks and calling master consensus peaks
├── full_4_merged # cisTopic pre-processing of merged dataset (169k cells), using master consensus peaks
├── full_5_cellranger # aligning full data to reference genome using cellranger, and validating ATACflow by comparing Cell Ranger and ATACflow outputs
├── general # general, sample-wide plots and definitions
├── libds_1_vsn_preprocessing # aligning downsampled FASTQs to the reference genome
├── public_1_cistopic_qc # analysing fragments files from public repositories, not used in manuscript
├── public_2_vsn_preprocessing # aligning full public mouse brain data to mouse genome
├── public_3_cistopic_qc # cisTopic pre-processing in SCREEN regions and calling consensus peaks
├── public_4_cistopic_consensus # cisTopic pre-processing in consensus peaks and calling master consensus peaks
└── public_downsample_series # downsampling analysis on public mouse cortex data
```

# How to interpret the directory structure
1. As new experiments were performed, sequencing data was deposited in `1_data_repository/original_fastq`. Each sample's sequencing data was then merged to a maximum of 3 files (barcode read, and two mates) and deposited in `1_data_repository/full_fastq`.
2. For each experiment, the full sequencing data was then aligned to the reference genome. Results are in `full_1_vsn_preprocessing`. Symlinks to `.bam` and `.fragments.tsv.gz` were placed in `1_data_repository/full_bams` and `1_data_repository/full_fragments`. We refer to VSN, as our pipeline at the time was still a part of [VSN](https://github.com/vib-singlecell-nf/vsn-pipelines), but now has its own repository [ATACflow](https://github.com/aertslab/ATACflow).
3. For each sample, we then filtered true cell barcodes from noise barcodes in `full_2_cistopic`. This filtering was performed using thresholds on TSS enrichment and number of unique fragments.
4. Since we then knew the number of cells present in each sample, we could downsample the full sequencing data to the same common read depth (40k reads/cell). This was performed using the notebook `1_data_repository/5_downsample_fastq.ipynb` and the downsampled FASTQs were deposited in `1_data_repository/libds_fastq`. `libds` stands for "library downsampled". For a long time, we then re-called cells in these FASTQ files and proceeded with analysis like this. For most samples, the number of cells called was very similar, but for some samples that were added later, there were large discrepancies. We were thus faced with a dilemma: either adapt our cell filtering algorithm so that for the new samples cell counts would be the same between full data and downsampled data, or simply take the list of filtered barcodes from the full data and re-do all the analysis on the downsampled data using this barcode list instead. We chose the latter approach. This new sampling strategy was then referred to as `fixedcells`, as the number and identity of cells was now fixed after identification in the full sequencing data.
5. We then re-aligned the downsampled FASTQs to GRCh38 in `libds_1_vsn_preprocessing` and Jaccard-process those files in `fixedcells_1_vsn_preprocessing`, as we did not not need to realign the downsampled sequencing data after switching to the `fixedcells` cell filtering strategy.
6. We performed cisTopic clustering, Freemuxlet donor assignment, Seurat cell type annotation and consensus peak calling in `fixedcells_2_cistopic`.
7. We re-count each sample's fragments in its own consensus peak set, re-do Seurat cell type annotation, and do all further downstream analysis (such as DAR calling and motif enrichment analysis) based on these count matrices. FRIP scores are also calculated using each sample's specific consensus peaks. Freemuxlet donor assignment was re-taken from the first pass done in `fixedcells_2_cistopic` because it is a bam-level analysis and independent of consensus peaks. We also re-calculated new consensus peak sets for each sample, and aggregated each of these second-pass consensus peak sets into one master peak set.
8. We recounted all fragments of all cells in all samples in the master peak set to generate a `fixedcells_merged` cisTopic object, and performed some analyses on the merged object in `fixedcells_4_merged`.
9. In `fixedcells_5_cell_downsampling`, we performed some analyses to investigate the effect of number of cells on some metrics, mostly Seurat cell type assignment and DAR calling. In order to do this, we subsampled each of the 47 individual `fixedcells` cisTopic objects to 2500, 2000, 1500, ... cells.
10. We attempted to do some analyses on the merged cisTopic object where each technology had the same number of cells (`fixedcells_6_merged_equalcells`), equal to the number of cells of the technology that had the least number of cells (s3-ATAC). However, at the same time, we were doing the downsampling analysis and saw that the number of cells *per cell type* also had an impact on the analysis... 
11. So, in `fixedcells_7_merged_equalcells_celltypefair`, we did the same, but now we took the same number of cells for each cell type for each technology, and in `fixedcells_8_individual_tech_cistopic_objects`, we simply split this merged object into 8 objects, one for each tech.
12. The same strategy was employed to create cisTopic objects for each technology that had the same number of cells for each cell type, but also from each donor within each cell type. Since s3-ATAC had so few cells compared to the rest, some concessions had to be made (the subsampling is done in `fixedcells_4_merged/9a_subset_malefemale.ipynb`, you can see the strategy used there).
13. In `full_5_cellranger`, we realigned all the 10x scATAC-seq data using Cell Ranger. We also performed our comparison with VSN there, and calibrated the Seurat scores using the multiome. We then filtered cells, and re-aligned the downsampled multiome data in `fixedcells_cellranger_arc`, and analysed the results (Venn diagram and correlations between ATAC and RNA counts).
14. In `fixedcells_downsample_series`, most of these analyses were performed on further read-downsampled data (35k, 30k, ... 5k reads/cell).
15. In `public_*` directories, all the public data was analysed, including a read downsampled analysis.

# Reproducing manuscript figures
You can find our ATACflow pipeline in [its own repository](https://github.com/aertslab/ATACflow). This pipeline can be used to realign data from all techniques assessed here to the reference genome.  

Here you can find where the code for each figure in the manuscript can be found:  

**Main figures:**  
1b: `general/fixedcells_merged_graphs.ipynb`  
1d: `general/scatterplots_bytech_kde_v2.py`  
1e-h: `general/fixedcells_boxplots.ipynb`  
1i-m: `fixedcells_downsample_series/4c_qc_plots.ipynb`  
2a-e: `general/fixedcells_general_statistics_scatterplots.ipynb`  
2f-i: `general/fixedcells_boxplots.ipynb`  
2j: `fixedcells_5_cell_downsampling/5b_DAR_scores.ipynb`  
2k-l: `fixedcells_8_individual_tech_cistopic_objects/4b_dar_scores.ipynb`  
2m: `fixedcells_7_merged_equalcells_celltypefair/4c_dar_traces.ipynb`  
2n-o: `fixedcells_8_individual_tech_cistopic_objects/7_peak_dar_overlap_samples.ipynb`  
3a: `fixedcells_8_individual_tech_cistopic_objects/8b_cistarget_analysis.ipynb`  
3b-c: `general/fixedcells_general_statistics_scatterplots.ipynb`  
3d-e: `fixedcells_9_individual_malefemale_celltypefair/5a_analyse_malefemale.ipynb`  
3f: `fixedcells_9_individual_malefemale_celltypefair/7b_male_female_tracks.ipynb`  
3g: `fixedcells_4_merged/5b_LISI.ipynb`  
3i: `fixedcells_4_merged/5b_LISI.ipynb`  

**Supplementary figures:**  
S1a: `full_5_cellranger/2b_validation_graphs.ipynb`  
S1b: `full_1_vsn_preprocessing/3_otsu_filtering.ipynb`  
S2-S3: `full_3_cistopic_consensus/9_plot_all_qc.ipynb`  
S4a-b: `fixedcells_cellranger_arc/2_cell_filtering.ipynb`  
S4c: `fixedcells_cellranger_arc/3_venn.ipynb`  
S5a-i: `general/fixedcells_general_statistics_gameshowell.ipynb`  
S6: `fixedcells_3_cistopic_consensus/3b_cell_type_analysis.ipynb`  
S7a: `fixedcells_downsample_series/5b_seurat_celltypes.ipynb`  
S7b-c: `fixedcells_downsample_series/7b_DARs_analysis.ipynb`  
S8a: `fixedcells_5_cell_downsampling/3_seurat_celltypes.ipynb`  
S8b-c: `fixedcells_5_cell_downsampling/5b_DAR_scores.ipynb`  
S9: `fixedcells_8_individual_tech_cistopic_objects/7_peak_dar_overlap_samples.ipynb`  
S10: `fixedcells_7_merged_equalcells_celltypefair/4d_dar_carrot.ipynb`  
S11a: `fixedcells_2_cistopic/2b_analyse_freemuxlet.ipynb`  
S11b: `fixedcells_3_cistopic_consensus/0_deteremine_male_vs_female.ipynb`  
S11c: `fixedcells_2_cistopic/2b_analyse_freemuxlet.ipynb`  
S12a: `public_4_cistopic_consensus/2_plot_all_qc.ipynb`  
S12b: `public_downsample_series/5_analyse_qc.ipynb`  
S13a-b: `full_4_merged/8_lisi.ipynb`  
S14: `1_data_repository/9_saturation_analysis.ipynb`  
S15: `fixedcells_3_cistopic_consensus/1b_count_fragments_in_blacklist.ipynb`  
S16: `full_5_cellranger/5_compare_rna_atac_seurat.ipynb`

Supplementary files:  
Supplementary table with quality control statistics: `general/fixedcells_general_statistics.ipynb`

# Contributing authors
All of these analyses were performed at the Stein Aerts lab by Florian De Rop, but they were largely based on a strong foundation laid by Christopher Flerin, who designed the initial analysis workflow. Gert Hulselmans also played a major role, as he designed [ATACflow](https://github.com/aertslab/ATACflow) (then still part of [VSN](https://github.com/vib-singlecell-nf/vsn-pipelines)) together with Christopher, and wrote most of the low-level scripts that work at the fragments and FASTQ level (calling `bwa-mem`, detecting and correcting barcodes, writing fragments files, calculating Jaccard indices, calling and speeding up Freemuxlet, subsampling BAM files, ...). This benchmark was supervised by Holger Heyn and Stein Aerts, who coordinated all work shown here and helped form major decisions at critical points.

All work shown here was done with the highest regard for fairness and transparency. If you have any questions, suggestions or criticisms, please contact us or open a github issue.
