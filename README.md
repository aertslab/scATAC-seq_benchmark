# scATAC-seq_benchmark
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
# Reproducing manuscript figures
Here you can find where the code for each figure in the manuscript can be found:  

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

Supplementary figures:  
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
Supplementary table with metadata: `general/fixedcells_general_statistics.ipynb`
