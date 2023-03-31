# scATAC-seq_benchmark
These are all the Jupyter notebooks and scripts that were used to analyse data and generate figures for our paper "Systematic benchmarking of scATAC-seq protocols" (De Rop et al., 2023). With these scripts, you should be able to reproduce everything found in our manuscript.

## Directory structure
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
