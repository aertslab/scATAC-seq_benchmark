#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
sample_id = args[1]
f_loom = args[2]
f_frag = args[3]
f_reference = args[4]
f_annotation = args[5]
f_out = args[6]

print(paste0("Processing sample ", args[1]))

# load pbmc object
# pbmc.integrated <- readRDS("../0_resources/seurat_references/pbmc_integrated.RDS")
pbmc.rna <- readRDS('../0_resources/seurat_references/pbmc_ssc_mat__integrated.rds')

################################################################################
# ATAC
################################################################################

### get data from loom:
atacloomcon <- Connect(filename = f_loom, mode = "r")
atacloomcon
atac_tmp <- as.Seurat(atacloomcon, assay='ATAC')
atacloomcon$close_all()

# correctly parse regions (default delims are '-','-')
regions = StringToGRanges(
    rownames(GetAssayData(atac_tmp, slot = "counts", assay='ATAC')),
    sep=c(':','-')
    )

# create chromatin assay
chromatinassay = CreateChromatinAssay(
    counts=GetAssayData(atac_tmp, slot = "counts", assay='ATAC'),
    # genome='hg38',
    fragments = f_frag,
    ranges=regions
    )
    #annotation=annotation)

atac <- CreateSeuratObject(counts = chromatinassay, assay='ATAC')

annotations <- readRDS(f_annotation)
Annotation(atac) <- annotations

# We exclude the first dimension as this is typically correlated with sequencing depth
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunSVD(atac)
atac <- RunUMAP(atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

##################################################
# Identify RNA-ATAC anchors

# quantify gene activity
gene.activities <- GeneActivity(atac, features = VariableFeatures(pbmc.rna))
# gene.activities <- gene.activities[rownames(gene.activities)!="",]

# add gene activities as a new assay
atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(atac) <- "ACTIVITY"
atac <- NormalizeData(atac)
atac <- ScaleData(atac, features = rownames(atac))

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = pbmc.rna,
    query = atac,
    features = VariableFeatures(object = pbmc.rna),
    reference.assay = "RNA",
    query.assay = "ACTIVITY",
    reduction = "cca")

# predict celltype
celltype.predictions <- TransferData(
    anchorset = transfer.anchors,
    refdata = pbmc.rna$CellType,
    weight.reduction = atac[["lsi"]],
    dims = 2:30)

pbmc.atac <- AddMetaData(atac, metadata = celltype.predictions)

md = celltype.predictions
pred_thr = 0.7

tmp = data.frame(
      composite_sample_id = paste0(rownames(md),'-',sample_id),
      barcode = rownames(md),
      sample_id = sample_id,
      cell_type = md$predicted.id,
      cell_type_pred_score = md$prediction.score.max
      )

tmp$cell_type_hiconf_70 = tmp$cell_type
tmp$cell_type_hiconf_70[tmp$cell_type_pred_score<pred_thr] = 'Unknown'

table(tmp$cell_type)
table(tmp$cell_type_hiconf_70)

write.table(tmp, file=f_out, sep='\t', row.names=FALSE, quote=FALSE)
print(paste0('Wrote file ', f_out))
