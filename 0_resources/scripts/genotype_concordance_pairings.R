f_fmx_dir = 'out_fmx/data/freemuxlet/'

samples = c(
    'Broad_1',
    'Broad_2',
    'Broad_mito_1',
    'Broad_mito_2',
    'CNAG_1',
    'CNAG_2',
    'Sanger_1',
    'Sanger_2',
    'Stanford_1',
    'Stanford_2',
    'VIB_1',
    'VIB_2',
    'VIB_Hydrop_1',
    'VIB_Hydrop_2',
    's3atac'
    )

f_fmx = sapply(samples, function(s) paste0(f_fmx_dir,s,'_freemuxlet.clust1.samples.gz'))

fmx = lapply(f_fmx, read.delim, sep='\t', stringsAsFactors=FALSE)

# assign unique barcode:
for(i in 1:length(fmx)) {
    fmx[[i]]$ubarcode = paste(names(fmx)[i], "#", fmx[[i]]$BARCODE, sep="")
    fmx[[i]]$replicate = names(fmx)[i]
}


################################################################################
# freemuxlet sample pairing
################################################################################

# these groupings are made manually:
# based on the CLUST? in the heatmap

sampleA = list(
    c('Broad_1',1),
    c('Broad_2',0),
    c('Broad_mito_1',1),
    c('Broad_mito_2',0),
    c('CNAG_1',1),
    c('CNAG_2',1),
    c('Sanger_1',1),
    c('Sanger_2',1),
    c('Stanford_1',0),
    c('Stanford_2',1),
    c('VIB_1',1),
    c('VIB_2',1),
    c('VIB_Hydrop_1',1),
    c('VIB_Hydrop_2',1),
    c('s3atac',1)
)

sampleB = lapply(sampleA, function(A) c(A[1], ifelse(A[2]==0,1,0)) )


fmxA = list()
for(i in 1:length(sampleA)) {
    rep = sampleA[[i]][1]
    clust = sampleA[[i]][2]
    tmp = fmx[[rep]]
    bestguess = paste(clust,',',clust,sep="")
    ix = which(tmp$BEST.GUESS==bestguess)
    tmp = tmp[ix,]
    tmp$replicate = rep
    tmp$sample = "sampleA"
    fmxA[[i]] = tmp
}
fmxB = list()
for(i in 1:length(sampleB)) {
    rep = sampleB[[i]][1]
    clust = sampleB[[i]][2]
    tmp = fmx[[rep]]
    bestguess = paste(clust,',',clust,sep="")
    ix = which(tmp$BEST.GUESS==bestguess)
    tmp = tmp[ix,]
    tmp$replicate = rep
    tmp$sample = "sampleB"
    fmxB[[i]] = tmp
}

fmxA = do.call("rbind",fmxA)
fmxB = do.call("rbind",fmxB)


# unified table:
fmx = do.call(rbind,fmx)

# check there are no unexpected intersections
intersect(fmxA$ubarcode,fmxB$ubarcode)

ixA = which(fmx$ubarcode %in% fmxA$ubarcode)
ixB = which(fmx$ubarcode %in% fmxB$ubarcode)
intersect(ixA,ixB)


fmx$sample = 'Undetermined'
fmx$sample[ixA] = "sampleA"
fmx$sample[ixB] = "sampleB"

write.table(fmx, file="out_fmx/genotype_concordance_unified.txt", sep='\t', quote=FALSE )





################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

