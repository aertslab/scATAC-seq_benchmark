#!/usr/bin/env python3

from itertools import combinations
import os

base_dir = 'out_fmx/data/freemuxlet'
samples = [
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
    's3atac',
    ]

f_fmx_vcfs = {
        x: os.path.join(base_dir, x + '_freemuxlet.clust1.vcf.gz') for x in samples
        }


# bcftools filters:
min_DP = 10
min_GQ = 10

fmx_samples = ['CLUST0', 'CLUST1']

intermediate_dir = 'out_fmx/intermediate'

# the file {sample_rename} contains a single line with the new sample name to use in place of the existing CLUST0, CLUST1, etc.

filtered_vcfs = []
for n,(fmxVCFname,fmxVCF) in enumerate(f_fmx_vcfs.items()):
    for sample in fmx_samples:
        f_vcf = f"{fmxVCFname}.{sample}.vcf.gz"
        filtered_vcfs.append(f_vcf)
        cmd = f"bcftools view -s {sample} -i 'MIN(FMT/DP)>{min_DP} & MIN(FMT/GQ)>{min_GQ}' {fmxVCF} -O v | bcftools reheader -s out_fmx/sample_rename | bcftools view -O z -o {intermediate_dir}/{f_vcf} && tabix -p vcf -f {intermediate_dir}/{f_vcf}"
        print(cmd)


for i,x in enumerate(filtered_vcfs):
    for j,y in enumerate(filtered_vcfs):
        v0 = x.replace(".vcf.gz","")
        v1 = y.replace(".vcf.gz","")
        cmd = f"vcf-compare -g {intermediate_dir}/{x} {intermediate_dir}/{y} > {intermediate_dir}/{v0}_{v1}.vcfcompare"
        print(cmd)

