#!/usr/bin/env python3


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


fmx_samples = ['CLUST0', 'CLUST1']
intermediate_dir = 'out_fmx/intermediate'

nvar = []
nrd = []
for r1ix,r1 in enumerate(samples):
    for s1ix,s1 in enumerate(fmx_samples):
        nvar_tmp2 = [f"{r1}.{s1}"]
        nrd_tmp2 = [f"{r1}.{s1}"]
        cnames = []
        for r2ix,r2 in enumerate(samples):
            for s2ix,s2 in enumerate(fmx_samples):
                f_comp = f"{intermediate_dir}/{r1}.{s1}_{r2}.{s2}.vcfcompare"
                with open(f_comp,'r') as f:
                    dat = f.read().splitlines()
                # number of variants compared ((fix this!)):
                # the VN rows contain overlap/unique counts per file
                # we need to check if there's an entry for both vcfs; this would indicate there are variants overlapping
                # absence of such an entry would indicate no overlapping variants
                tmp_nvar = [ y for y in dat if y.startswith("VN") ]
                tmp_nvar_split = [ x.split('\t') for x in tmp_nvar ]
                v_overlap = 0
                for i,x in enumerate(tmp_nvar_split):
                    vcf_cnt = sum([ "vcf.gz" in y for y in x ])
                    if(vcf_cnt==2):
                        # number of variants overlapping:
                        v_overlap = x[1]
                nvar_tmp2.append(str(v_overlap))

                ## # old method of getting number of variants overlapping: (doesn't work correctly)
                ## tmp_sn = [ y.split('\t')[2] for y in dat if y.startswith("SN") ][:3]
                ## tmp_nvar = sum([ int(y) for y in tmp_sn ])
                ## nvar_tmp2.append(str(tmp_nvar))

                # non reference discord rate:
                #if(v_overlap==0):
                #    tmp_nrd = "NA"
                #else:
                tmp_nrd = [ y for y in dat if y.startswith("SN\tNon-reference") ]
                tmp_nrd = tmp_nrd[0].split('\t')[-1]
                nrd_tmp2.append(tmp_nrd)
                cnames.append(f"{r2}.{s2}")
        nrd.append(nrd_tmp2)
        nvar.append(nvar_tmp2)

# write to file:
f_nvar = open('out_fmx/genotype_concordance_nvar.txt','w')
f_nrd = open('out_fmx/genotype_concordance_nrd.txt','w')
print('\t'.join(['']+cnames), file=f_nvar)
print('\t'.join(['']+cnames), file=f_nrd)
for i,x in enumerate(nvar):
    print('\t'.join(x), file=f_nvar)
    print('\t'.join(nrd[i]), file=f_nrd)
f_nvar.close()
f_nrd.close()

print("Created two files:\n\tgenotype_concordance_nvar.txt\n\tgenotype_concordance_nrd.txt\n")

