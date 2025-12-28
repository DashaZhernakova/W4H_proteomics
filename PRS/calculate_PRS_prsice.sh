
##
### target data QC
##

# cd /mnt/sannaLAB-Temp/dasha/PRS/data/microarray
# 
# plink2 \
#     --bfile W4H_forMergingWith16S_geno01_hwe5e-6 \
#     --chr 1-22,X,Y \
#     --maf 0.01 \
#     --hwe 1e-6 \
#     --geno 0.01 \
#     --mind 0.01 \
#     --make-bed \
#     --out microarray.QC 
#     #--set-all-var-ids @:# 
# 
# sed -i "s:GSA-::g" microarray.QC.bim
# 
# sh ../../scripts/liftover.sh microarray.QC hg19 ../ref_data/hg38ToHg19.over.chain
# 
# plink2 \
# --bfile microarray.QC.hg19 \
# --rm-dup force-first list \
# --make-bed \
# --out microarray.QC.hg19.rmdup

##
### Base data QC
##
cd /mnt/sannaLAB-Temp/dasha/PRS/data/

protein=$1
prot_id=`awk -F'\t' -v prot="$protein" '$3 == prot {print $1}' /mnt/sannaLAB-Temp/dasha/PRS/data/olink_protein_names.txt`

echo "Processing $protein, $prot_id"

base_data=pQTLs/invn_${prot_id}_UKBB_PPP_F_only_regenie.gz 
target_data=microarray/microarray.QC.hg19.rmdup
ld_ref=ref_data/1kg_mrcieu/EUR
out_prefix=../results_p1e-3/PRSice/${protein}


cut -f2 ${target_data}.bim > ${target_data}.snplist
zcat $base_data | awk -F'\t' '
    NR == FNR { snps[$1] = 1; next }  # Load SNP IDs from .bim into memory
    FNR == 1 { print; next }          # Print header
    $2 in snps' ${target_data}.snplist - > ${base_data%gz}in_microarray.txt


### PRSice ###

# try1
/mnt/sannaLAB-Temp/dasha/tools/PRSice/PRSice_linux \
--base ${base_data%gz}in_microarray.txt \
--snp ID --chr CHROM --bp GENPOS --A1 ALLELE1 --A2 ALLELE0 --stat BETA --beta --pvalue Pval \
--target $target_data \
--ld $ld_ref \
--score avg \
--pheno olink_clean_CVD+INF_rm_outliers_4sd_average_across_timepoints.txt \
--pheno-col $protein \
--ignore-fid \
--out ${out_prefix}

#try2
/mnt/sannaLAB-Temp/dasha/tools/PRSice/PRSice_linux \
--base ${base_data%gz}in_microarray.txt \
--snp ID --chr CHROM --bp GENPOS --A1 ALLELE1 --A2 ALLELE0 --stat BETA --beta --pvalue Pval \
--target $target_data \
--ld $ld_ref \
--clump-p 1e-3 \
--score avg \
--pheno olink_clean_CVD+INF_rm_outliers_4sd_average_across_timepoints.txt \
--pheno-col $protein \
--ignore-fid \
--out ${out_prefix} \
--extract ${out_prefix}.valid


rm ${base_data%gz}in_microarray.txt




