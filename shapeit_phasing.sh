#!/bin/bash
for CHR in {1..22}; do
    UNPHASED_VCF=sample1.chr$CHR.vcf.gz
    UNPHASED_VCF_NOMAS=sample1.chr$CHR.noMultiAllelicSites.vcf.gz
    GEN_MAP=genetic_map_chr${CHR}_combined_b37.txt
    REF_HAPS=1000GP_Phase3_chr$CHR.hap.gz
    REF_LEGEND=1000GP_Phase3_chr$CHR.legend.gz
    REF_SAMPLE=1000GP_Phase3.sample
    LOG_ALIGN=sample1.chr$CHR.alignments
    EXCLUDE_LIST=$LOG_ALIGN.strand.exclude
    LOG_MAIN=sample1.chr$CHR.main
    
    PHASED_HAPS=sample1.chr$CHR.phased.haps.gz
    PHASED_SAMPLE=sample1.chr$CHR.phased.samples
    PHASED_VCF=sample1.chr$CHR.onlyPhased.vcf
    
    LOG_CONVERT=sample1.chr$CHR.convert
    FINAL_VCF=sample1.chr$CHR.phased.vcf.gz
    
    
	#Preparation
	bcftools view -M 2 -O z $UNPHASED_VCF > $UNPHASED_VCF_NOMAS
	shapeit -check -V $UNPHASED_VCF_NOMAS -M $GEN_MAP --input-ref $REF_HAPS $REF_LEGEND $REF_SAMPLE --output-log $LOG_ALIGN

	#Main run
	shapeit -V $UNPHASED_VCF_NOMAS -M $GEN_MAP --input-ref $REF_HAPS $REF_LEGEND $REF_SAMPLE -O $PHASED_HAPS $PHASED_SAMPLE --exclude-snp $EXCLUDE_LIST --no-mcmc --output-log $LOG_MAIN
	
    shapeit -convert --input-haps $PHASED_HAPS $PHASED_SAMPLE --output-vcf $PHASED_VCF --output-log $LOG_CONVERT

	#Zipping and indexing
	bcftools view -O z $PHASED_VCF > $PHASED_VCF.gz
	bcftools index -f $PHASED_VCF

	#Merging phased and unphased vcfs, keeping all unphased sites from the original vcf, but replacing the phased ones.
	bcftools merge --force-samples $UNPHASED_VCF $PHASED_VCF | awk 'BEGIN {ofs=”\t”}
        $0 ~ /^#CHROM/ {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}
        $0 !~ /^#/ {
            if(substr($11, 1, 3) != "./.")
                $10 = $11
            print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10
        }' | bcftools view -O z > $FINAL_VCF
done
