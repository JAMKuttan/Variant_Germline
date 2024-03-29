#!/bin/bash
#snsp.sh

module load samtools/1.6 bedtools/2.26.0
ID=$1

baseDir="`dirname \"$0\"`"
perl ${baseDir}/filter_giab.pl $ID

bedtools intersect -header -a ${ID}.PASS.vcf -b /project/shared/bicf_workflow_ref/GRCh38/utswv2_cds.bed | bedtools intersect -v -header -a stdin -b /project/shared/bicf_workflow_ref/GRCh38/HLA_HG38.bed | bedtools intersect -header -a stdin -b /project/shared/bicf_workflow_ref/GRCh38/giab3_platinum2.hg38.highConf.bed | bgzip > ${ID}.t.vcf.gz
bcftools reheader -h /project/shared/bicf_workflow_ref/GRCh38/UTSWV2.vcfheader.txt -o ${ID}.utswcoding.vcf.gz ${ID}.t.vcf.gz
tabix ${ID}.utswcoding.vcf.gz
bedtools multiinter -i ${ID}.utswcoding.vcf.gz /project/shared/bicf_workflow_ref/GRCh38/giab3.utswcoding.vcf.gz  /project/shared/bicf_workflow_ref/GRCh38/platinum_v2.utswcoding.vcf.gz -names union giab platinum |cut -f 1,2,3,5 | bedtools sort -i stdin | bedtools merge -c 4 -o distinct | bgzip > ${ID}.utswcoding.multiinter.bed.gz 
tabix ${ID}.utswcoding.multiinter.bed.gz
bcftools annotate -a  ${ID}.utswcoding.multiinter.bed.gz --columns CHROM,FROM,TO,PlatRef -h /project/shared/bicf_workflow_ref/PlatRef.header ${ID}.utswcoding.vcf.gz > ${ID}.eval.vcf
perl /project/PHG/PHG_Clinical/clinseq_workflows/scripts/calc_snsp.pl ${ID}.eval.vcf
