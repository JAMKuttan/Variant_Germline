#!/usr/bin/env nextflow

params.input = "$baseDir"
params.fastqs="$params.input/*.fastq.gz"
params.design="$params.input/design.txt"
params.output = "$baseDir/output"
params.aligner='bwa'

params.genome="/project/shared/bicf_workflow_ref/human/GRCh38"
params.capture="$params.genome/clinseq_prj/hemepanelV3.bed"
params.pairs="pe"
params.cancer='skip'
params.callers = 'all'
params.markdups='picard'
params.callsvs="skip"


reffa=file("$params.genome/genome.fa")
dbsnp="$params.genome/dbSnp.vcf.gz"
indel="$params.genome/GoldIndels.vcf.gz"

fastqs=file(params.fastqs)
design_file = file(params.design)
dbsnp=file(dbsnp)
knownindel=file(indel)
index_path = file(params.genome)
capture_bed = file(params.capture)

alignopts = ''
if (params.markdups == 'fgbio_umi') {
   alignopts='-u'
}

def fileMap = [:]

fastqs.each {
    final fileName = it.getFileName().toString()
    prefix = fileName.lastIndexOf('/')
    fileMap[fileName] = it
}
def prefix = []
def read1_files = []
def read2_files = []
def mername = []
new File(params.design).withReader { reader ->
    def hline = reader.readLine()
    def header = hline.split("\t")
    prefixidx = header.findIndexOf{it == 'SampleID'};
    fidx = header.findIndexOf{it == 'FamilyID'};
    sidx = header.findIndexOf{it == 'SubjectID'};
    oneidx = header.findIndexOf{it == 'FqR1'};
    twoidx = header.findIndexOf{it == 'FqR2'};
    if (twoidx == -1) {
       twoidx = oneidx
       }
    if (sidx == -1) {
       sidx = prefixidx
    }
    if (fidx == -1) {
       fidx = sidx
    }
    if (twoidx == -1) {
       twoidx = oneidx
    }
    while (line = reader.readLine()) {
    	   def row = line.split("\t")
    if (fileMap.get(row[oneidx]) != null) {
	prefix << tuple(row[fidx],row[prefixidx],fileMap.get(row[oneidx]),fileMap.get(row[twoidx]))
	   }

}
}

if( ! prefix) { error "Didn't match any input files with entries in the design file" }

Channel
  .from(prefix)
  .set { read }

process trim {
  errorStrategy 'ignore'
  publishDir "$params.output/$famid/$pair_id", mode: 'copy'

  input:
  set famid,pair_id, file(read1), file(read2) from read

  output:
  set famid,pair_id, file("${pair_id}.trim.R1.fastq.gz"),file("${pair_id}.trim.R2.fastq.gz") into trimread
   set famid,pair_id, file("${pair_id}.trim.R1.fastq.gz"),file("${pair_id}.trim.R2.fastq.gz") into genotypefq
  file("${pair_id}.trimreport.txt") into trimstat

  script:
  """
  bash $baseDir/process_scripts/preproc_fastq/trimgalore.sh -p ${pair_id} -a ${read1} -b ${read2}
  perl $baseDir/scripts/parse_trimreport.pl ${pair_id}.trimreport.txt *trimming_report.txt
  """
}

process align {
  errorStrategy 'ignore'
  publishDir "$params.output/$famid/$pair_id", mode: 'copy'

  input:
  set famid,pair_id, file(fq1), file(fq2) from trimread

  output:
  set pair_id, file("${pair_id}.bam"),file("${pair_id}.bam.bai") into svbam
  set pair_id, file("${pair_id}.bam"),file("${pair_id}.bam.bai") into cnvbam
  set famid,pair_id, file("${pair_id}.bam") into aligned

  script:
  """
  bash $baseDir/process_scripts/alignment/dnaseqalign.sh -r $index_path -a $params.aligner -p $pair_id -x $fq1 -y $fq2 $alignopts
  """
 }

// process hlacalls {
//   errorStrategy 'ignore'
//   publishDir "$params.output/$pair_id", mode: 'copy'

//   input:
//   set famid,pair_id, file(fq1), file(fq2) from genotypefq

//   output:
//   file("*.hisat_hla.*") into genotype
//   when:
//   params.genome =~ /human/ 
//   script:
//   """
//   bash $baseDir/process_scripts/alignment/hisat_genotype.sh -p $pair_id -x $fq1 -y $fq2
//   """
// }


process markdups {
  publishDir "$params.output/$famid/$pair_id", mode: 'copy'
  input:
  set famid,pair_id, file(sbam) from aligned

  output:
  set famid,pair_id, file("${pair_id}.dedup.bam") into qcbam
  set famid,pair_id, file("${pair_id}.dedup.bam"),file("${pair_id}.dedup.bam.bai") into deduped

  script:
  """
  bash $baseDir/process_scripts/alignment/markdups.sh -a $params.markdups -b $sbam -p $pair_id
  """
}

process seqqc {
  errorStrategy 'ignore'
  publishDir "$params.output/$famid/$pair_id", mode: 'copy'

  input:
  set famid,pair_id, file(sbam) from qcbam

  output:
  file("${pair_id}.flagstat.txt") into alignstats
  file("${pair_id}.ontarget.flagstat.txt") into ontarget
  file("${pair_id}.dedupcov.txt") into dedupcov
  file("${pair_id}.meanmap.txt") into meanmap
  file("${pair_id}.libcomplex.txt") into libcomplex
  file("${pair_id}.hist.txt") into insertsize
  file("${pair_id}.alignmentsummarymetrics.txt") into alignmentsummarymetrics
  file("*fastqc*") into fastqc
  set pair_id, file("${pair_id}.ontarget.bam"),file("${pair_id}.ontarget.bam.bai") into ontargetbam
  set pair_id,file("${pair_id}.ontarget.bam"),file("${pair_id}.ontarget.bam.bai") into genocovbam
  file("${pair_id}.genomecov.txt") into genomecov
  file("${pair_id}.covhist.txt") into covhist
  file("*coverage.txt") into capcovstat
  file("${pair_id}.mapqualcov.txt") into mapqualcov

  script:
  """
  module load samtools/gcc/1.8
  bash $baseDir/process_scripts/alignment/bamqc.sh -c $capture_bed -n dna -d 1 -r $index_path -b $sbam -p $pair_id
  mv ${pair_id}.genomecov.txt ${pair_id}.dedupcov.txt
  mv ${pair_id}.covhist.txt ${pair_id}.covuniqhist.txt
  mv ${pair_id}_lowcoverage.txt ${pair_id}_lowcoverageuniq.txt
  mv ${pair_id}_exoncoverage.txt ${pair_id}_exoncoverageuniq.txt
  bash $baseDir/process_scripts/alignment/bamqc.sh -c $capture_bed -n dna -r $index_path -b $sbam -p $pair_id
  
  """
}

process parse_stat {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  file(txt) from alignstats.toList()
  file(lc) from libcomplex.toList()
  file(is) from insertsize.toList()
  file(gc) from genomecov.toList()
  file(on) from ontarget.toList()
  file(tr) from trimstat.toList()
  file(mq) from mapqualcov.toList()
  file(de) from dedupcov.toList()
  file(mm) from meanmap.toList()
  file(asmet) from alignmentsummarymetrics.toList()
  
  output:
  file('*sequence.stats.txt')
  file('*.png')
  
  script:
  """
  module load R/3.2.1-intel
  perl $baseDir/scripts/parse_seqqc.pl *.genomecov.txt
  perl $baseDir/scripts/covstats.pl *.mapqualcov.txt
  Rscript $baseDir/scripts/plot_hist_genocov.R
  """
}

process gatkbam {
  errorStrategy 'ignore'
  publishDir "$params.output/$famid/$pair_id", mode: 'copy'

  input:
  set famid,pair_id, file(sbam), file(idx) from deduped

  output:
  set famid,file("${pair_id}.final.bam"),file("${pair_id}.final.bai") into gbam
  
  script:
  """
  bash $baseDir/process_scripts/variants/gatkrunner.sh -a gatkbam -b $sbam -r ${index_path} -p $pair_id
  """
}

process pindel {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid", mode: 'copy'
  input:
  set subjid,file(ssbam),file(ssidx) from svbam
  output:
  file("${subjid}.pindel_tandemdup.pass.vcf.gz") into tdvcf
  file("${subjid}.pindel_indel.pass.vcf.gz") into pindelvcf
  file("${subjid}.dna.genefusion.txt") into gf
  script:
  """
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/variants/pindel.sh -r ${index_path} -p ${subjid}
  perl $baseDir/process_scripts/variants/filter_pindel.pl -d ${subjid}.pindel_tandemdup.vcf.gz -s ${subjid}.pindel_sv.vcf.gz -i ${subjid}.pindel_indel.vcf.gz
  module load samtools/gcc/1.8 snpeff/4.3q
  bgzip ${subjid}.pindel_indel.pass.vcf
  bgzip ${subjid}.pindel_tandemdup.pass.vcf
  grep '#CHROM' ${subjid}.pindel_sv.pass.vcf > ${subjid}.dna.genefusion.txt
  cat ${subjid}.pindel_sv.pass.vcf | \$SNPEFF_HOME/scripts/vcfEffOnePerLine.pl |java -jar \$SNPEFF_HOME/SnpSift.jar extractFields - CHROM POS END ANN[*].EFFECT ANN[*].GENE ANN[*].HGVS_C ANN[*].HGVS_P GEN[*] |grep -E 'CHROM|gene_fusion' |uniq >> ${subjid}.dna.genefusion.txt
  """
}

gbam
   .groupTuple(by:0)		
   .into { fbbam; strelkabam; platbam; gkbam }

process freebayes {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/individual_callers", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from fbbam
  output:
  set subjid,file("${subjid}.fb.vcf.gz") into fbvcf
  set subjid,file("${subjid}.fb.ori.vcf.gz") into fbori
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a freebayes
  bash $baseDir/process_scripts/variants/norm_annot.sh -r $index_path -p ${subjid}.fb -v ${subjid}.freebayes.vcf.gz
  mv ${subjid}.freebayes.vcf.gz ${subjid}.fb.ori.vcf.gz
  mv ${subjid}.fb.norm.vcf.gz ${subjid}.fb.vcf.gz
  """
}
process gatk {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/individual_callers", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from gkbam
  output:
  set subjid,file("${subjid}.gatk.vcf.gz") into gatkvcf
  set subjid,file("${subjid}.gatk.ori.vcf.gz") into gatkori
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a gatk
  bash $baseDir/process_scripts/variants/norm_annot.sh -r $index_path -p ${subjid}.gatk -v ${subjid}.gatk.vcf.gz
  mv ${subjid}.gatk.vcf.gz ${subjid}.gatk.ori.vcf.gz
  mv ${subjid}.gatk.norm.vcf.gz ${subjid}.gatk.vcf.gz
  """
}
process strelka {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/individual_callers", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from strelkabam
  output:
  set subjid,file("${subjid}.strelka2.vcf.gz") into strelkavcf
  set subjid,file("${subjid}.strelka2.ori.vcf.gz") into strelkaori
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a strelka2
  bash $baseDir/process_scripts/variants/norm_annot.sh -r $index_path -p ${subjid}.strelka2 -v ${subjid}.strelka2.vcf.gz
  mv ${subjid}.strelka2.vcf.gz ${subjid}.strelka2.ori.vcf.gz
  mv ${subjid}.strelka2.norm.vcf.gz ${subjid}.strelka2.vcf.gz
  """
}
process platypus {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid/individual_callers", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from platbam
  output:
  set subjid,file("${subjid}.platypus.vcf.gz") into platvcf
  set subjid,file("${subjid}.platypus.ori.vcf.gz") into platori
  script:				       
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a platypus
  bash $baseDir/process_scripts/variants/norm_annot.sh -r $index_path -p ${subjid}.platypus -v ${subjid}.platypus.vcf.gz
  mv ${subjid}.platypus.vcf.gz ${subjid}.platypus.ori.vcf.gz
  mv ${subjid}.platypus.norm.vcf.gz ${subjid}.platypus.vcf.gz
  """
}
Channel
  .empty()
  .mix(fbvcf,platvcf,gatkvcf,strelkavcf)
  .groupTuple(by:0)
  .set { vcflist}

process integrate {
  errorStrategy 'ignore'
  publishDir "$params.output/$subjid", mode: 'copy'
  input:
  set subjid,file(vcf) from vcflist
  output:
  file("${subjid}.union.vcf.gz") into unionvcf
  file("${subjid}.annot.vcf.gz") into annotvcf
  script:
  """
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/variants/union.sh -r $index_path -p $subjid
  bash $baseDir/process_scripts/variants/annotvcf.sh -p $subjid -r $index_path -v ${subjid}.union.vcf.gz
  """
}
