#!/usr/bin/env nextflow

params.input = "$baseDir"
params.fastqs="$params.input/*.fastq.gz"
params.design="$params.input/design.txt"
params.output = "$baseDir/output"

params.genome="/project/shared/bicf_workflow_ref/human/GRCh38"
params.capture="$params.genome/UTSWV2.bed"
params.pairs="pe"
params.cancer='skip'
params.callers = 'all'
params.markdups='sambamba'
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
  publishDir "$params.output", mode: 'copy'

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
  publishDir "$params.output/$pair_id", mode: 'copy'

  input:
  set famid,pair_id, file(fq1), file(fq2) from trimread

  output:
  set pair_id, file("${pair_id}.bam"),file("${pair_id}.bam.bai") into svbam
  set pair_id, file("${pair_id}.bam"),file("${pair_id}.bam.bai") into cnvbam
  set famid,pair_id, file("${pair_id}.bam") into aligned

  script:
  """
  bash $baseDir/process_scripts/alignment/dnaseqalign.sh -r $index_path -p $pair_id -x $fq1 -y $fq2
  """
 }

process hlacalls {
  errorStrategy 'ignore'
  publishDir "$params.output/$pair_id", mode: 'copy'

  input:
  set famid,pair_id, file(fq1), file(fq2) from genotypefq

  output:
  file("*.hisat_hla.*") into genotype
  when:
  params.genome =~ /human/ 
  script:
  """
  bash $baseDir/process_scripts/alignment/hisat_genotype.sh -p $pair_id -x $fq1 -y $fq2
  """
}


process markdups {
  //publishDir "$params.output/$pair_id", mode: 'copy'

  input:
  set famid,pair_id, file(sbam) from aligned

  output:
  set pair_id, file("${pair_id}.dedup.bam") into qcbam
  set famid,pair_id, file("${pair_id}.dedup.bam") into deduped

  script:
  """
  bash $baseDir/process_scripts/alignment/markdups.sh -a $params.markdups -b $sbam -p $pair_id
  """
}

process seqqc {
  errorStrategy 'ignore'
  publishDir "$params.output/$pair_id", mode: 'copy'

  input:
  set pair_id, file(sbam) from qcbam

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
  publishDir "$params.output/$pair_id", mode: 'copy'

  input:
  set famid,pair_id, file(sbam), file(idx) from deduped

  output:
  set famid,file("${pair_id}.final.bam"),file("${pair_id}.final.bai") into gbam
  
  script:
  """
  bash $baseDir/process_scripts/variants/gatkrunner.sh -a gatkbam -b $sbam -r ${index_path} -p $pair_id
  """
}

gbam
   .groupTuple(by:0)		
   .into { ssbam; sambam; hsbam; strelkabam; platbam; gkbam }

 
process svcall {
  errorStrategy 'ignore'
  publishDir "$params.output/$pair_id", mode: 'copy'
  input:
  set pair_id,file(ssbam),file(ssidx) from svbam
  output:
  file("${pair_id}.sv.vcf.gz") into svvcf

  when:
  params.callsvs != "skip"
  
  script:
  """
  bash $baseDir/process_scripts/variants/svcalling.sh -a $params.callsvs -b $ssbam -r $index_path -p $pair_id
  """
}

process mpileup {
  errorStrategy 'ignore'
  publishDir "$baseDir/output", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from sambam
  
  output:
  set subjid,file("${subjid}.sam.vcf.gz") into samvcf
  when:
  params.callers == 'all' || params.callers == 'mpileup'
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a mpileup
  """
}
process hotspot {
  errorStrategy 'ignore'
  publishDir "$baseDir/output", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from hsbam
  output:
  set subjid,file("${subjid}.hotspot.vcf.gz") into hsvcf
  when:
  params.cancer == "detect"
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a hotspot
  """
}
process speedseq {
  errorStrategy 'ignore'
  publishDir "$baseDir/output", mode: 'copy'
  input:
  set subjid,file(gbam),file(gidx) from ssbam
  output:
  set subjid,file("${subjid}.ssvar.vcf.gz") into ssvcf
  when:
  params.callers == 'all' || params.callers == 'speedseq'
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a speedseq
  """
}
process strelka2 {
  errorStrategy 'ignore'
  publishDir "$baseDir/output", mode: 'copy'
  input:
  set subjid,file(gbam),file(gidx) from strelkabam
  output:
  set subjid,file("${subjid}.strelka2.vcf.gz") into strelkavcf
  when:
  params.callers == 'all' || params.callers == 'strelka2'
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a strelka2
  """
}
process platypus {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from platbam

  output:
  set subjid,file("${subjid}.platypus.vcf.gz") into platvcf
  when:
  params.callers == 'all' || params.callers == 'platypus'
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a platypus
  """
}

process gatk {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from gkbam

  output:
  set subjid,file("${subjid}.platypus.vcf.gz") into gatkvcf
  when:
  params.callers == 'gatk'
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a gatk
  """
}

if (params.cancer == "detect") {
   Channel
	.empty()
  	.mix(ssvcf,strelkavcf,samvcf,platvcf,hsvcf)
	.groupTuple(by:0)
	.into { vcflist}
}
else {
   Channel
	.empty()
  	.mix(ssvcf,strelkavcf,samvcf,platvcf)
	.groupTuple(by:0)
	.into { vcflist}

}


process integrate {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'
  input:
  set subjid,file(vcfs) from vcflist
    
  output:
  file("${subjid}.annot.vcf.gz") into annotvcf
  script:
  """
  bash $baseDir/process_scripts/variants/union.sh -r $index_path -p $subjid
  bash $baseDir/process_scripts/variants/annotvcf.sh -p $subjid -r $index_path -v ${subjid}.union.vcf.gz
  """
}
