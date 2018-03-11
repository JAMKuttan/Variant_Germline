#!/usr/bin/env nextflow

params.input = "$baseDir"
params.fastqs="$params.input/*.fastq.gz"
params.design="$params.input/design.txt"
params.output = "$baseDir/output"

params.fastqs="$params.input/*.fastq.gz"
params.design="$params.input/design.txt"

params.genome="/project/shared/bicf_workflow_ref/GRCh38"
params.capture="$params.genome/UTSWV2.bed"
params.pairs="pe"
params.markdups='sambamba'
params.cancer='skip'

reffa=file("$params.genome/genome.fa")
dbsnp="$params.genome/dbSnp.vcf.gz"
indel="$params.genome/GoldIndels.vcf.gz"

fastqs=file(params.fastqs)
design_file = file(params.design)
dbsnp=file(dbsnp)
knownindel=file(indel)
index_path = file(params.genome)
capture_bed = file(params.capture)


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
    familyidx = header.findIndexOf{it == 'FamilyID'};
    oneidx = header.findIndexOf{it == 'FqR1'};
    twoidx = header.findIndexOf{it == 'FqR2'};
    if (twoidx == -1) {
       twoidx = oneidx
       }
   if (familyidx == -1) {
       familyidx = prefixidx
       }    
    while (line = reader.readLine()) {
    	   def row = line.split("\t")
    if (fileMap.get(row[oneidx]) != null) {
	prefix << tuple(row[familyidx],row[prefixidx],fileMap.get(row[oneidx]),fileMap.get(row[twoidx]))
	   }

}
}

if( ! prefix) { error "Didn't match any input files with entries in the design file" }

Channel
  .from(prefix)
  .set { read }

process trimpe {
  errorStrategy 'ignore'
  publishDir "$baseDir/output", mode: 'copy'
  input:
  set famid,pair_id, file(read1), file(read2) from read
  output:
  set famid,pair_id, file("${pair_id}.trim.R1.fastq.gz"),file("${pair_id}.trim.R2.fastq.gz") into trimread
  set famid,pair_id, file("${pair_id}.trim.R1.fastq.gz"),file("${pair_id}.trim.R2.fastq.gz") into fusionfq
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
  """
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/alignment/dnaseqalign.sh -r $index_path -p $pair_id -x $fq1 -y $fq2
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
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/alignment/markdups.sh -a picard_umi -b $sbam -p $pair_id
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
  source /etc/profile.d/modules.sh
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
  source /etc/profile.d/modules.sh
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
  source /etc/profile.d/modules.sh
  module load gatk/3.7 samtools/1.4.1
  bash $baseDir/process_scripts/variants/gatkrunner.sh -a gatkbam -b $sbam -r ${index_path} -p $pair_id
  """
}

gbam
   .groupTuple(by:0)		
   .into { ssbam; sambam; hsbam; strelkabam; platbam }

 
process svcall {
  errorStrategy 'ignore'
  publishDir "$params.output/$pair_id", mode: 'copy'
  input:
  set pair_id,file(ssbam),file(ssidx) from svbam
  output:
  file("${pair_id}.delly.vcf.gz") into dellyvcf
  file("${pair_id}.sssv.sv.vcf.gz") into svvcf
  file("${pair_id}.sv.vcf.gz") into svintvcf
  file("${pair_id}.sv.annot.txt") into svannot
  script:
  """
  source /etc/profile.d/modules.sh	
  bash $baseDir/process_scripts/variants/svcalling.sh -b $ssbam -r $index_path -p $pair_id
  """
}

process mpileup {
  errorStrategy 'ignore'
  //publishDir "$baseDir/output", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from sambam
  
  output:
  set subjid,file("${subjid}.sam.vcf.gz") into samvcf
  script:
  """
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a mpileup
  """
}
process hotspot {
  errorStrategy 'ignore'
  //publishDir "$baseDir/output", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from hsbam
  output:
  set subjid,file("${subjid}.hotspot.vcf.gz") into hsvcf
  when:
  params.cancer == "detect"
  script:
  """
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a hotspot
  """
}
process speedseq {
  errorStrategy 'ignore'
  //publishDir "$baseDir/output", mode: 'copy'
  input:
  set subjid,file(gbam),file(gidx) from ssbam
  output:
  set subjid,file("${subjid}.ssvar.vcf.gz") into ssvcf
  script:
  """
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a speedseq
  """
}
process strelka2 {
  errorStrategy 'ignore'

  input:
  set subjid,file(gbam),file(gidx) from strelkabam
  output:
  set subjid,file("${subjid}.strelka2.vcf.gz") into strelkavcf
  script:
  """
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a strelka2
  """
}
process platypus {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'

  input:
  set subjid,file(gbam),file(gidx) from platbam

  output:
  set subjid,file("${subjid}.platypus.vcf.gz") into platvcf
  script:
  """
  source /etc/profile.d/modules.sh
  bash $baseDir/process_scripts/variants/germline_vc.sh -r $index_path -p $subjid -a platypus
  
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
  //publishDir "$params.output", mode: 'copy'
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
