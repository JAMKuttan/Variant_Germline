#!/usr/bin/env nextflow

// Default parameter values to run tests
params.genome="/project/apps_database/iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/"
params.capture="/project/BICF/BICF_Core/shared/refdata/GRCh38/gencode.genes.v24.chr.bed"
params.dbsnp="/project/apps_database/dbSNP/organisms/human_9606/VCF/GATK/All_20160407.vcf.gz"
params.indel="/project/BICF/BICF_Core/shared/refdata/Mills_G1K_indels.b38.vcf.gz"
params.alignment=1
params.fastqs="$baseDir/../test_data/*.fastq.gz"
params.pairs="pe"
params.design="$baseDir/../test_data/design.txt"
design_file = file(params.design)

fastqs=file(params.fastqs)
targetbed=file(params.capture)
dbsnp=file(params.dbsnp)
knownindel=file(params.indel)

// params genome is the directory
// base name for the index is always genome
index_path = file(params.genome)
index_name = "genome"


// Pair handling, helper function taken from rnatoy
// which is covered by the GNU General Public License v3
// https://github.com/nextflow-io/rnatoy/blob/master/main.nf
def fileMap = [:]

fastqs.each {
    final fileName = it.getFileName().toString()
    prefix = fileName.lastIndexOf('/')
    fileMap[fileName] = it
}
def prefix = []
new File(params.design).withReader { reader ->
    def hline = reader.readLine()
    def header = hline.split("\t")
    prefixidx = header.findIndexOf{it == 'SampleID'};
    oneidx = header.findIndexOf{it == 'FullPathToFqR1'};
    twoidx = header.findIndexOf{it == 'FullPathToFqR2'};
    if (twoidx == -1) {
       twoidx = oneidx
       }      
    while (line = reader.readLine()) {
    	   def row = line.split("\t")
	   if (fileMap.get(row[oneidx]) != null) {
	      prefix << tuple(row[prefixidx],fileMap.get(row[oneidx]),fileMap.get(row[twoidx]))
	   }
	  
} 
}

if( ! prefix) { error "Didn't match any input files with entries in the design file" }

if (params.pairs == 'pe') {
Channel
  .from(prefix)
  .set { read_pe }
Channel
  .empty()
  .set { read_se } 
}
if (params.pairs == 'se') {
Channel
  .from(prefix)
  .into { read_se }
Channel
  .empty()
  .set { read_pe }
}

//
// Trim raw reads using trimgalore
//
process trimpe {
  input:
  set pair_id, file(read1), file(read2) from read_pe
  output:
  set pair_id, file("${read1.baseName.split("\\.", 2)[0]}_val_1.fq.gz"), file("${read2.baseName.split("\\.", 2)[0]}_val_2.fq.gz") into trimpe
  script:
  """
  module load trimgalore/0.4.1 cutadapt/1.9.1
  trim_galore --paired -q 25 --illumina --gzip --length 35 --no_report_file ${read1} ${read2}
  """
}
process trimse {
  input:
  set pair_id, file(read1) from read_se
  output:
  set pair_id, file("${read1.baseName.split("\\.", 2)[0]}_trimmed.fq.gz") into trimse
  script:
  """
  module load trimgalore/0.4.1 cutadapt/1.9.1
  trim_galore -q 25 --illumina --gzip --length 35 --no_report_file ${read1}
  """
}
//
// Align trimmed reads to genome index with bwa
// Sort and index with samtools
// QC aligned reads with fastqc
// Alignment stats with samtools
//
process alignpe {

  //publishDir $outDir, mode: 'copy'
  cpus 4

  input:
  set pair_id, file(fq1), file(fq2) from trimpe
  output:
  set pair_id, file("${pair_id}.bam") into alignpe
  set file("${pair_id}.bam") into ssbampe
  set file("${pair_id}.discordants.bam") into disbam
  set file("${pair_id}.splitters.bam") into splitbam
  
  when:
  params.pairs == 'pe'
  params.alignment == 1
  script:
  """
  module load speedseq/20160506
  speedseq align -R '@RG\tLB:tx\tPL:illumina\tID:${pair_id}\tPU:barcode\tSM:${pair_id}' -o ${pair_id} -t 32 ${index_path}/${index_name}.fa ${fq1} ${fq2}
  java -Xmx4g -jar \$PICARD/picard.jar CollectInsertSizeMetrics INPUT=${pair_id}.bam HISTOGRAM_FILE=${pair_id}.hist.ps REFERENCE_SEQUENCE=${index_path}/${index_name}.fa OUTPUT=${pair_id}.hist.txt
  """
}
process alignse {

  //publishDir $outDir, mode: 'copy'
  cpus 4

  input:
  set pair_id, file(fq1) from trimse
  output:
  set pair_id, file("${pair_id}.bam") into alignse
  set file("${pair_id}.bam") into ssbamse
  when:
  params.pairs == 'se'
  params.alignment == 1
  script:
  """
  module load bwa/intel/0.7.12 samtools/intel/1.3 picard/1.127
  bwa mem -M -R '@RG\tLB:tx\tPL:illumina\tID:${pair_id}\tPU:barcode\tSM:${pair_id}' -t 32 ${index_path}/${index_name}.fa ${fq1} > output.sam
  samtools view -b -u -S -o output.unsort.bam output.sam
  samtools sort -o output.dups.bam output.unsort.bam
  java -Xmx4g -jar \$PICARD/picard.jar MarkDuplicates INPUT=output.dups.bam REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true METRICS_FILE=${pair_id}.dups OUTPUT=${pair_id}.bam
  """
}

Channel
  .empty()
  .mix(alignse, alignpe)
  .tap { aligned2 }
  .set { aligned }

Channel
  .empty()
  .mix(ssbamse, ssbampe)
  .set { ssbam }

//
// Calculate Metrics of Quality of Alignment
//
process seqqc {

  memory '4GB'
  publishDir "$baseDir/output", mode: 'copy'

  input:
  set pair_id, file(sbam) from aligned2
  output:
  file("${pair_id}.flagstat.txt") into alignstats
  file("${pair_id}.libcomplex.txt") into libcomplex
  file("${pair_id}.hist.txt") into insertsize
  file("${pair_id}.genomecov.txt") into genomecov
  script:
  """
  module load bedtools/2.25.0 picard/1.127 samtools/intel/1.3 fastqc/0.11.2 
  fastqc -f bam ${sbam}
  samtools flagstat ${sbam} > ${pair_id}.flagstat.txt
  java -Xmx4g -jar \$PICARD/picard.jar EstimateLibraryComplexity INPUT=${sbam} OUTPUT=${pair_id}.libcomplex.txt
  bedtools coverage -sorted -hist -g ${index_path}/${index_name}.genomefile.txt -b ${sbam} -a ${capture_bed} | grep ^all >  ${pair_id}.genomecov.txt
  """
}

// //
// // Summarize all flagstat output
// //
process parse_stat {

  publishDir "$baseDir/output", mode: 'copy'

  input:
  file(txt) from alignstats.toList()
  file(lc) from libcomplex.toList()
  file(is) from insertsize.toList()
  file(gc) from genomecov.toList()

  output:
  file('sequence.stats.txt')
  script:
  """
  perl $baseDir/scripts/parse_seqqc.pl *.flagstat.txt
  """
}

// //
// // Read summarization with subread
// //
process gatkbam {

  publishDir "$baseDir/output", mode: 'copy'
  cpus 4

  input:
  set pair_id, file(dbam) from aligned

  output:
  set file("${pair_id}.final.bam") into gatkbam
  set file("${pair_id}.final.bam") into sambam
  set file("${pair_id}.final.bam") into platbam
  
  """
  module load module subread/1.5.0-intel gatk/3.3-0
  java -Xmx4g -jar $GATK_JAR -T RealignerTargetCreator -known ${knownindel} -R ${index_path}/${index_name} -o ${pair_id}.bam.list -I ${dbam}
  java -Xmx4g -jar $GATK_JAR -I ${dbam} -R ${index_path}/${index_name}.fa --filter_mismatching_base_and_quals -T IndelRealigner -targetIntervals ${pair_id}.bam.list -o ${pair_id}.realigned.bam
  java -Xmx4g -jar $GATK_JAR -l INFO -R ${index_path}/${index_name} --knownSites ${dbsnp} -I ${pair_id}.realigned.bam -T BaseRecalibrator -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o ${pair_id}.recal_data.grp
  java -Xmx4g -jar $GATK_JAR -T PrintReads -R ${index_path}/${index_name} -I ${pair_id}.realigned.bam -BQSR ${pair_id}.recal_data.grp -o ${pair_id}.final.bam
  """
}

// //
// // Run GATK
// //
process gatk {

  publishDir "$baseDir/output", mode: 'copy'
  cpus 32

  input:
  set file(gbam) from gatkbam.toList()
  output:
  set file("final.gatkpanel.vcf.gz") into gatkvcf
     script:
  """
  module load module gatk/3.3-0 vcftools/0.1.11 bedtools/2.25.0 samtools/intel/1.3

  java -jar $GATK_JAR -R ${index_path}/${index_name}.fa -D ${dbsnp} -T HaplotypeCaller -stand_call_conf 30 -stand_emit_conf 10.0 -A FisherStrand -A QualByDepth -A VariantType -A DepthPerAlleleBySample -A HaplotypeScore -A AlleleBalance -I ${gbam.join(' -I ')}  -o final.gatk.vcf -nt 1 -nct 32
  vcf-annotate -n --fill-type final.gatk.vcf | java -jar $SNPEFF_HOME/SnpSift.jar filter '((QUAL >= 10) & (QD > 2) & (FS <= 60) & (MQ > 40) & (DP > 10))' | bedtools intersect -header -a stdin -b ${targetbed} |bgzip > final.gatkpanel.vcf};

  """
}

process mpileup {

  publishDir "$baseDir/output", mode: 'copy'
  cpus 4

  input:
  set file(gbam) from sambam.toList()
  output:
  set file("final.sampanel.vcf.gz") into samvcf
  script:
  """
  module load module samtools/intel/1.3 vcftools/0.1.11 bedtools/2.25.0 bcftools/intel/1.3
  ls *.bam > final.bam.list  
  samtools mpileup -t 'AD,ADF,ADR,INFO/AD,SP' -ug -Q20 -C50 -f ${index_path}/${index_name}.fa -b final.bam.list | bcftools call --format-fields gq,gp -vmO z -o final.sam.vcf.gz
  vcf-annotate -n --fill-type final.sam.vcf.gz ate -n --fill-type final.sam.vcf.gz | java -jar $SNPEFF_HOME/SnpSift.jar filter '((QUAL >= 10) & (MQ >= 20) & (DP >= 10))' |bedtools intersect -header -a stdin -b ${targetbed} |bgzip > final.sampanel.vcf.gz
  """
}

// process printHello {
//    input:
//    file r1
//    file r2
//    output:
//    stdout into result
//    """
//    echo ${index_path} ${index_name}
//    """
// }
// result.println()
