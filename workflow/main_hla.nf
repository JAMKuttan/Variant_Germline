#!/usr/bin/env nextflow

// Default parameter values to run tests
params.fastqs="$baseDir/../test2/*.fastq.gz"
params.design="$baseDir/../test2/design.txt"

params.genome="/project/shared/bicf_workflow_ref/GRCh38"
params.capture="$params.genome/MedExomeHemOnc.bed"
params.pairs="pe"

dbsnp="$params.genome/dbSnp.vcf.gz"
indel="$params.genome/GoldIndels.vcf.gz"
design_file = file(params.design)

fastqs=file(params.fastqs)
knownindel=file(indel)
gatkref=file("$params.genome/genome.fa")
index_path = file(params.genome)
index_name = "genome"
capture_bed = file(params.capture)
dbsnp=file(dbsnp)

seqprg = file("/project/BICF//BICF_Core/s166458/seqprg/")

snpeff_vers = 'GRCh38.82';
if (params.genome == '/project/shared/bicf_workflow_ref/GRCm38') {
   snpeff_vers = 'GRCm38.82';
}
if (params.genome == '/project/shared/bicf_workflow_ref/GRCh37') {
   snpeff_vers = 'GRCh37.75';
}

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

process trimpe {
  input:
  set pair_id, file(read1), file(read2) from read_pe
  output:
  set pair_id, file("${read1.baseName.split("\\.fastq")[0]}_val_1.fq.gz"), file("${read2.baseName.split("\\.fastq")[0]}_val_2.fq.gz") into trimpe
  file("${pair_id}.trimreport.txt") into trimstatpe
  script:
  """
  module load trimgalore/0.4.1 cutadapt/1.9.1
  trim_galore --paired -q 25 --illumina --gzip --length 35 ${read1} ${read2}
  perl $baseDir/scripts/parse_trimreport.pl ${pair_id}.trimreport.txt ${read1}_trimming_report.txt ${read2}_trimming_report.txt
  """
}
process trimse {
  input:
  set pair_id, file(read1) from read_se
  output:
  set pair_id, file("${read1.baseName.split("\\.fastq.gz")[0]}_trimmed.fq.gz") into trimse
  file("${pair_id}.trimreport.txt") into trimstatse
  script:
  """
  module load trimgalore/0.4.1 cutadapt/1.9.1
  trim_galore -q 25 --illumina --gzip --length 35 ${read1}
  perl $baseDir/scripts/parse_trimreport.pl ${pair_id}.trimreport.txt ${read1}_trimming_report.txt
  """
}

process alignpe {

  publishDir "$baseDir/output", mode: 'copy'
  cpus 32

  input:
  set pair_id, file(fq1), file(fq2) from trimpe
  output:
  set pair_id, file("${pair_id}.bam"), file("${pair_id}.bam.bai") into alignpe
  set pair_id,file("${pair_id}.bam"),file("${pair_id}.discordants.bam"),file("${pair_id}.splitters.bam") into svbam
  file("${pair_id}.hist.txt") into insertsize
  when:
  params.pairs == 'pe'
  script:
  """
  module load seqtk/1.2-r94 bwa/intel/0.7.15 samtools/intel/1.3 speedseq/20160506 picard/1.127
  seqtk mergepe ${fq1} ${fq2} | bwa mem -M -p -t 30 -R '@RG\\tID:${pair_id}\\tLB:tx\\tPL:illumina\\tPU:barcode\\tSM:${pair_id}' ${index_path}/genome.fa - 2> log.bwamem | ${seqprg}/bwa.kit/k8 /cm/shared/apps/bwa/intel/0.7.15/bwakit/bwa-postalt.js -p ${pair_id}.hla ${index_path}/genome.fa.alt | samtools view -1 - > output.unsort.bam
  ${seqprg}/bwa.kit/run-HLA ${pair_id}.hla > ${pair_id}.hla.top 2> ${pair_id}.log.hla
  touch ${pair_id}.hla.HLA-dummy.gt
  cat ${pair_id}.hla.HLA*.gt | grep ^GT | cut -f2- > ${pair_id}.hla.all
  sambamba sort -t 30 -o output.dups.bam output.unsort.bam
  sambamba markdup -t 20 output.dups.bam ${pair_id}.bam
  sambamba view -h output.unsort.bam | samblaster -M -a --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 -d discordants.sam -s splitters.sam > temp.sam
  gawk '{ if (\$0~"^@") { print; next } else { \$10="*"; \$11="*"; print } }' OFS="\\t" splitters.sam | samtools  view -S -b - | samtools sort -o ${pair_id}.splitters.bam -
  gawk '{ if (\$0~"^@") { print; next } else { \$10="*"; \$11="*"; print } }' OFS="\\t" discordants.sam | samtools  view -S  -b - | samtools sort -o ${pair_id}.discordants.bam -
  java -Xmx4g -jar \$PICARD/picard.jar CollectInsertSizeMetrics INPUT=${pair_id}.bam HISTOGRAM_FILE=${pair_id}.hist.ps REFERENCE_SEQUENCE=${index_path}/${index_name}.fa OUTPUT=${pair_id}.hist.txt
  """
}
process alignse {

  //publishDir $outDir, mode: 'copy'
  cpus 32

  input:
  set pair_id, file(fq1) from trimse
  output:
  set pair_id, file("${pair_id}.bam"),file("${pair_id}.bam.bai") into alignse

  when:
  params.pairs == 'se'
  script:
  """
  module load bwa/intel/0.7.15 samtools/intel/1.3 picard/1.127 speedseq/20160506
  bwa mem -M -p -t 30 -R '@RG\tLB:tx\tPL:illumina\tID:${pair_id}\tPU:barcode\tSM:${pair_id}' ${index_path}/genome.fa ${fq1} 2> log.bwamem | ${seqprg}/bwa.kit/k8 /cm/shared/apps/bwa/intel/0.7.15/bwakit/bwa-postalt.js -p ${pair_id}.hla ${index_path}/genome.fa.alt | samtools view -1 - > output.unsort.bam
  ${seqprg}/bwa.kit/run-HLA ${pair_id}.hla > ${pair_id}.hla.top 2> ${pair_id}.log.hla
  touch ${pair_id}.hla.HLA-dummy.gt
  cat ${pair_id}.hla.HLA*.gt | grep ^GT | cut -f2- > ${pair_id}.hla.all
  sambamba sort -t 30 -o output.dups.bam output.unsort.bam
  sambamba markdup -t 20 output.dups.bam ${pair_id}.bam
  """
}

Channel
  .empty()
  .mix(alignse, alignpe)
  .set { aligned }

Channel
  .empty()
  .mix(trimstatse, trimstatpe)
  .set { trimstat }

// Calculate Metrics of Quality of Alignment
process seqqc {

  publishDir "$baseDir/output", mode: 'copy'

  input:
  set pair_id, file(sbam),file(idx) from aligned
  output:
  file("${pair_id}.flagstat.txt") into alignstats
  file("${pair_id}.libcomplex.txt") into libcomplex
  file("${pair_id}.ontarget.flagstat.txt") into ontarget
  set pair_id, file("${pair_id}.ontarget.bam"),file("${pair_id}.ontarget.bam.bai") into targetbam
  file("${pair_id}.ontarget.bam") into ssbam
  file("${pair_id}.ontarget.bam.bai") into ssidx
  file("${pair_id}.genomecov.txt") into genomecov
  set file("${pair_id}_fastqc.zip"),file("${pair_id}_fastqc.html") into fastqc

  script:
  """
  module load bedtools/2.25.0 picard/1.127 samtools/intel/1.3 fastqc/0.11.2 speedseq/20160506
  fastqc -f bam ${sbam}
  sambamba flagstat -t 30 ${sbam} > ${pair_id}.flagstat.txt
  sambamba view -t 30 -f bam -L  ${capture_bed} -o ${pair_id}.ontarget.bam ${sbam}
  sambamba flagstat -t 30 ${pair_id}.ontarget.bam > ${pair_id}.ontarget.flagstat.txt
  bedtools coverage -sorted -hist -g ${index_path}/genomefile.txt -b ${pair_id}.ontarget.bam -a ${capture_bed} | grep ^all >  ${pair_id}.genomecov.txt
  perl $baseDir/scripts/subset_bam.pl ${pair_id}.ontarget.bam ${pair_id}.ontarget.flagstat.txt
  java -Xmx8g -jar \$PICARD/picard.jar EstimateLibraryComplexity INPUT=${pair_id}.subset50M.bam OUTPUT=${pair_id}.libcomplex.txt
  """
}

process parse_stat {

  publishDir "$baseDir/output", mode: 'copy'

  input:
  file(txt) from alignstats.toList()
  file(lc) from libcomplex.toList()
  file(is) from insertsize.toList()
  file(gc) from genomecov.toList()
  file(on) from ontarget.toList()
  file(tr) from trimstat.toList()
  
  output:
  file('sequence.stats.txt')

  script:
  """
  perl $baseDir/scripts/parse_seqqc.pl *.libcomplex.txt
  """
}

process svcall {

  publishDir "$baseDir/output", mode: 'copy'
  cpus 32
  input:
  set pair_id,file(ssbam),file(discordbam),file(splitters) from svbam
  output:
  
  script:
  """
  module load samtools/intel/1.3 bedtools/2.25.0 bcftools/intel/1.3 snpeff/4.2 speedseq/20160506 vcftools/0.1.14
  speedseq sv -t 30 -o ${pair_id}.sssv -R ${index_path}/${index_name}.fa -B ${ssbam} -D ${discordbam} -S ${splitters}
  perl $baseDir/scripts/parse_svresults.pl -i ${pair_id}.sssv.sv.vcf.gz -r ${index_path}
  """
}

process gatkbam {

  publishDir "$baseDir/output", mode: 'copy'

  input:
  set pair_id, file(dbam), file(idx) from targetbam

  output:
  set pair_id,file("${pair_id}.final.bam"),file("${pair_id}.final.bai") into gatkbam
  file("${pair_id}.final.bam") into sambam
  file("${pair_id}.final.bai") into samidx
  file("${pair_id}.final.bam") into platbam
  file("${pair_id}.final.bai") into platidx
  script:
  """
  module load gatk/3.5 samtools/intel/1.3
  samtools index ${dbam}
  java -Xmx32g -jar \$GATK_JAR -T RealignerTargetCreator -known ${knownindel} -R ${gatkref} -o ${pair_id}.bam.list -I ${dbam} -nt 30 -nct 1
  java -Xmx32g -jar \$GATK_JAR -I ${dbam} -R ${gatkref} --filter_mismatching_base_and_quals -T IndelRealigner -targetIntervals ${pair_id}.bam.list -o ${pair_id}.realigned.bam
  java -Xmx32g -jar \$GATK_JAR -l INFO -R ${gatkref} --knownSites ${dbsnp} -I ${pair_id}.realigned.bam -T BaseRecalibrator -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o ${pair_id}.recal_data.grp -nt 1 -nct 30
  java -Xmx32g -jar \$GATK_JAR -T PrintReads -R ${gatkref} -I ${pair_id}.realigned.bam -BQSR ${pair_id}.recal_data.grp -o ${pair_id}.final.bam -nt 1 -nct 8
  """
}

process gatkgvcf {

  publishDir "$baseDir/output", mode: 'copy'
  cpus 30

  input:
  set pair_id,file(gbam),file(gidx) from gatkbam
  output:
  file("${pair_id}.gatk.g.vcf") into gatkgvcf
  script:
  """
  module load gatk/3.5 bedtools/2.25.0 snpeff/4.2 vcftools/0.1.14 
  java -Xmx32g -jar \$GATK_JAR -R ${gatkref} -D ${dbsnp} -T HaplotypeCaller -stand_call_conf 30 -stand_emit_conf 10.0 -A FisherStrand -A QualByDepth -A VariantType -A DepthPerAlleleBySample -A HaplotypeScore -A AlleleBalance -I ${gbam} -o ${pair_id}.gatk.g.vcf -nct 2
  """
}
process gatk {

  publishDir "$baseDir/output", mode: 'copy'
  cpus 4

  input:
  file(gvcf) from gatkgvcf.toList()

  output:
  file("final.gatkpanel.vcf.gz") into gatkvcf
  script:
  """
  module load python/2.7.x-anaconda gatk/3.5 bedtools/2.25.0 snpeff/4.2 vcftools/0.1.14
  java -Xmx32g -jar \$GATK_JAR -R ${gatkref} -D ${dbsnp} -T GenotypeGVCFs -o final.gatk.vcf -nt 4 --variant ${(gvcf as List).join(' --variant ')}
  java -Xmx32g -jar \$GATK_JAR -R ${gatkref} -T VariantFiltration -V final.gatk.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o final.filtgatk.vcf 
  vcf-annotate -n --fill-type final.filtgatk.vcf | java -jar \$SNPEFF_HOME/SnpSift.jar filter '((QUAL >= 10) & (MQ > 40) & (DP >= 10))' |bedtools intersect -header -a stdin -b ${capture_bed} |bgzip > final.gatkpanel.vcf.gz
  """
}

process mpileup {

  publishDir "$baseDir/output", mode: 'copy'
  cpus 32

  input:
  file(gbam) from sambam.toList()
  file(gidx) from samidx.toList()
  output:
  file("final.sampanel.vcf.gz") into samvcf
  script:
  """
  module load python/2.7.x-anaconda samtools/intel/1.3 bedtools/2.25.0 bcftools/intel/1.3 snpeff/4.2 vcftools/0.1.14
  ls *.bam > final.bam.list
  cut -f 1 ${index_path}/genomefile.chr.txt | xargs -I {} -n 1 -P 32 sh -c "samtools mpileup -t 'AD,ADF,ADR,INFO/AD,SP' -ug -Q20 -C50 -f ${index_path}/${index_name}.fa -b final.bam.list -r {} | bcftools call --format-fields gq,gp -vmO z -o final.sam.{}.vcf.gz"
  vcf-concat final.sam.*.vcf.gz | vcf-sort |vcf-annotate -n --fill-type | java -jar \$SNPEFF_HOME/SnpSift.jar filter '((QUAL >= 10) & (MQ >= 20) & (DP >= 10))' |bedtools intersect -header -a stdin -b ${capture_bed} |bgzip > final.sampanel.vcf.gz
  """
}

process speedseq {

  publishDir "$baseDir/output", mode: 'copy'
  cpus 32
  input:
  file(gbam) from ssbam.toList()
  file(gidx) from ssidx.toList()
  output:
  file("final.sspanel.vcf.gz") into ssvcf
  file("final.complex.vcf.gz") into complexvcf
  script:
  """
  module load python/2.7.x-anaconda samtools/intel/1.3 bedtools/2.25.0 bcftools/intel/1.3 snpeff/4.2 speedseq/20160506 vcftools/0.1.14
  speedseq var -q 10 -w ${capture_bed} -t 32 -o final.ssvar ${index_path}/${index_name}.fa ${gbam}
  vcf-annotate -n --fill-type -n final.ssvar.vcf.gz | java -jar \$SNPEFF_HOME/SnpSift.jar filter '((QUAL >= 10) & (DP >= 10))' |bgzip > final.sspanel.vcf.gz
  java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar filter "(TYPE='complex')" final.sspanel.vcf.gz |bgzip> final.complex.vcf.gz
  """
}
process platypus {

  publishDir "$baseDir/output", mode: 'copy'
  cpus 32

  input:
  file(gbam) from platbam.toList()
  file(gidx) from platidx.toList()

  output:
  file("final.platpanel.vcf.gz") into platvcf

  script:
  """
  module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 platypus/gcc/0.8.1 bcftools/intel/1.3 samtools/intel/1.3 vcftools/0.1.14
  Platypus.py callVariants --minMapQual=10 --mergeClusteredVariants=1 --nCPU=30 --bamFiles="${(gbam as List).join(',')}" --refFile=${index_path}/${index_name}.fa --output=final.platypus.vcf
  vcf-annotate -n --fill-type -n final.platypus.vcf | java -jar \$SNPEFF_HOME/SnpSift.jar filter "((QUAL >= 10) & (QD > 2) & (FILTER = 'PASS'))" |bedtools intersect -header -a stdin -b ${capture_bed} |bgzip > final.platpanel.vcf.gz
  """
}

process integrate {
  publishDir "$baseDir/output", mode: 'copy'

  input:
  file(pvcf) from platvcf
  file(ssvcf)from ssvcf
  file(samvcf) from samvcf
  file(gvcf) from gatkvcf
  file(complex) from complexvcf

  output:
  file("final.integrated.vcf.gz") into finalvcf
  script:
  """
  module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 platypus/gcc/0.8.1  bcftools/intel/1.3 samtools/intel/1.3 jags/4.2.0
  module load vcftools/0.1.14
  zcat ${ssvcf} | bedtools intersect -v -header -a stdin -b ${complex} | bgzip > ss.shuff.vcf.gz
  vcf-shuffle-cols -t ${ssvcf} ${pvcf} | vcf-sort |bedtools intersect -v -header -a stdin -b ${complex} | bgzip > plat.shuff.vcf.gz
  vcf-shuffle-cols -t ${ssvcf} ${gvcf} | vcf-sort | bedtools intersect -v -header -a stdin -b ${complex} | bgzip > gatk.shuff.vcf.gz
  vcf-shuffle-cols -t ${ssvcf} ${samvcf} | bedtools intersect -v -header -a stdin -b ${complex} | bgzip > sam.shuff.vcf.gz
  tabix ss.shuff.vcf.gz
  tabix plat.shuff.vcf.gz
  tabix gatk.shuff.vcf.gz
  tabix sam.shuff.vcf.gz
  vcf-compare ss.shuff.vcf.gz sam.shuff.vcf.gz gatk.shuff.vcf.gz plat.shuff.vcf.gz > vcf_compare.out
  vcf-isec -f --prefix integrate ss.shuff.vcf.gz sam.shuff.vcf.gz gatk.shuff.vcf.gz plat.shuff.vcf.gz
  perl $baseDir/scripts/baysic_blc.pl -c ${complex} -e $baseDir -f ${index_path}/${index_name}.fa ss.shuff.vcf.gz sam.shuff.vcf.gz gatk.shuff.vcf.gz plat.shuff.vcf.gz
  
  """
}
process annotate {
  publishDir "$baseDir/output", mode: 'copy'

  input:
  file(intvcf) from finalvcf

  output:
  file("annot.vcf.gz") into annotvcf
  script:
  if (params.genome == '/project/shared/bicf_workflow_ref/GRCh38')
  """
  module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 platypus/gcc/0.8.1  jags/4.2.0
  java -Xmx10g -jar \$SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c \$SNPEFF_HOME/snpEff.config ${snpeff_vers} final.integrated.vcf.gz | java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar annotate ${index_path}/dbSnp.vcf.gz -  | java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar annotate ${index_path}/clinvar.vcf.gz - | java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar annotate ${index_path}/ExAC.vcf.gz - | java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar annotate ${index_path}/cosmic.vcf.gz - | java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar dbnsfp -v -db ${index_path}/dbNSFP.txt.gz - | java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar gwasCat -db ${index_path}/gwas_catalog.tsv - |bgzip > annot.vcf.gz
  """
  else
  """
  java -Xmx10g -jar \$SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c \$SNPEFF_HOME/snpEff.config ${snpeff_vers} final.integrated.vcf.gz |bgzip > annot.vcf.gz
  """
}
