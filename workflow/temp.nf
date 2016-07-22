#!/usr/bin/env nextflow

// Default parameter values to run tests
params.genome="/project/shared/bicf_workflow_ref/GRCh38"
capture="$params.genome/gencode.genes.v24.chr.bed"
dbsnp="$params.genome/dbSnp.vcf.gz"
indel="$params.genome/Mills_G1K_indels.b38.vcf.gz"
params.gatkref="/project/apps_database/iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa"
params.alignment=1
params.fastqs="$baseDir/../test_data/*.fastq.gz"
params.pairs="pe"
params.design="$baseDir/../test_data/design.txt"
design_file = file(params.design)

bedgenome = "$params.genome/genomefile.txt"
btoolsgenome = file(bedgenome)

fastqs=file(params.fastqs)
knownindel=file(indel)
gatkref=file(params.gatkref)
// params genome is the directory
// base name for the index is always genome
index_path = file(params.genome)
index_name = "genome"
capture_bed = file(capture)
dbsnp=file(dbsnp)

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
process write_names {
  input:
  set pair_id, file(read1), file(read2) from read_pe
  output:
  set "${read1.baseName.split("\\.fastq")[0]}_val_1.fq.gz" into results
  script:
  """
  ln -s ${read1} "${read1.baseName.split("\\.fastq")[0]}_val_1.fq.gz"
  """
  }

process sample_correlation {    
    
    input:
    file input_files from results.toList()
    
    exec:
    println "${(input_files as List).join(',')}"

}