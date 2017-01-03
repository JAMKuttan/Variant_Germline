#!/usr/bin/perl -w
#integration_vcf_snpindel_highcov.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %opt = ();
my $results = GetOptions (\%opt,'refdir|r=s','gver|g=s','outdir|o=s');

$fasta = $opt{refdir}."/genome.fa";

my @vcffiles = `ls *.vcf.gz`;
chomp(@vcffiles);
my %vcfs;
foreach $vcf (@vcffiles) {
    my ($prefix,$caller,$ext) = split(/\./,$vcf);
    $vcfs{$prefix}{$caller} = $vcf;
}
foreach $sid (keys %vcfs) {
    open SH, ">$sid\_integrate.sh" or die $!;
    print SH "#!/bin/bash\n#SBATCH --job-name varcall_bysamp\n#SBATCH -N 1\n#SBATCH -t 1-0:0:00\n";
    print SH "#SBATCH -o job_%j.out\n#SBATCH -e job_%j.err\n";
    
    print SH "module load gatk/3.5 python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 platypus/gcc/0.8.1  bcftools/intel/1.3 samtools/intel/1.3 jags/4.2.0\n";
    print SH "module load vcftools/0.1.14\n";
    print SH "cd ",$opt{outdir}."/".$sid,"\n";
    print SH "tabix ",$vcfs{$sid}{sam},"\n";
    print SH "tabix ",$vcfs{$sid}{platypus},"\n";
    print SH "tabix ",$vcfs{$sid}{gatk},"\n";
    print SH "tabix ",$vcfs{$sid}{ssvar},"\n";
    print SH "tabix ",$vcfs{$sid}{hotspot},"\n";
    print SH "java -Xmx32g -jar \$GATK_JAR -R $fasta -T CombineVariants --variant:gatk ".$vcfs{$sid}{gatk}." --variant:mpileup ".$vcfs{$sid}{sam}." --variant:samhotspot ".$vcfs{$sid}{hotspot}." --variant:freebayes ".$vcfs{$sid}{ssvar}." --variant:platypus ".$vcfs{$sid}{platypus}." -genotypeMergeOptions PRIORITIZE -priority mpileup,gatk,freebayes,platypus,samhotspot -o $sid\.int.vcf\n";
    print SH "bedtools multiinter -i ".join(" ",$vcfs{$sid}{gatk},$vcfs{$sid}{sam},$vcfs{$sid}{ssvar},$vcfs{$sid}{platypus},$vcfs{$sid}{hotspot})." -names gatk sam ssvar platypus hotspot |cut -f 1,2,3,5 | bedtools sort -i stdin >  $sid\_integrate.bed\n";
    print SH "bgzip $sid\_integrate.bed\n";
    print SH "tabix $sid\_integrate.bed.gz\n";
    print SH "bcftools annotate -a $sid\_integrate.bed.gz --columns CHROM,FROM,TO,CallSet -h $opt{refdir}/CallSet.header $sid\.int.vcf | bgzip > $sid.union.vcf.gz\n";
    if ($opt{refdir} =~ m/GRCh/) {
	print SH "bcftools annotate -a $opt{refdir}/ExAC.vcf.gz --columns CHROM,POS,AC_POPMAX,AN_POPMAX,POPMAX  $sid.union.vcf.gz | bgzip > $sid.exac.vcf.gz\n";
	print SH "java -Xmx10g -jar \$SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c \$SNPEFF_HOME/snpEff.config $opt{gver} $sid.exac.vcf.gz | java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar annotate $opt{refdir}/dbSnp.vcf.gz -  | java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar annotate $opt{refdir}/clinvar.vcf.gz - | java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar annotate $opt{refdir}/cosmic.vcf.gz - | java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar dbnsfp -v -db $opt{refdir}/dbNSFP.txt.gz - | java -Xmx10g -jar \$SNPEFF_HOME/SnpSift.jar gwasCat -db $opt{refdir}/gwas_catalog.tsv - |bgzip > $sid\.annot.vcf.gz\n";
    }else {
	print SH "java -Xmx10g -jar \$SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c \$SNPEFF_HOME/snpEff.config $opt{gver} $sid.union.vcf.gz |bgzip > $sid\.annot.vcf.gz\n";
    }
    print SH "bcftools stats $sid\.annot.vcf.gz > $sid\.stats.txt\n";
    print SH "plot-vcfstats -s -p $sid\_statplot $sid\.stats.txt\n";
    close SH;
}
