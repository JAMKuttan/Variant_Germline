#!/usr/bin/perl -w
#combine_variants.pl

my $prefix = shift @ARGV;

my @prioritize = ('sam','ssvar','platypus','gatk','hotspot');

my @prev = ('nocaller');

open OUT, ">combineVariants.sh" or die $!;

foreach $caller (@prioritize) {
  my $rm = join('\\|', @prev);
  print OUT "grep $caller $prefix\_integrate.bed | grep -v '$rm' > $caller.bed\n";
  print OUT "bedtools intersect -header -a $prefix\.$caller\.vcf.gz -b $caller.bed |bgzip > $caller\.part.vcf.gz\n";
  print OUT "tabix $caller\.part.vcf.gz\n";
  
  push @prev, $caller;
}
print OUT "bcftools concat -a -O v -o $prefix\.int.vcf *.part.vcf.gz\n";
close OUT;
system("sh combineVariants.sh");
