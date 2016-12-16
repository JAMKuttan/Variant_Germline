use strict;
use warnings;
use Statistics::R;

my ($vcfstat,$outprefix) = @ARGV;
my %parameterHash = ('width' => 800, 'height' => 600, 'pointsize' => 10, 'las' => 1);
$parameterHash{$_->[0]} = $_->[1] foreach(map {[split(/=/, $_, 2)]} @parameterList);
open IN, "<$vcfstat" or die $!;
my %typeCountHash;
while (my $line = <IN>) {
    if ($line =~ m/^SN/) {
	$line =~ s/number of |://;
	my ($sn,$id,$type,$num) = split(/\t/,$line);
	$typeCount{$type} = $num;
    }elsif ($line =~ m/^TSTV/) {
	my ($tstv,$id,$allelect,$snps,$transitions,$transversions,$numindels,$repeatcons,$repeatincos,$na);
	$tstv{'Transition'} = $transitions;
	$tstv{'Transversion'} = $transversion;
    }
}
$pngFile = $outprefix."snptype.png";
my @typeList = ('SNPs', 'indels', 'MNPs', 'others');
{
    my $R = Statistics::R->new();
    $R->run('counts <- c()');
    foreach my $type (@typeList) {
	$R->run(sprintf('counts["%s"] <- %d', $type, $tstv{$type}));
    }
    $R->run(sprintf('png(filename = "%s", width = %d, height = %d, units = "px", pointsize = %d, bg = "white", type = "cairo")', $pngFile, @parameterHash{'width', 'height', 'pointsize'}));
    $R->run(sprintf('barplot(counts, border = NA, las = %d)', $parameterHash{'las'}));
    $R->run('dev.off()');
    $R->stop();
}
$pngFile = $outprefix."tstv.png";
my @typeList = ('Transition', 'Transversion');
{
    my $R = Statistics::R->new();
    $R->run('counts <- c()');
    foreach my $type (@typeList) {
	$R->run(sprintf('counts["%s"] <- %d', $type, $typeCountHash{$type}));
    }
    $R->run(sprintf('png(filename = "%s", width = %d, height = %d, units = "px", pointsize = %d, bg = "white", type = "cairo")', $pngFile, @parameterHash{'width', 'height', 'pointsize'}));
    $R->run(sprintf('barplot(counts, border = NA, las = %d)', $parameterHash{'las'}));
    $R->run(sprintf('title(main = "TS/TV = %.3f")', $typeCountHash{'Transition'} / $typeCountHash{'Transversion'}));
    $R->run('dev.off()');
    $R->stop();
}
