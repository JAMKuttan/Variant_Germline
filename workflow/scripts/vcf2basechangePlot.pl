use strict;
use warnings;
use Statistics::R;

my ($vcfFile, $pngFile, @parameterList) = @ARGV;
my %parameterHash = ('width' => 800, 'height' => 600, 'pointsize' => 10, 'las' => 1);
$parameterHash{$_->[0]} = $_->[1] foreach(map {[split(/=/, $_, 2)]} @parameterList);
my %basechangeCountHash = ();
{
	my @columnList = ();
	open(my $reader, $vcfFile);
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^##/);
		next if($line =~ s/^#// && (@columnList = split(/\t/, $line)));
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		if($tokenHash{'FILTER'} eq 'PASS') {
			my $refBase = $tokenHash{'REF'};
			foreach my $altBase (split(/,/, $tokenHash{'ALT'})) {
				$basechangeCountHash{"$refBase>$altBase"} += 1 if(length($refBase) == 1 && length($altBase) == 1);
			}
		}
	}
	close($reader);
}
my @baseList = ('A', 'G', 'C', 'T');
{
	my $R = Statistics::R->new();
	$R->run('counts <- c()');
	foreach my $refBase (@baseList) {
		foreach my $altBase (@baseList) {
			if($refBase ne $altBase) {
				$R->run(sprintf('counts["%s"] <- %d', "$refBase>$altBase", $basechangeCountHash{"$refBase>$altBase"}));
			}
		}
	}
	$R->run(sprintf('png(filename = "%s", width = %d, height = %d, units = "px", pointsize = %d, bg = "white", type = "cairo")', $pngFile, @parameterHash{'width', 'height', 'pointsize'}));
	$R->run(sprintf('barplot(counts, border = NA, las = %d)', $parameterHash{'las'}));
	$R->run('dev.off()');
	$R->stop();
}
