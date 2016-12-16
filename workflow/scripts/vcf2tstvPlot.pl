use strict;
use warnings;
use Statistics::R;

my ($vcfFile, $pngFile, @parameterList) = @ARGV;
my %parameterHash = ('width' => 800, 'height' => 600, 'pointsize' => 10, 'las' => 1);
$parameterHash{$_->[0]} = $_->[1] foreach(map {[split(/=/, $_, 2)]} @parameterList);
my %basechangeTypeHash = ();
$basechangeTypeHash{$_} = 'Transition' foreach('A>G', 'G>A', 'C>T', 'T>C');
$basechangeTypeHash{$_} = 'Transversion' foreach('A>T', 'T>A', 'C>G', 'G>C');
$basechangeTypeHash{$_} = 'Transversion' foreach('A>C', 'C>A', 'T>G', 'G>T');
my %typeCountHash = ();
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
				if(length($refBase) == 1 && length($altBase) == 1) {
					$typeCountHash{$_} += 1 if(defined($_ = $basechangeTypeHash{"$refBase>$altBase"}));
				}
			}
		}
	}
	close($reader);
}
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
