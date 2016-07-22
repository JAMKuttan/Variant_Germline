#!/usr/bin/perl -w
#uploadqc.pl

open OUT, ">sequence.stats.txt" or die $!;
print OUT join("\t",'Sample','total','pairs','maprate','propair','percdups',
	       'medinsert','avginsert','stdinsert','perc20x','perc50x',
	       'perc100x','perc200x','perc500x'),"\n";
my @statfiles = @ARGV;

foreach $sfile (@statfiles) {
    $sfile =~ m/(\S+)\.flagstat.txt/;
    my $prefix = $1;
    my %hash;
    open FLAG, "<$prefix\.flagstat.txt" or die $!;
    my ($total, $read1ct,$read2ct,$maprate,$concorrate);
    while (my $line = <FLAG>) {
	chomp($line);
	if ($line =~ m/(\d+) \+ \d+ in total/) {
	    $hash{total} = $1;
	}elsif ($line =~ m/(\d+) \+ \d+ read1/) {
	    $hash{pairs} = $1;
	}elsif ($line =~ m/(\d+) \+ \d+ mapped/) {
	    $hash{maprate} = 100*sprintf("%.4f",$1/$hash{total});
	}elsif ($line =~ m/(\d+) \+ \d+ properly paired/) {
	    $hash{propair} = 100*sprintf("%.4f",$1/$hash{total});
	}
    }
    open DUP, "<$prefix\.libcomplex.txt" or die $!;
    while (my $line = <DUP>) {
	chomp($line);
	if ($line =~ m/## METRICS/) {
	    $header = <DUP>;
	    $nums = <DUP>;
	    chomp($header);
	    chomp($nums);
	    my @stats = split(/\t/,$header);
	    my @row = split(/\t/,$nums);
	    my %info;
	    foreach my $i (0..$#stats) {
		$info{$stats[$i]} = $row[$i];
	    }
	    $hash{percdups} = sprintf("%.4f",$info{PERCENT_DUPLICATION});
	}
    }
    open DUP, "<$prefix\.hist.txt" or die $!;
    while (my $line = <DUP>) {
	chomp($line);
	if ($line =~ m/## METRICS/) {
	    $header = <DUP>;
	    $nums = <DUP>;
	    chomp($header);
	    chomp($nums);
	    my @stats = split(/\t/,$header);
	    my @row = split(/\t/,$nums);
	    my %info;
	    foreach my $i (0..$#stats) {
		$info{$stats[$i]} = $row[$i];
	    }
	    $hash{medinsert} = sprintf("%.0f",$info{MEDIAN_INSERT_SIZE});
	    $hash{avginsert} = sprintf("%.0f",$info{MEAN_INSERT_SIZE});
	    $hash{stdinsert} = sprintf("%.0f",$info{STANDARD_DEVIATION});
	}
    }
    my %cov;
    open COV, "<$prefix\.genomecov.txt" or die $!;
    while (my $line = <COV>) {
      chomp($line);
      my ($all,$depth,$bp,$total,$percent) = split(/\t/,$line);
      $cov{$depth} = $percent;
    }
    my @depths = sort {$a <=> $b} keys %cov;
    my @perc = @cov{@depths};
    my @cum_sum = cumsum(@perc);
    $hash{'perc20x'} = 100*sprintf("%.4f",1-$cum_sum[20]);
    $hash{'perc50x'} = 100*sprintf("%.4f",1-$cum_sum[50]);
    $hash{'perc100x'} = 100*sprintf("%.4f",1-$cum_sum[100]);
    $hash{'perc200x'} = 100*sprintf("%.4f",1-$cum_sum[200]);
    $hash{'perc500x'} = 100*sprintf("%.4f",1-$cum_sum[500]);
    print OUT join("\t",$prefix, $hash{total},$hash{pairs},$hash{maprate},
		   $hash{propair},$hash{percdups},$hash{medinsert},
		   $hash{avginsert},$hash{stdinsert},$hash{'perc20x'},$hash{'perc50x'},
		   $hash{'perc100x'},$hash{'perc200x'},$hash{'perc500x'}),"\n";
  }

sub cumsum {
  my @nums = @_;
  my @cumsum = ();
  my $mid = 0;
  for my $i (0..$#nums) {
    $mid += $nums[$i];
    push(@cumsum,$mid);
  }
  return @cumsum;
}
