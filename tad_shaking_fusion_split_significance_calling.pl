#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
# use Data::Printer;

#usage information
my $usage = <<USAGE;
tad_shaking_fusion_split_significance_calling.pl V1.0, written by corephi
This program is similar to linux bash "cut -f", the differece
is this can not only output by ordered col numbers, but also 
can output by ordered col names
----------------------------------------------------------
Usage: tad_shaking_fusion_split_significance_calling.pl -a AvsB.bed -b AvsB.tad -o outfile
 Options:
  -d    differencial boundary significance bed formated file name, 
            the 4th column is P or Q value
  -s    tad shaking bin bed file, default STDIN
  -o    output file name, default STDOUT
  -b    bin size default, 40000 
Note: This scrpit can also be used in pipeline.
     

USAGE
my $diff_file    = '-';
my $shake_file    = '-';
my $bin_size = 40000;
my $outfile    = '-';
die $usage
  unless GetOptions(
    "d:s" => \$diff_file,
    "s:s" => \$shake_file,
    "o:s" => \$outfile,
  );

die "d and s can not get from STDIN at the same time\n"
  if $diff_file eq '-' && $shake_file eq '-';
  
my $out_fh;
if ( $outfile && $outfile ne '-' ) {
    open $out_fh, ">", $outfile or die "cannot open file $outfile:$!\n";
}
else {
    $out_fh = *STDOUT;
}

my $diff_fh;
if ( $diff_file && $diff_file ne '-' ) {
    die "a file does not exists\n" unless -e $diff_file;
    open $diff_fh, "<", $diff_file or die "cannot open file $diff_file:$!\n";
}
else {
    $diff_fh = *STDIN;
}

my $shake_fh;
if ( $shake_file && $shake_file ne '-' ) {
    die "a file does not exists\n" unless -e $shake_file;
    open $shake_fh, "<", $shake_file or die "cannot open file $shake_file:$!\n";
}
else {
    $shake_fh = *STDIN;
}

#####################################################################################################
#Reading bin size
#####################################################################################################
my %diff = ();
while (<$diff_fh>) {
	s/\r?\n//;
	my @lines = split /\s+/;
	my ($chr, $start, $end, $score) = @lines;
	my $locus = "$chr:$start-$end";
	$diff{$locus} = $score;
}
close $diff_fh;

while (<$shake_fh>) {
	s/\r?\n//;
	my @lines = split /\s+/;
	my ($chr, $start, $end, $name) = @lines;
	my $locus = "$chr:$start-$end";
	my $pvalue = get_p_value($locus);
	print $out_fh "$chr\t$start\t$end\t$pvalue\n";
}
close $shake_fh;
close $out_fh;

sub get_p_value {
	my $bin = shift;
	my ($chr, $start, $end) = split /:|-/, $bin;
	my $start_l = $start - $bin_size;
	my $end_l = $start;
	my $bin_l = "$chr:$start_l-$end_l";
	
	my $start_r = $end ;
	my $end_r = $end + $bin_size;
	my $bin_r = "$chr:$start_r-$end_r";
	
	my @scores = ();
	push @scores, $diff{$bin} if exists $diff{$bin} ;
	push @scores, $diff{$bin_l} if exists $diff{$bin_l};
	push @scores, $diff{$bin_r} if exists $diff{$bin_r};
	
	@scores = sort {$a <=> $b} @scores;
	return @scores ? $scores[0] : 1;
	
}