#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
use Data::Printer;

#usage information
my $usage = <<USAGE;
grep_file_by_id V1.2, written by corephi
This program is equal to linux bash "grep -f", the differece
is: this program use hash to inprove speed.
----------------------------------------------------------
Usage: grep_file_by_id.pl -a afile -b bfile -o outfile
 Options:
  -a    afile name, the first line must be id, default STDIN
  -b    idlist file name, it must be one line per id
  -c    calculated tad number, default 5
  -s    shaking tad number
  -h 	help


istead.

USAGE
my $afile        = '';
my $bfile        = '';
my $cal_num      = 4;
my $shaking_num          = 4;
my $help_flag  = 0;
my $outfile      = '-';
die $usage
  unless GetOptions(
    "a:s" => \$afile,
    "b:s" => \$bfile,
    "c:i" => \$cal_num,
    "s:i" => \$shaking_num,
    "o:s" => \$outfile,
    "h"   => \$help_flag,
  );
die $usage if $help_flag;
die "a and b can not get from STDIN at the same time\n"
  if $afile && $bfile && $afile eq '-' && $bfile eq '-';
die "a and b can not get from STDIN at the same time\n"
  unless $afile || $bfile;

my $bfile_fh;
if ( $bfile && $bfile ne '-' ) {
    die "b file does not exists\n" unless -e $bfile;
    open $bfile_fh, "<", $bfile or die "cannot open file $bfile:$!\n";
}
else {
    $bfile_fh = *STDIN;
}

my $outfile_fh;
if ( $outfile && $outfile ne '-' ) {
    open $outfile_fh, ">", $outfile or die "cannot open file $outfile:$!\n";
}
else {
    $outfile_fh = *STDOUT;
}

#append content to a file
my $afile_fh;
if ( $afile && $afile ne '-' ) {
    die "a file does not exists\n" unless -e $afile;
    open $afile_fh, "<", $afile or die "cannot open file $afile:$!\n";
}
else {
    $afile_fh = *STDIN;
}
if ($help_flag) {
    my $header = readline $afile_fh;
    print $outfile_fh $header if $help_flag;
}


#store a file to hash
my $last_chr = 0;
my $index = 0;
my %chr_index = ();
my %shaking = ();
while (<$afile_fh>) {
    my $line = $_;
    next if /^#/;
    s/\r\n/\n/;
    my ($chr, $start, $end, $name) = split /\t/;
	my $locus = "$chr:$start-$end";
	
	$index = 1 if $chr ne $last_chr;
	if (exists $shaking{$locus}) {
		next; 
	}else{
	}	
	$chr_index{$chr}{$index} = $locus;
	$shaking{$locus} = 0;
	# p $index;
	$last_chr = $chr;
	$index++;
}

close $afile_fh;
#store bfile in to hash
while (<$bfile_fh>) {
    next if /^#/;
    s/\r\n/\n/;
    chomp;
    my ($chr, $start, $end, $name) = split /\t/;
	my $locus = "$chr:$start-$end";
	$shaking{$locus} = 1;	
}
close $bfile_fh;
# p %shaking;

my %calcued_shaking = ();
foreach my $chr (sort keys %chr_index) {
	my @bins = sort {$a <=> $b} keys %{$chr_index{$chr}};
	next if @bins < $cal_num;
	foreach my $index (1 .. @bins - $cal_num + 1) {
		my @slide_bins = ($index .. $index + $cal_num - 1);
		if (is_cluster($chr, \@slide_bins)) {
			foreach my $index (@slide_bins) {
				$calcued_shaking{$chr}{$index} = 1;
			}
		}
	}
}



foreach my $chr (sort keys %calcued_shaking) {
	foreach my $index (sort {$a <=> $b} keys %{$calcued_shaking{$chr}}) {
		my ($chr, $start, $end) = split /:|-/, $chr_index{$chr}{$index};
		print "$chr\t$start\t$end\tClutster_$.\n";
	}
}

sub is_cluster {
	my ($chr, $bin_rf) = @_;
	my @bins = @$bin_rf;
	my $is_cluster = 1;
	my @shaking_bins = ();
	foreach my $bin (@bins) {
		my $locus = $chr_index{$chr}{$bin};
		# p $locus;
		push @shaking_bins, $bin if $shaking{$locus};
	}
	# p @shaking_bins;
	$is_cluster = 0 if @shaking_bins < $shaking_num;
	$is_cluster = 0 unless $shaking{$chr_index{$chr}{$bins[0]}};
	$is_cluster = 0 unless $shaking{$chr_index{$chr}{$bins[$#bins]}};
	# p $is_cluster;
	return $is_cluster;
}

