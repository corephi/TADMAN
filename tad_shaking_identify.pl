#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
# use Data::Printer;

#usage information
my $usage = <<USAGE;
seperate_bed.pl V1.0, written by corephi
This program is similar to linux bash "cut -f", the differece
is this can not only output by ordered col numbers, but also 
can output by ordered col names
----------------------------------------------------------
Usage: seperate_bed.pl -a cond1.tad -b cond2.tad -o outfile
 Options:
  -a    input file name, default STDIN
  -b    input file name, default STDIN
  -o    output file name, default STDOUT
  -b    bin size default, 40000 
  -r    rename the new track name;
Note: This scrpit can also be used in pipeline.
cat in.txt | seperate_bed.pl -c 3,2,1 
     

USAGE
my $a_file    = '';
my $b_file    = '';
my $bin_number = 3;
my $bin_size = 40000;
my $outfile    = '-';
my $rename = 0;
die $usage
  unless GetOptions(
    "a:s" => \$a_file,
    "b:s" => \$b_file,
    "r" => \$rename,
    "o:s" => \$outfile,
  );

die "a and b can not get from STDIN at the same time\n"
  if $a_file && $b_file eq '-';
  
my $out_fh;
if ( $outfile && $outfile ne '-' ) {
    open $out_fh, ">", $outfile or die "cannot open file $outfile:$!\n";
}
else {
    $out_fh = *STDOUT;
}

my $a_fh;
if ( $a_file && $a_file ne '-' ) {
    die "a file does not exists\n" unless -e $a_file;
    open $a_fh, "<", $a_file or die "cannot open file $a_file:$!\n";
}
else {
    $a_fh = *STDIN;
}

my $b_fh;
if ( $b_file && $b_file ne '-' ) {
    die "a file does not exists\n" unless -e $b_file;
    open $b_fh, "<", $b_file or die "cannot open file $b_file:$!\n";
}
else {
    $b_fh = *STDIN;
}

#####################################################################################################
#Reading bin size
#####################################################################################################
my %bins = ();
while (<$a_fh>) {
	s/\r?\n//;
	my @lines = split /\s+/;
	my ($chr, $start, $end, $name) = @lines;
	my $locus = "$chr:$start-$end";
	$bins{$locus}{CondA} = $name;
}
close $a_fh;

while (<$b_fh>) {
	s/\r?\n//;
	my @lines = split /\s+/;
	my ($chr, $start, $end, $name) = @lines;
	my $locus = "$chr:$start-$end";
	$bins{$locus}{CondB} = $name;
}
close $b_fh;

#####################################################################################################
#Findding TAD shaking
#####################################################################################################
foreach my $bin (keys %bins) {
	my ($chr, $start, $end) = split /:|-/, $bin;
	my $start_l = $start - $bin_size;
	my $end_l = $end - $bin_size;
	my $bin_l = "$chr:$start_l-$end_l";

	my $start_r = $start + $bin_size;
	my $end_r = $end + $bin_size;
	my $bin_r = "$chr:$start_l-$end_l";
	
	my $bin =  "$chr:$start-$end";
	my $type = check_around_bind($bin);
	print $out_fh "$chr\t$start\t$end\t$type\n";
	
}

sub check_around_bind {
	my $bin = shift;
	my ($chr, $start, $end) = split /:|-/, $bin;
	my $start_l = $start - $bin_size;
	my $end_l = $start;
	my $bin_l = "$chr:$start_l-$end_l";
	
	my $start_r = $end ;
	my $end_r = $end + $bin_size;
	my $bin_r = "$chr:$start_r-$end_r";
	# warn "$bin\t$bin_l\t$bin_r\n";
	my $type = "others";
	my $a_tid = exists $bins{$bin}{CondA} ? $bins{$bin}{CondA} : '';
	my $b_tid = exists $bins{$bin}{CondB} ? $bins{$bin}{CondB} : '';
	
	
	if (exists $bins{$bin_l} && exists $bins{$bin_r}) {
		my $l_a_tid = exists $bins{$bin_l}{CondA} ? $bins{$bin_l}{CondA} : '';
		my $l_b_tid = exists $bins{$bin_l}{CondB} ? $bins{$bin_l}{CondB} : '';
		
		my $r_a_tid = exists $bins{$bin_r}{CondA} ? $bins{$bin_r}{CondA} : '';
		my $r_b_tid = exists $bins{$bin_r}{CondB} ? $bins{$bin_r}{CondB} : '';
		
		my @condA_tad_ids = grep {$_} ($a_tid, $l_a_tid, $r_a_tid );

		my %condA_tad = map {$_ => 1 } @condA_tad_ids;
		my $condA_tad_num = keys %condA_tad;
		
		my @condB_tad_ids = grep {$_} ($b_tid, $l_b_tid, $r_b_tid );
		my %condB_tad = map {$_ => 1 } @condB_tad_ids;
		my $condB_tad_num = keys %condB_tad;
		
		if (@condA_tad_ids == 3 && @condB_tad_ids == 3) {
			# p %condA_tad;
			# p %condB_tad;
			
			if ($condA_tad_num == 1 && $condB_tad_num == 1) {
				$type = "NonChange";
			}elsif (($condA_tad_num == 1 &&  $condB_tad_num == 2 ) || ($condA_tad_num == 2 &&  $condB_tad_num == 1 ) ) {
				$type = "split";
			}elsif( $condA_tad_num == 2 &&  $condB_tad_num == 2 ) {
				if ( ($l_a_tid eq $a_tid  && $a_tid ne $r_a_tid) && ($l_b_tid ne $b_tid && $b_tid eq $r_b_tid) ||
					($l_a_tid ne $a_tid  && $a_tid eq $r_a_tid) && ($l_b_tid eq $b_tid && $b_tid ne $r_b_tid) ) {
						$type = 'shaking';
					}
				else{
						$type = 'NonChange';
				}
				
			}else{
			}			
		}else{
			$type = "others";
		}
		
	}else{
		$type = "others";
	}
	return $type;
}