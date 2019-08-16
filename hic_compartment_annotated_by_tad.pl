#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use Cwd qw[abs_path];
use IPC::Cmd qw[can_run run];
use threads;
use Thread::Queue;
# use Data::Printer;
use File::Temp qw/ :POSIX /;


my $usage = <<USAGE;
SYSNOPSIS
hic_compartment_annotated_by_tad.pl [options] A.bed B.bed C.bed

 Options:
   -i --in-bed        bed formated compartment file, 6 column
   -c --column        signal column used for retrive
   -f --folder        use folder to infer basename
   -o --out-file      output file
   -m --mode          'TAD' or "BIN" mode, default "TAD"
   -h --help          print this help.
USAGE

# my $in_compartment = '-';
my $in_compartment = '../Compartment/HiC.ESC.bed';
my $out_bed = '-';
my $mode = 'BIN';
my $tad_bed = '../../TADs/BEDs/ESC.TAD.bed';
my $help        = 0;
die $usage
  unless GetOptions(
    "i|in-bed=s"       => \$in_compartment,
    "o|out-bed:s" => \$out_bed,
    "t|tad-bed:s" => \$tad_bed,
    "m|mode:s" => \$mode,
    "help"          => \$help,
  );
die $usage if $help;


#check and open TAD bed
die "TAD file does not exits:$!" unless -e $tad_bed;
open my $tad_bed_fh, "<", $tad_bed or die "cannot open file $tad_bed:$!\n";


#open out put file  
my $out_bed_fh;
if ( $out_bed && $out_bed ne '-' ) {
    open $out_bed_fh, ">", $out_bed or die "cannot open file $out_bed:$!\n";
}
else {
    $out_bed_fh = *STDOUT;
}

#open compartment bed file
my $in_compartment_fh;
if ( $in_compartment && $in_compartment ne '-' ) {
    die "a file does not exists\n" unless -e $in_compartment;
    open $in_compartment_fh, "<", $in_compartment or die "cannot open file $in_compartment:$!\n";
}
else {
    $in_compartment_fh = *STDIN;
}

#open tad bed file
my $in_tad_fh;
if ( $tad_bed && $tad_bed ne '-' ) {
    die "a file does not exists\n" unless -e $tad_bed;
    open $in_tad_fh, "<", $tad_bed or die "cannot open file $tad_bed:$!\n";
}
else {
    $in_tad_fh = *STDIN;
}

my $bed_tmp_file =  tmpnam();
system "cut -f 1-6 $in_compartment > $bed_tmp_file ";

my $commmand = "bedtools intersect -a $bed_tmp_file -b $tad_bed -wo";
my @results = split /\r?\n/, `$commmand`;


####################################################################################################################
#store the Maping info
####################################################################################################################
my %tads = ();
my %compartments = ();
my %tads2compartments = ();
my %compartments2tads = ();


foreach (@results) {
	s/\r?\n//;
	my ($chr_c, $start_c, $end_c, $name_c, $score_c, $strand_c, $chr_t, $start_t, $end_t, $name_t ) = split /\t/;
	$tads2compartments{$name_t}{$name_c} = $score_c;
	$compartments2tads{$name_c} = $name_t;  #为避免冲突，限定：一个Compartments只能属于1个TAD
	
	$tads{$name_t}{chr} = $chr_t;
	$tads{$name_t}{start} = $start_t;
	$tads{$name_t}{end} = $end_t;

	$compartments{$name_c}{chr} = $chr_c;
	$compartments{$name_c}{start} = $start_c;
	$compartments{$name_c}{end} = $end_c;
	$compartments{$name_c}{score} = $end_c;
	$compartments{$name_c}{strand} = $end_c;
}

####################################################################################################################
# calculate compartment in tad format
####################################################################################################################

foreach my $name_t (keys %tads) {
	my @scores = ();
	foreach my $name_c (keys %{$tads2compartments{$name_t}}) {
		push @scores, $tads2compartments{$name_t}{$name_c};
	}
	
	my $sum = 0;
	foreach my $score (@scores) {
		$sum += $score;
	}
	
	my $avg = $sum / @scores;
	
	my $compartment = $avg > 0 ? "A" : "B";
	
	$tads{$name_t}{compartment} = $compartment;
	$tads{$name_t}{score} = $avg;
}

#rewrite to compartments
foreach my $name_c (keys %compartments) {
	my $name_t = $compartments2tads{$name_c};
	my $compartment = $tads{$name_t}{compartment};
	$compartments{$name_c}{compartment} = $compartment;
}




############################################################
#Writing files
############################################################
if ($mode eq 'BIN') {
	while (<$in_compartment_fh>) {
		next if /^#/;
		s/\r?\n//;
		my ($chr_c, $start_c, $end_c, $name_c, $score_c, $strand_c) = split /\t/;
		my $compartment_c  = '.';
		if (exists  $compartments{$name_c}) {
			$compartment_c = $compartments{$name_c}{compartment}
		}else{
			warn "$name_c has gap\n";
			$compartment_c = $score_c > 0 ? "A" : "B";	
		}
		print $out_bed_fh "$chr_c\t$start_c\t$end_c\t$name_c\t$score_c\t$strand_c\t$compartment_c\n";
	}	
}else{
	while (<$in_tad_fh>) {
		next if /^#/;
		s/\r?\n//;
		my ( $chr_t, $start_t, $end_t, $name_t ) = split /\t/;
		my $compartment_t = $tads{$name_t}{compartment};
		my $score_t = $tads{$name_t}{score};
		print $out_bed_fh "$chr_t\t$start_t\t$end_t\t$name_t\t$score_t\t.\t$compartment_t\n";
	}	
}

