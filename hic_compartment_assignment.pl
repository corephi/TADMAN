#!/usr/bin/perl -w
use strict;
use 5.010;
use File::Basename;
use IPC::Cmd qw[can_run run];
use File::Copy;
use File::Temp qw[tempdir];
use Getopt::Long;
use threads;
use Thread::Queue;
# use Data::Printer;

my $usage = <<USAGE;

SYSNOPSIS
hic_compartment_assignment.pl V1.0, written by corephi

This program is used to compartmentate chromatin.

----------------------------------------------------------
More scripts? Join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, 
please mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------

hic_compartment_assignment.pl [options] file|glob

 Options:	
   -g|--genome-size      genome size file
   -e|--eigen-bed        hic compartment eigen value, bed6 format
   -s|--signal-bed       bed format sigal, such as H3K9me3, ATAC-Seq
   -d|--signal-col       signal column used in bed. 1-based, default 5
   -f|--eigen-col        eigen column used in bed. 1-based, default 4
   -l|--lower-A          whether signal lower is compartment A, such as H3K9me3,
                           default, off, ie ATAC-Seq;
   -w|--window           window-size used for compartmentation
                           default 10000
   -o|--out-file         output file
   -h|--help             print this usage.
   

USAGE
my $outfile             = '-';
my $genome_size_file     = '';
my $hic_eigen_file       = '-';
my $singal_bed_file      = '';
my $window_size          = 100000;
my $signal_col           = 4;
my $eigen_col           = 4;
my $normalization_method = "KR";
my $help                 = 0;
my $lower_A              = 0;
die $usage
  unless GetOptions(
    "g|genome-size:s" => \$genome_size_file,
    "e|eigen-bed:s"   => \$hic_eigen_file,
    "f|eigen-col=i"  => \$eigen_col,
    "d|signal_col-col=i"  => \$signal_col,
    "s|signal-bed:s"  => \$singal_bed_file,
    "o|out-file:s"      => \$outfile,
    "w|window=i"      => \$window_size,
    "h|help"          => \$help,
    "l|lower-A"       => \$lower_A,
  );
##########################################################################################
#Checking the parameter infomation
##########################################################################################
die $usage if $help;

#check the running environment
can_run('bedtools')    or die 'bedtools is not installed!';

#check the input files
die "genome size file:$genome_size_file does not exists!\n"
  unless -e $genome_size_file;
  
  
die "signal bed file:$singal_bed_file does not exists!\n"
  unless -e $singal_bed_file;

#prepare the tmp folder
my $tmp_folder = tempdir( CLEANUP => 1 );
  
  
##########################################################################################
#make the window
##########################################################################################
warn "Make genome windows...\n";
my $window_file = $tmp_folder . "/Genome.${window_size}.bin.bed";
my $window_command =
"bedtools makewindows -g $genome_size_file -w $window_size >$window_file";
system($window_command);

##########################################################################################
#Calculate the coverage by bin
##########################################################################################

#sort the original signal bed files
warn "Sorting Original singal bed files\n";
my $signal_file_basename =
  basename( $singal_bed_file, ".bed", ".bdg", ".narrowPeak", ".broadPeak" );
my $sorted_signal_bed_file = "$tmp_folder/$signal_file_basename.sorted.bed";
my $sort_command =
"bedtools sort -i $singal_bed_file -faidx $genome_size_file > $sorted_signal_bed_file";
system($sort_command);

#calculate the coverage
warn "Caluclaing the coverage by each windows...\n";
my $coverage_bed_file = "$tmp_folder/SignalCoverage.bed";
my $map_command =
"bedtools map -g $genome_size_file -a $window_file -b $sorted_signal_bed_file -c $signal_col -o mean >$coverage_bed_file";
system($map_command);

##########################################################################################
#Read the SignalCoverage
##########################################################################################
warn("Reading calculate the SignalCoverage...\n");

#to keep the order
open my $genome_size_fh, "<", $genome_size_file
  or die "Cannot open genome size file\n";
my %genome  = ();
my @seq_ids = ();    #to keep the orders as genome.size
while (<$genome_size_fh>) {
    s/\r?\n//;
    my ( $seq_id, $len ) = split /\t/;
    next if $seq_id =~ /chrM|chrMT|MT/igx;
    $genome{$seq_id}{len} = $len;
    push @seq_ids, $seq_id;
}


#read the coverage
open my $coverage_fh, "<", $coverage_bed_file
  or die "cannot open window file:$coverage_bed_file\n";
my %data        = ();
my $wid         = 0;
while (<$coverage_fh>) {
    next if /^#/;
    s/\r?\n//;
    my ( $seq_id, $start, $end, $score ) = split /\t/;
    next if $seq_id =~ /chrM|chrMT|MT/igx;
    my $wid    = "$seq_id:$start-$end";
    next if $seq_id =~ /chrM|chrMT|MT/igx;
    $data{$seq_id}{$wid}{name}   = $wid;
    $data{$seq_id}{$wid}{seq_id} = $seq_id;
    $data{$seq_id}{$wid}{start}  = $start;    #0-based
    $data{$seq_id}{$wid}{end}    = $end;

    $data{$seq_id}{$wid}{signal} = $score;
}
close $coverage_fh;


##########################################################################################
#Read the EigenVector
##########################################################################################
my $eigen_fh;
warn("Reading the EigenVector...\n");

if ( $hic_eigen_file && $hic_eigen_file ne '-' ) {
    die "HiC compartment eigen file:$hic_eigen_file  does not exists\n" unless -e $hic_eigen_file;
    open $eigen_fh, "<", $hic_eigen_file or die "cannot open file $hic_eigen_file:$!\n";
}
else {
    $eigen_fh = *STDIN;
}

while (<$eigen_fh>) {
    next if /^#/;
    s/\r?\n//;
    my ( $seq_id, $start, $end ) = my @lines = split /\t/;
    next if $seq_id =~ /chrM|chrMT|MT/igx;
    my $wid    = "$seq_id:$start-$end";
	my $eigen_value = $lines[$eigen_col - 1];

	if ( $eigen_value eq "NaN" ) {
		$data{$seq_id}{$wid}{eigen}       = '.';
		$data{$seq_id}{$wid}{compartment} = '.';
	}
	else {
		$data{$seq_id}{$wid}{eigen} = $eigen_value;
		$data{$seq_id}{$wid}{compartment} = $eigen_value > 0 ? "1" : "2";
	}

}
close $eigen_fh;

##########################################################################################
#Asign the compartment
##########################################################################################

foreach my $seq_id (@seq_ids) {
    next if $seq_id =~ /chrM|chrMT|MT/igx;
    my @compartment1 = ();
    my @compartment2 = ();

    foreach my $wid ( keys %{ $data{$seq_id} } ) {
        my $t_compartment = $data{$seq_id}{$wid}{compartment};
        die "$seq_id\t$wid\n" unless $t_compartment;
        my $t_score = $data{$seq_id}{$wid}{signal};
        if ( $t_compartment eq '1' ) {
            push @compartment1, $t_score unless $t_score eq '.';
        }
        else {
            push @compartment2, $t_score unless $t_score eq '.';
        }

    }
    my $avg_compart1 = avg(@compartment1);
    my $avg_compart2 = avg(@compartment2);

    if ( $avg_compart1 > $avg_compart2 ) {
        foreach my $wid ( keys %{ $data{$seq_id} } ) {
            my $t_compartment = $data{$seq_id}{$wid}{compartment};
            next if $t_compartment eq '.';
            if ($lower_A) {
                $data{$seq_id}{$wid}{compartment} =
                  $t_compartment eq '1' ? "B" : "A";
            }
            else {
                $data{$seq_id}{$wid}{compartment} =
                  $t_compartment eq '1' ? "A" : "B";
            }
        }
    }
    else {
        foreach my $wid ( keys %{ $data{$seq_id} } ) {
            my $t_compartment = $data{$seq_id}{$wid}{compartment};
            next if $t_compartment eq '.';
            if ($lower_A) {
                $data{$seq_id}{$wid}{compartment} =
                  $t_compartment eq '1' ? "A" : "B";
            }
            else {
                $data{$seq_id}{$wid}{compartment} =
                  $t_compartment eq '1' ? "B" : "A";
            }
        }
    }
}


##########################################################################################
#out put the tmp results
##########################################################################################
my $assiged_bed_tmp = $tmp_folder."Assigned.Compartment.bed";
open my $assiged_bed_tmp_fh, ">", $assiged_bed_tmp
	or die "Cannot create file:$assiged_bed_tmp $!\n";
foreach my $seq_id (@seq_ids) {
    foreach my $wid ( sort keys %{ $data{$seq_id} } ) {
        my $name        = $data{$seq_id}{$wid}{name};
        my $seq_id      = $data{$seq_id}{$wid}{seq_id};
        my $start       = $data{$seq_id}{$wid}{start};
        my $end         = $data{$seq_id}{$wid}{end};
        my $singal      = $data{$seq_id}{$wid}{signal};
        my $eigen       = $data{$seq_id}{$wid}{eigen};
        my $compartment = $data{$seq_id}{$wid}{compartment};

		# change the EigenVector to positive to Compartment A, and negative to Compartment B
        if ( $compartment eq "A" && $eigen < 0 ) {
            $eigen = $eigen * -1;
        }
        elsif ( $compartment eq "B" && $eigen > 0 ) {
            $eigen = $eigen * -1;
        }

        print $assiged_bed_tmp_fh
          "$seq_id\t$start\t$end\t$name\t$eigen\t.\t$singal\t$compartment\n";
    }

}
close $assiged_bed_tmp_fh;


my $assiged_sorted_bed_tmp = $tmp_folder."Assigned.Sorted.Compartment.bed";

my $assigned_sorted_command = "bedtools sort -i $assiged_bed_tmp -faidx $genome_size_file > $assiged_sorted_bed_tmp";
system($assigned_sorted_command);


my $outfile_fh;
if ( $outfile && $outfile ne '-' ) {
    open $outfile_fh, ">", $outfile or die "cannot open file $outfile:$!\n";
}
else {
    $outfile_fh = *STDOUT;
}
print $outfile_fh
  "#Chr\tStart\tEnd\tName\tEigenValue\tStrand\tSingal\tCompartment\n";

  
open my $assiged_sorted_bed_tmp_fh, "<", $assiged_sorted_bed_tmp
	or die "Cannot open file:$assiged_sorted_bed_tmp $!\n";
	
while(<$assiged_sorted_bed_tmp_fh>){
	print $outfile_fh $_;
}  
  
  

sub avg {
    my $sum = 0;
    foreach my $num (@_) {
        $sum += $num;
    }
    return $sum/@_;
    # return $sum;
}
