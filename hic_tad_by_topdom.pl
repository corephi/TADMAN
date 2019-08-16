#!/usr/bin/perl -w
use strict;
use 5.010;
use File::Basename;
use IPC::Cmd qw[can_run run];
use Cwd qw [abs_path];
use File::Copy;
use Getopt::Long;
use threads;
use Thread::Queue;
# use Data::Printer;

my $usage = <<USAGE;

SYSNOPSIS
hic_tads.pl V1.1, written by corephi

This program is used to compartmentate chromatin.

----------------------------------------------------------
More scripts? Join "bioinformatics*CN" QQ group: 276151571
If you have any probel or suggestions about this program, 
please mail to: hpyu\@genetics.ac.cn
----------------------------------------------------------

hic_tads.pl [options] file|glob

 Options:	
   -b|--bed           BED file with bins coordinates.
   -c|--hic           HiC-Pro 3 columns normalized matrix
   -w|--window        window-size used for TopDom, default 5,
                      it can also be set as 3,5,7,9,12,15
   -o|--output        output folder
   -p|--progress      progress, default 88
   -h|--help          print this usage.
   
USAGE
my $out_folder = dirname 'TADs';
my $progress   = 88;
my $bed_file   = '';
my $hic_file   = '';

# my $window = "3,5,7,9,12,15";
my $window = "5";
my $help   = 0;
die $usage
  unless GetOptions(
    "b|bed=s"      => \$bed_file,
    "c|hic=s"      => \$hic_file,
    "w|window:s"   => \$window,
    "o|output:s"   => \$out_folder,
    "p|progress=i" => \$progress,
    "h|help"       => \$help,
  );
##########################################################################################
#Checking the parameter infomation
##########################################################################################
die $usage if $help;

#check the running environment
can_run('TopDom.R')         or die 'TopDom.R is not installed!';
can_run('sparseToDense.py') or die 'sparseToDense.py is not installed!';

#check the input files
die "BIN bed file:$bed_file does not exists!\n" unless -e $bed_file;
$bed_file = abs_path($bed_file);
die "normalized hic contact file:$hic_file does not exists!\n"
  unless -e $hic_file;
$hic_file = abs_path($hic_file);

#window-size
my @windows = split /,/, $window;
die "No window size:@windows\n" unless @windows;

#mkdir
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;
my $log_folder = $out_folder . "/logs";
mkdir $log_folder unless -d $log_folder;

$out_folder = abs_path($out_folder);
$log_folder = abs_path($log_folder);

##########################################################################################
#Read the chromosome files
##########################################################################################
warn "Reading BIN BED file...\n";
open my $bin_fh, "<", $bed_file
  or die "Cannot open file:$!\n";
my %seq_ids = ();
while (<$bin_fh>) {
    s/\r?\n//;
    my ( $seq_id, ) = split /\t/;
    next if $seq_id =~ /chrM|chrMT|MT/igx;
    $seq_ids{$seq_id} = 1;
}
##########################################################################################
# sparseToDense
# Reformat the 3 column iced interaction files to N+3 matrix
##########################################################################################
warn "Converting HiC-Pro matrix to TopDom matrix...\n";
my $matrix_folder = abs_path( $out_folder . "/Matrix" );
mkdir $matrix_folder unless -d $matrix_folder;
chdir $matrix_folder;
my @tmp = split /\//, $out_folder;
my $prefix = $tmp[-1];
my $sparseToDense_command =
"sparseToDense.py -d -c -o $prefix -b $bed_file $hic_file 2>$log_folder/sparseToDense.log";
run_command($sparseToDense_command);
chdir $out_folder;

##########################################################################################
#TopDom
##########################################################################################
warn "Identify topDom...\n";
my $topDom_folder = $out_folder . "/TopDom";
mkdir $topDom_folder unless -d $topDom_folder;
$topDom_folder = abs_path($topDom_folder);
chdir $topDom_folder;
my @topdom_commands = ();

my %topdom_files = ();
foreach my $window_size (@windows) {
    foreach my $seq_id ( sort keys %seq_ids ) {
        next if $seq_id =~ /chrM|chrMT|MT/igx;
        my $matrix_file = "$matrix_folder/${seq_id}_$prefix";
        die "Converted matrix file does not exists\n" unless -e $matrix_file;
        my $stderr = "$log_folder/TopDom_$seq_id.$window_size.log";
        my $stdout = "$log_folder/TopDom_$seq_id.$window_size.txt";
        my $matrix_command =
"TopDom.R -i $matrix_file -o $seq_id.$window_size -w $window_size 1> $stdout 2>$stderr";
        $topdom_files{$seq_id}{matrix} = abs_path($matrix_file);
        $topdom_files{$seq_id}{pcc}{$seq_id}{$window_size} =
          abs_path("$seq_id.$window_size.pcc");
        $topdom_files{$seq_id}{bed}{$seq_id}{$window_size} =
          abs_path("$seq_id.$window_size.bed");
        $topdom_files{$seq_id}{domain}{$seq_id}{$window_size} =
          abs_path("$seq_id.$window_size.domain");
        push @topdom_commands, $matrix_command;
    }

}

run_parallel( $progress, @topdom_commands );
chdir $out_folder;

##########################################################################################
# Merging PCC
##########################################################################################
warn "Merging PCCs...\n";
my @pcc_files = ();
foreach my $window_size (@windows) {
    foreach my $seq_id ( sort keys %seq_ids ) {
        next if $seq_id =~ /chrM|chrMT|MT/igx;
        my $pcc_file = "$topdom_files{$seq_id}{pcc}{$seq_id}{$window_size} ";
        push @pcc_files, $pcc_file;
    }
}
my $pcc_file          = abs_path("$out_folder/PCC.txt");
my $merge_pcc_command = "cat @pcc_files | sort -r | uniq > $pcc_file";
run_command($merge_pcc_command);

##########################################################################################
# Choose the window size
##########################################################################################
warn "Automatically choose TADs...\n";
my $box_plot_command = "PCC_Boxplot.R -i $pcc_file -o PCC";
run_command($box_plot_command);
my $pcc_median_file = abs_path("$out_folder/PCC.median");

open my $pcc_median_fh, "<", $pcc_median_file
  or die "Cannot open PCC median file:$!\n";
readline $pcc_median_fh;
my %pcc_median;
while (<$pcc_median_fh>) {
    s/\r?\n//;
    my ( $window_size, $median ) = split /\t/;
    $pcc_median{$median} = $window_size;
}
my @medians = sort { $a <=> $b } keys %pcc_median;
my $choosed_windows_size = $pcc_median{ $medians[-1] };
warn "Choosed window size is : $choosed_windows_size\n";

########################################################################################
#Merging resutls
########################################################################################
warn "Outputing the final restuls...\n";
my @tad_bed_files = ();
foreach my $seq_id ( sort keys %seq_ids ) {
    next if $seq_id =~ /chrM|chrMT|MT/igx;
    my $tad_bed_file =
      $topdom_files{$seq_id}{bed}{$seq_id}{$choosed_windows_size};
    push @tad_bed_files, $tad_bed_file;
}

my $merge_bed_command = "cat @tad_bed_files> $out_folder/TADs.bed";
run_command($merge_bed_command);

sub run_command {
    my $command = shift;
    warn "\t$command\n";
    my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
      run( command => $command, verbose => 0 );
    if ($success) {
        warn "\tDone: $command!\n";
    }
    else {
        my @stderrs = @$stderr_buf;
        warn "Something went wrong:\n\t\t$command\n\t\t@stderrs";
    }
}

#===  FUNCTION  ================================================================
#         NAME: run_parallel
#      PURPOSE: given commands, run them in multiple threads
#   PARAMETERS: $process_num:
#				$missions: commonds
#      RETURNS: NA
#===============================================================================
sub run_parallel {
    my ( $process_num, @missions ) = @_;
    my $stream = Thread::Queue->new( @missions, undef );
    my $mission_num = scalar @missions;

    #assgn the task
    my @running = ();
    my @Threads;
    while ( @Threads < @missions ) {
        @running = threads->list(threads::running);

        if ( @running < $process_num ) {
            my $command = $stream->dequeue();
            my $thread = threads->new( \&run_command, $command );
            push( @Threads, $thread );
            my $tid = $thread->tid;
        }
        @running = threads->list(threads::running);
        foreach my $thr (@Threads) {
            if ( $thr->is_running() ) {
                my $tid = $thr->tid;
            }
            elsif ( $thr->is_joinable() ) {
                my $tid = $thr->tid;
                $thr->join;
            }
        }

        @running = threads->list(threads::running);
    }

    #join the threads
    while (@running) {
        foreach my $thr (@Threads) {
            $thr->join if ( $thr->is_joinable() );
        }
        @running = threads->list(threads::running);
        sleep(3);
    }
    return 0;
}

