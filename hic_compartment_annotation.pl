#!/usr/bin/perl -w

use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use Cwd qw[abs_path];
use IPC::Cmd qw[can_run run];
use Scalar::Util qw[looks_like_number];
use threads;
use Thread::Queue;
# use Data::Printer;

my $usage = <<USAGE;
SYSNOPSIS
hic_compartmente_annotation.pl [options] A.bed B.bed C.bed

 Options:
   -a --assembly      assembly, default mm10
   -b --bin           bin-size, default 100000
   -c --column        signal column used for retrive
   -f --folder        use folder to infer basename
   -o --outfolder     output folder
   -p|--progress      progress, default 16
   -h --help          print this help.
USAGE

my $out_folder  = './';
my $assembly    = 'mm10';
my $folder_flag = 0;
my $col         = 5;
my $progress    = 16;
my $bin_size    = 100000;
my $help        = 0;
die $usage
  unless GetOptions(
    "a|assembly:s"  => \$assembly,
    "b|bin=i"       => \$col,
    "o|outfolder:s" => \$out_folder,
    "f|folder"      => \$folder_flag,
    "c|column=i"    => \$bin_size,
    "p|progress=i"  => \$progress,
    "help"          => \$help,
  );
die $usage if $help;
my @files = @ARGV;
die $usage unless @files;

######################################
#parameters
######################################
my @beds = ();

#filter the input beds
foreach my $glob_str (@ARGV) {

    #pick up directory
    my @bed_tmp = map { abs_path($_) } grep { -s $_ } glob $glob_str;

    #pick up beds
    push @beds, @bed_tmp if @bed_tmp;
}

die "No beds were found, please check the the glob string\n"
  unless @beds;

#Calculate the basename
my %bed_files   = ();
my %tmp_cell    = ();
my %tmp_markers = ();
my %marker      = ();
foreach my $bed ( sort @beds ) {
    my $basename = $bed;
    if ($folder_flag) {
        my @paths = split /\//, $bed;
        $basename = $paths[-2];
    }
    else {
        $basename = basename( $bed, ".bed" );
    }
    my ( $cell_type, $marker ) = split /_|\./, $basename;
    $bed_files{$cell_type}{$marker} = abs_path($bed);

    $tmp_cell{$cell_type} = 1;
    $tmp_markers{$marker} = 1;
}
my @cell_types = sort keys %tmp_cell;
my @markers    = sort keys %tmp_markers;

$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;

######################################
#BinarizedSignal
######################################
warn "Binarize compartment A/B\n";
my %data;
foreach my $cell_type ( sort keys %bed_files ) {
    foreach my $marker ( sort keys %{ $bed_files{$cell_type} } ) {
        my $bed_file = $bed_files{$cell_type}{$marker};
        open my $bed_fh, "<", $bed_file
          or die "Cannot open bed file:$bed_file\n@!";

        my $last_seq_id = '';
        my $wid         = 0;
        while (<$bed_fh>) {
            next if /^#/;
            s/\r?\n//;
            my @tmp    = split /\t/;
            my $seq_id = $tmp[0];
            my $score  = $tmp[ $col - 1 ];
            if ( $last_seq_id eq $seq_id ) {
                $wid++;
            }
            else {
                $wid = 1;
            }
            if ( $score eq '.' || $score eq "NA" || $score eq 'NaN' ) {
                $score = 2;
            }
            else {
                $score = $score > 0 ? 1 : 0;
            }
            $data{$cell_type}{$seq_id}{$wid}{$marker} = $score;

            $last_seq_id = $seq_id;
        }
        close $bed_fh;
    }
}

#Outputting BinarizedSignal
my $bi_signal_folder = $out_folder . "/BinarizedSignal";
mkdir $bi_signal_folder unless -d $bi_signal_folder;
foreach my $cell_type ( sort keys %data ) {
    foreach my $seq_id ( sort keys %{ $data{$cell_type} } ) {

        my $data_file =
          $bi_signal_folder . "/${cell_type}_${seq_id}._binary.txt";
        open my $data_fh, ">", $data_file
          or die "Cannot open file:$data_file\n@!";

        #print header
        print $data_fh "$cell_type\t$seq_id\n";
        my @markers = sort keys %{ $bed_files{$cell_type} };
        print $data_fh join "\t", @markers;
        print $data_fh "\n";

        foreach
          my $wid ( sort { $a <=> $b } keys %{ $data{$cell_type}{$seq_id} } )
        {
            my @scores = ();
            foreach my $marker (@markers) {
                my $score = $data{$cell_type}{$seq_id}{$wid}{$marker};
                push @scores, $score;
            }
            print $data_fh join "\t", @scores;
            print $data_fh "\n";
        }

    }
}

######################################
#Segmentation
######################################
warn "Segregating the files\n";
my @learn_commands = ();
my $seg_folder     = "$out_folder/Segmentation/";
mkdir $seg_folder unless -e $seg_folder;

foreach my $num_stat ( 2 .. 2**@beds ) {
    my $command =
"ChromHMM LearnModel -b $bin_size $bi_signal_folder $seg_folder $num_stat $assembly";
    push @learn_commands, $command;
}
run_parallel( $progress, @learn_commands );
######################################
#StatePruning
######################################
my $prun_folder           = " $out_folder/StatePruning/";
my $state_pruning_command = "ChromHMM StatePruning $seg_folder $prun_folder";
run_command($state_pruning_command);

######################################
#Assign the status
######################################
my $assign_folder = $out_folder . "/Assignment";
mkdir $assign_folder unless -d $assign_folder;
my @assign_commands = ();
foreach my $cell_type ( sort keys %bed_files ) {
    foreach my $marker ( sort keys %{ $bed_files{$cell_type} } ) {
        my $bed = $bed_files{$cell_type}{$marker};
        foreach my $num_stat ( 2 .. 2**@beds ) {
            my $bed_results =
              "$assign_folder/${cell_type}_${marker}_${num_stat}.bed";
            my $bed_input = "$seg_folder/${cell_type}_${num_stat}_segments.bed";
            my $command =
"bedtools intersect -a $bed -b $bed_input -wo | cut -f 1-8,12 > $bed_results";
            push @assign_commands, $command;
        }
    }
}

run_parallel( $progress, @assign_commands );

######################################
#Heatmap
######################################
my $heatmap_folder = $out_folder . "/Heatmap";
mkdir $heatmap_folder unless -d $heatmap_folder;

#Store the Data;
my %heatmap     = ();
my %tmp_bin_ids = ();
foreach my $cell_type ( sort keys %bed_files ) {
    foreach my $marker ( sort keys %{ $bed_files{$cell_type} } ) {
        foreach my $num_stat ( 2 .. 2**@beds ) {
            my $bed_in =
              "$assign_folder/${cell_type}_${marker}_${num_stat}.bed";
            open my $bed_fh, "<", $bed_in
              or die "Cannot open file $bed_in:$!";
            while (<$bed_fh>) {
                s/\r?\n//;
                my @tmp    = split /\t/;
                my $id     = $tmp[3];
                my $eigen  = $tmp[4];
                my $status = $tmp[8];
                $heatmap{$cell_type}{$marker}{$num_stat}{$id}{Eigen} =
                  $eigen eq '.' ? 'NaN' : $eigen;
                $heatmap{$cell_type}{$marker}{$num_stat}{$id}{Status} = $status;
                $tmp_bin_ids{$id} = 1;

            }
        }
    }
}
my @bin_ids = sort keys %tmp_bin_ids;

#Output the Data;

#Total
foreach my $num_stat ( 2 .. 2**@beds ) {

    open my $total_fh, ">", $heatmap_folder . "/Total.$num_stat.txt";
    my @labels = ("ID");
    foreach my $cell_type (@cell_types) {
        foreach my $marker (@markers) {
            push @labels, $cell_type . "_" . $marker;
        }
    }

    push @labels, "Status";
    print $total_fh join "\t", @labels;
    print $total_fh "\n";

    foreach my $id (@bin_ids) {
        my %status;
        my @lines = ($id);
        foreach my $cell_type (@cell_types) {
            foreach my $marker (@markers) {
                my $eigen =
                  $heatmap{$cell_type}{$marker}{$num_stat}{$id}{Eigen};
                my $status =
                  $heatmap{$cell_type}{$marker}{$num_stat}{$id}{Status};
                push @lines, $eigen;
                $status{$status} = 1;
            }
        }
        my @status = keys %status;
        p @status if @status > 1;
        die "Status errer on $id\n" if @status > 1;
        push @lines, $status[0];

        print $total_fh join "\t", @lines;
        print $total_fh "\n";
    }
    close $total_fh;

}

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

