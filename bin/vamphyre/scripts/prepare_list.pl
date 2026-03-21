#!/usr/bin/perl -w
################################################################################
#   Program Prepare_list.pl                                                    #
#   Description: Prepare files for VAMPhyRE analysis. Replaces FASTA header    #
#                with a single identifier and creates a list for renaming      #  
# Programmer: Alfonso Mendez-Tenorio                                           #
# Creation Date: June 18th, 2017                                               #
################################################################################
use strict;


my ($seqfile, $destfile, $filelist, $basename);

$basename = ">temp";
#$count = 0;


my $numargvs = $#ARGV + 1;
if ($numargvs != 2) {
    print "Prepare_file.pl: Prepares files for VAMPhyRE analysis\n";
    print "Notes: (1): Every file must contain a single genome sequence\n";
    print "       (2): The first argument must be the list of files to be processed\n";
    print "       (3): The second argument must be the list with new and old headers\n";
    print "       (4): The original files will be reheaded\n";
    print "\n";
    print "\nUsage: Prepare_list.pl list_of_files list_of_headers\n\n";
    exit;
}

$filelist = $ARGV[0];
$destfile = $ARGV[1];
#$filelist = "sample.txt";
#$destfile = "headers.txt";

unless (open (ARCH, $filelist))
    {
        print "The file $filelist cannot be open\n";
        print "Program terminated\n";
        exit;
    }

my @files = <ARCH>;



close ARCH;
#$seqfile = $ARGV[0];
#$destfile = $ARGV[1];

my $numfiles = scalar @files;

open (DEST, ">$destfile");
for (my $j = 0; $j<$numfiles; $j++) {
        $seqfile = $files[$j];
        chomp $seqfile;
        print "Processing $seqfile ... ";

        unless (open (ARCH, $seqfile))
        {
            print "The file $seqfile cannot be open\n";
            print "Program terminated\n";
            exit;
        }
        my @lines = <ARCH>;
        close ARCH;
         
         
        my $numlines = scalar @lines;
        my $header = $lines[0];
        chomp $header;
        my $newheader;
        if (($j > -1) && ($j<10)) {$newheader = $basename . '000' . $j} elsif
           (($j > 9) && ($j<100)) {$newheader = $basename . '00'  . $j} elsif
           (($j > 99) && ($j<1000)) {$newheader = $basename . '0'   . $j};   
        print DEST "$newheader;$header\n";
        print "reheaded...\n";

        open (ARCH, ">$seqfile");
        print ARCH  "$newheader\n";
        for (my$i=1;$i<$numlines; $i++){
            if ($lines[$i] =~ /^>/) {
               print ARCH "X\n";
            }  else {
               print ARCH "$lines[$i]";
            }   
        }
        close ARCH;
        
}
close DEST;

print "$numfiles files processed.\n";
print "\nFinished!\n";

exit;

