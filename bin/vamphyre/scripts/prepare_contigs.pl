#!/usr/bin/perl -w
################################################################################
#   Program Prepare_contigs.pl                                                 #
#   Description: Joins the contig sequences contained in single file using a   #
#                character X as delimiter.                                     # 
# Programmer: Alfonso Mendez-Tenorio                                           #
# Creation Date: June 18th, 2017                                               #
################################################################################
use strict;
use File::Basename;

my ($seqfile, $destfile, $filelist);

my $numargvs = $#ARGV + 1;
if ($numargvs != 1) {
    print "Prepare_contigs.pl: Concatenates contigs contained in single files\n";
    print "Notes: (1): Every file must contain all the contigs to be joined\n";
    print "       (2): The character \"X\" is placed between contigs\n"; 
    print "       (3): The argument most be the list of files to be processed\n\n"; 
    print "\nUsage: Prepare_contigs.pl [list of files]\n\n";
    exit;
}

$filelist = $ARGV[0];

#$filelist = "listadraft.txt";

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

for (my $j = 0; $j<$numfiles; $j++) {
        $seqfile = $files[$j];
        chomp $seqfile;

        my $basename = basename($seqfile);
        my $pathfile = dirname($seqfile);
        print ("$basename\n");
        print ("$pathfile\n");
        
        #$destfile = "conc_" . $seqfile;
        $destfile = $pathfile . '/' . "conc_" . $basename; 
        print "Processing $seqfile ... saved as $destfile       ...OK\n";

        unless (open (ARCH, $seqfile))
        {
            print "The file $seqfile cannot be open\n";
            print "Program terminated\n";
            exit;
        }
        my @lines = <ARCH>;
        close ARCH;

        my $numlines = scalar @lines;

        open (ARCH, ">$destfile");
        print ARCH  "$lines[0]";
        for (my$i=1;$i<$numlines; $i++){
            if ($lines[$i] =~ /^>/) {
               print ARCH "X\n";
            }  else {
               print ARCH "$lines[$i]";
            }   
        }
        close ARCH;
}
print "$numfiles files processed.\n";
print "\nFinished!\n";

exit;
