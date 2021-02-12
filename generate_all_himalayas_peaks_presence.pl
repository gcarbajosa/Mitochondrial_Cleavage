#!/usr/bin/perl

use strict;
use warnings;

my $outfile = "peaks_presence_full_table.txt";

# Read list of files into a hash
my $fileListFilename = $ARGV[0];

# Open a file containing the list of files you want to process
# For example:
# eve1.txt
# eve2.txt
# ...

# Define variables
my $filelist = undef;
my $peakCoverage = undef;
my $filecount = 0;

# Load all the file names
open FILELIST , "$fileListFilename";
while (<FILELIST>) {

	chomp;
	my @F = split /\s/;	
	$filelist->{$F[0]}=1;	
}
close FILELIST;

# Go through all the files in the filelist
for my $him_filename (keys %$filelist ) {
	
	$filecount++;
	
	open EVEFILE , $him_filename or die "Cannot open $him_filename:$!";
	
	while (<EVEFILE>) {
		
        # Split the fields
		my @F = split /\s/;
		
		# Add to the count of hits on the mitochondrial chromosome position
		$peakCoverage->{$F[1]}{'count'}++;
	}
}

# Open file to write results
open WRITE , ">$outfile";

# Check the coverage for each positions
for my $peak (sort {$a <=> $b} keys %$peakCoverage ) {

	my $peakCount = $peakCoverage->{$peak}{'count'};
	my $peakCoverageRatio = $peakCount/$filecount;
	print WRITE "$peak\t$peakCount\t$peakCoverageRatio\t$filecount\n";
	
}

close WRITE;

exit; 

