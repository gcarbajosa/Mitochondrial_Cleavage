#!/usr/bin/perl
# A script to generate a matrix with the clevage ratios to be used as input in, for example, quantitaive trait loci analysis with PLINK
# Requires
# ARGV[0] filename of the peaks presence positions table (required)
# ARGV[1] filename with a list of ratio files you want to use as input (required)
# ARGV[2] threshold for the proportion you want to use as threshold (optional, default is 0.5)
# ARGV[3] filename for the ouput file (optional, default is peaks_presence_matrix.05.txt)
# $ARGV[4] read coverage (optional, default is 20)

use strict;
use warnings;

my $posTable = $ARGV[0] or "peaks_presence_full_table.txt";
my $ratio_filenames_list = $ARGV[1] or "ratio_filenames_list.txt";
my $threshold = $ARGV[2] or 0.5; # I am using 0.5 as default
my $outfile = $ARGV[3] or "peaks_presence_matrix.05.txt";
my $readCoverage = $ARGV[4] or 20;

# Read list of files into a hash
my $filelist = undef;
my $peakCoverage = undef;
my $posToKeep =undef;

# Load the table using a threshold to keep the positions you want
# Open the peaks presence table and filter by the specificed threshold 
open POS_TABLE , $posTable;
while (<POS_TABLE>) {
	
	my @F = split /\s/;
	
	next unless $F[2] > $threshold;
	$posToKeep->{$F[0]}++;
	#print;
	
}
close POS_TABLE;

# Open file to write results
open WRITE , ">$outfile";

# Print the positions to head the file
for my $pos (sort {$a <=> $b} keys %$posToKeep ) {

	print WRITE "\t$pos";
	
}
print WRITE "\n";

# Load all the ratio file names
open FILELIST , $ratio_filenames_list;
while (<FILELIST>) {

	chomp;
	my @F = split /\s/;
	
	$filelist->{$F[0]}=1;	
}

close FILELIST;

#  Loop through the filenames list to fill in the matrix with the requires cleavage ratios
for my $rat_filename (keys %$filelist ) {
	
	print "RATFILE: ",$rat_filename,"\n";
	
	print WRITE $file;
	
	open RATFILE , $rat_filename or die "Cannot open $rat_filename:$!";
	
	while (<RATFILE>) {
		# print to check
		#print;
		chomp;
		my @F = split /\s/;

		print "$_\n" if $posToKeep->{$F[1]};
		
		if ($posToKeep->{$F[1]}) {
			# Store if the ratio is there and if the coverage is good enough
			if ($F[-1] >= $readCoverage) {
				print WRITE "\t$F[3]";}
			else {
				# else print a blank
				print WRITE "\t-9";
			}
		}
		
	}
	
	close RATFILE;
	
	print WRITE "\n";
}

close WRITE;

exit;

