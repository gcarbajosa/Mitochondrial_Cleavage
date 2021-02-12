#!/usr/bin/env perl
# A script to calculate the highest peaks within a threshold selected range of ratios
# It reads a file with the ratios and getd the subset of regions ("mountain ranges" or "himalayas") that are above the threshold

use strict;
use warnings;

use List::Util qw( max min );
use Getopt::Long qw(GetOptions);

# This are the defaults to test but you can call it with other options
my $ratiosFile = "test_ratiosFile.txt";
my $everestOutfile = 'test_everest.txt';
my $himalayaOutfile = 'test_himalaya.txt';
my $peakThres = 0.1;
my $gapAllowed = 5;
my $covThres = 20;

GetOptions('peakThres=i' => \$peakThres,
	   'covThres=i' => \$covThres,
	   'gapAllowed=i' => \$gapAllowed,
	   'ratiosFile=s'   => \$ratiosFile,
	   'everestOutfile=s'   => \$everestOutfile,
	   'himalayaOutfile=s'   => \$himalayaOutfile
	   ) or die ("Error in command line arguments\n");
      

## Read a file with ratios and store the regions above the selected threshold
my $aboveRegions = undef;

open RATIOSFILE, $ratiosFile or die "Cannot open file to read: $!\n";

while (my $line=<RATIOSFILE>) {
	chomp $line;
	print $line,"\n";
	my ($chr , $start , $end , $ratio, $starts_ratio , $ends_ratio , $coverage) = split "\t", $line;

	# Check here that the coverage is above what you want.
	if ($ratio >=  $peakThres and  $coverage >= $covThres) {
		# Store the set of ratios that are above the threshold
		$aboveRegions->{$start}=$ratio;
	}
}

# Be nice and close the file
close RATIOSFILE;

## Start by assigning regions (a.k.a. "himalayas")

# Initialize region counter
my $himalayasCount = 0;
my $himalayasRegions = undef;

# Get last position stored
my $lastPosStored = max(keys %$aboveRegions);

# Go one by one from 1 to last position stored
for (my $i = 1; $i <= $lastPosStored; $i++) {
	
	# if there is a peak in that position
	if ($aboveRegions->{$i}) {
		
		# check if there is another peak stored already within the range
		# Check the current range
		if ($himalayasRegions->{$himalayasCount}) {
			
			# If so assign the position to the range
			my $posMinusGap = $i-$gapAllowed;
			
			my $withinRangeFlag = 0;
			
			for (my $j = $posMinusGap; $j <= $i+$gapAllowed; $j++) {
				
				if ($himalayasRegions->{$himalayasCount}{$j}) {
					$himalayasRegions->{$himalayasCount}{$i}=$aboveRegions->{$i};
					# Move on to next position
					$withinRangeFlag = 1;
					;
				}
			}
			
			unless ($withinRangeFlag) {
				# If it was not within previous range create new one
				$himalayasCount++;
				$himalayasRegions->{$himalayasCount}{$i}=$aboveRegions->{$i};
			}
			
		} else {
			# or create a new range and assign the position if this is the first range
			$himalayasCount++;
			$himalayasRegions->{$himalayasCount}{$i}=$aboveRegions->{$i};
		}
		
	} else {
		# if there is no peak try the next one...
		next;
	}				
}

# Open files to write results
open(HIMALAYAFILE, '>', $himalayaOutfile) or die "Cannot open file to read: $!\n";
open(EVERESTFILE, '>', $everestOutfile) or die "Cannot open file to read: $!\n";

## Get the highest peak that represent the mountain range ("everest" for the given "himalaya") 
# Go through the stored ranges
for my $range (sort {$a <=> $b} keys %{$himalayasRegions}) {

	my ($range_start) = min(keys $himalayasRegions->{$range});
	my ($range_end) = max(keys $himalayasRegions->{$range});
	
	my $top_ratio = 0;
	
	# Get the maximum peak for the range and print it out
	my ($everest_key, $everest_ratio) = ();
	
	for my $peak (sort {$a <=> $b} keys %{$himalayasRegions->{$range}}) {
		
		if ($top_ratio < $himalayasRegions->{$range}{$peak}) {
			
			$everest_key = $peak;
			$everest_ratio = $himalayasRegions->{$range}{$peak};
			$top_ratio = $everest_ratio;
		} else {

			next;
		}
	}
	
	# Print Himalaya with top ratio in bed format
	print HIMALAYAFILE "chrM\t$range_start\t$range_end\t$everest_ratio\n";
	
	# Print everest with top ratio in bed format
	print EVERESTFILE "chrM\t$everest_key\t$everest_key\t$everest_ratio\n";
	
}

# Be nice and close the file
close HIMALAYAFILE;
close EVERESTFILE;

# That's it!
exit;
