#!/usr/bin/env perl
# A script to calculate linear cleavage ratios using a bam file with PE reads as input

use strict;
use warnings;

use List::Util qw( max );
use POSIX qw(strftime);

# GMT formatted appropriately for your locale:
my $datestring = strftime "%a %b %e %H:%M:%S %Y", gmtime;

# print START time
printf("\nSTART date and time - $datestring\n\n");

## Define variables
# Variable to store the inferred LINEAR cleavage sites
my $lin_cleavage = undef;

# Output files
# Ratios
my $outfile_cov = 'linear_cleavage_ratio_PE_fast.cov.txt';
my $outfile_start = 'linear_cleavage_ratio_PE_fast.start.txt';
my $outfile_end = 'linear_cleavage_ratio_PE_fast.end.txt';
my $outfile_both = 'linear_cleavage_ratio_PE_fast.both.txt';
my $outfile_all = 'linear_cleavage_ratio_PE_fast.all.txt';

# Input
my $in_lin_file = shift; # $ARGV[0] is going to be the file you process, for example '/users/my_user/bam_files/my_sample_trimmed.sort.MT.q10.bam';
my $prefix = shift; # $ARGV[1] is going to be the name prefix for the file you are running. This prefix will also be the name of the folder with the results (e.g. my_cleavage_out)
my $outpath = $prefix;

# 1st check that everything is fine up to this point
print"IN: $in_lin_file\nPREFIX: $prefix\nOUTPATH: $outpath\n\n";

# Create output folder if not there
if (!-e $outpath) {mkdir($outpath)};

$outfile_cov = $prefix."_".$outfile_cov;
$outfile_start = $prefix."_".$outfile_start;
$outfile_end = $prefix."_".$outfile_end;
$outfile_both = $prefix."_".$outfile_both;
$outfile_all = $prefix."_".$outfile_all;

# End here if the results are there already
if (-s "$outpath/$outfile_all") {print "File $outpath/$outfile_all already created so no need to run it again\n";exit;};

# 2nd check that everything is fine up to this point
print "OUT_COV: $outfile_cov\nOUT_START: $outfile_start\nOUT_END: $outfile_end OUT_BOTH: $outfile_both\n\n";

print "LIN IN $in_lin_file\n\n";

## Read linear data
# Read only unique mappers determined by quality 255 to avoid NUMTS
open LINBAM, "samtools view -f 3 -q 255 $in_lin_file |" or die "Cannot open file to read: $!\n";

# 1 Find and store and count fragments starts and ends
# Loop through the file to get the starts and ends of the fragments and count them
while (my $line=<LINBAM>) {
	
	# Split data
	my @F = split "\t", $line;
	
	# Get the POS and length
	my ($ID,$pos,$pnext,$length) = @F[0,3,7,8];
	
	# Since this is PE you want to use the read mapped to forward strand as reference
	# therefore filter reads with negative or 0 length
	next if $length <= 0;	

	# Define start and end of the reads/fragments 
	my $start = $pos;
	my $end = $pos + $length -1;
	
	# Add the start and end to the count of cleavage sites
	$lin_cleavage->{'start'}{$start}++;
	$lin_cleavage->{'end'}{$end}++;
	
	# Store coverage from start to end
	for (my $range = $start; $range < $end+1; $range++) {
		
		$lin_cleavage->{'cov'}{$range}++;
	}	
}

## Close linear data
close LINBAM;

# Open file to print results
open OUTCOV, ">$outpath/$outfile_cov" or die "Cannot open file to write coverage: $!\n";
open OUTBOTH, ">$outpath/$outfile_both" or die "Cannot open file to write both: $!\n";
open OUTSTART, ">$outpath/$outfile_start" or die "Cannot open file to write start: $!\n";
open OUTEND, ">$outpath/$outfile_end" or die "Cannot open file to write ends: $!\n";
open OUTALL, ">$outpath/$outfile_all" or die "Cannot open file to write all: $!\n";

my $end_MT_chr = max keys %{$lin_cleavage->{'end'}};

# Calculate cleavage ratio for each of the mitochondrial chromosome positions
foreach my $pos (1..$end_MT_chr) {
	
	# Initialize variables for start,end and both ratios, etc
	my $st_ratio = undef;
	my $en_ratio = undef;
	my $both_ratio = undef;
	my $st_count = 0;
	my $en_count = 0;
	my $en_count_posMinus1 = 0;
	my $both_cov_count = 0;
	
	# END SHOULD BE THE POSITION YOU ARE CHECKING -1 (i.e. THE PREVIOUS POSITION)]
	# Define previous position
	my $posMinus1 = $pos-1;

	# Get coverage for this and previous position
	my $cov_count = ($lin_cleavage->{'cov'}{$pos} or 0);
	my $posMinus1_cov_count = ($lin_cleavage->{'cov'}{$posMinus1} or 0);
	
	# Get counts for the start and end positions 
	$st_count = ($lin_cleavage->{'start'}{$pos} or 0);
	$en_count = ($lin_cleavage->{'end'}{$pos} or 0);
	$en_count_posMinus1 = ($lin_cleavage->{'end'}{$posMinus1} or 0);

	# Calculate ratios
	# STARTS (#of starts overlapping pos/#of reads overlapping pos)
	$st_ratio = $st_count/$cov_count if $cov_count;
	
	# ENDS (#of ends overlapping pos/#of reads overlapping pos)
	$en_ratio = $en_count/$cov_count if $cov_count;
	
	# BOTH (#of starts overlapping pos + #of ends overlapping pos-1/#of reads overlapping pos+ #of reads overlapping pos-1)
	# !! HERE THE ends are calculated based on pos-1 !!! 
	# This way you calculate ratios of ends adjacent to starts (e.g in pos 2 end is 1 and start is 2)
	($both_ratio, $both_cov_count) = get_both_ratio($en_count_posMinus1 , $st_count ,$cov_count);

	# Print start, end, and both ratios and coverage	
	print OUTSTART "chrM\t$pos\t$pos\t$st_ratio\n";
	print OUTEND "chrM\t$pos\t$pos\t$en_ratio\n";
	print OUTBOTH "chrM\t$posMinus1\t$pos\t$both_ratio\n" if $posMinus1; # Check that first position is not 0!
	print OUTCOV "chrM\t$posMinus1\t$pos\t$both_cov_count\n" if $posMinus1; # Check that first position is not 0!
	print OUTALL "chrM\t$posMinus1\t$pos\t$both_ratio\t$st_ratio\t$en_ratio\t$both_cov_count\n" if $posMinus1; # Check that first position is not 0!
}

# Close nicely
close OUTCOV;
close OUTEND;
close OUTSTART;
close OUTBOTH;

# print END TIME
$datestring = strftime "%a %b %e %H:%M:%S %Y", gmtime;
printf("\n\nEND date and time - $datestring\n");

exit;

################################# SUBS #################################

sub get_both_ratio {
	# Get the list of counts
	my ($end_1_count , $start2_count ,$cov2_count)= @_;
	my $bridge_count = ($cov2_count - $start2_count);
	my $ratio =undef;
	if ($end_1_count + $start2_count + $bridge_count) {
		# Calculate ratio
		$ratio = ($end_1_count + $start2_count)/($end_1_count + $start2_count + $bridge_count);
	} else {
		# Or return 0 if there is no read coverage for those positions
		$ratio = 0;
	}

	return ($ratio, $bridge_count);
}


