#! /usr/bin/perl
#
####################################################################################################
# Read the template file returned by nextflow and get the corresponding PDB files to be extracted.
# Then use t-coffee reformating to extract the corresponding chain ID.
# ARG[0] : the trace file coming from BenchFam
# ARG[1] : the PATH for the scratch directory
####################################################################################################

use strict;
use File::Copy;

my $file = $ARGV[0];
my $path = $ARGV[1];

my @scratch1;
my @scratch2;
my @scratch3;
my @scratch4;
my @scratch5;
my @scratch6;
my @scratch7;
my @scratch8;
my @pfam1;
my @pfam2;
my @step;
my @status;

my $temp1;
my $temp2;
my $temp3;
my $temp4;
my $temp5;
my $temp6;
my $temp7;
my $temp8;

my $pfam_tmp;

my $ia;
my $ib;
my $ic;
my $id;

open (TEMP,"$file");
open (SHLL,">run_copy_scratch.sh");
while (my $line1=<TEMP>)
{
	$ia++;
	$line1=~/(\w+)(\s+)(\S+)(\s+)(\w+)(\s+)(\w+)(\s+)(\S+)(\s+)(\w+)/;
#        $pfam_tmp="./".substr($9,1,7);
        $pfam_tmp=substr($9,1,7);
	$temp1="$path"."$3"."*/*irmsd";
	$temp2="$path"."$3"."*/*aln";
	$temp3="$path"."$3"."*/*lib";
	$temp4="$path"."$3"."*/3Dmcoffee.aln";
	$temp5="$path"."$3"."*/*pdb";
	$temp6="$path"."$3"."*/super.pml";
	$temp7="$path"."$3"."*/modified.fasta";
	$temp8="$path"."$3"."*/modified.template";


	push(@scratch1,"$temp1");
	push(@scratch2,"$temp2");
	push(@scratch3,"$temp3");
	push(@scratch4,"$temp4");
	push(@scratch5,"$temp5");
	push(@scratch6,"$temp6");
	push(@scratch7,"$temp7");
	push(@scratch8,"$temp8");
	push(@pfam1,"$pfam_tmp");
	push(@step,"$7");
	push(@status,"$11");

	if ($7 eq '3_align_libraries' &&  $11 eq 'COMPLETED')
	{	
		$ic++;
#		print "$pfam_tmp  $7 \n";
		push (@pfam2,"$pfam_tmp");	
	}
}

for ($ib=0;$ib<$ia;$ib++)
{
	for ($id=0;$id<$ic;$id++)
	{ 
#		print "$pfam1[$ib] / $pfam2[$id] / $step[$ib] / $scratch1[$ib] \n";
#		if ($pfam1[$ib] eq $pfam2[$id] && $step[$ib] eq '2_extract_PDB') 
		if ($pfam1[$ib] eq $pfam2[$id])
		{
			if ($step[$ib] eq '2_extract_PDB')
			{
				print SHLL "mkdir BENCHFAM/$pfam1[$ib] \n";
				print SHLL "mkdir BENCHFAM/$pfam1[$ib]/PDB \n";
				print SHLL "cp $scratch5[$ib] BENCHFAM/$pfam1[$ib]/PDB \n";
				print SHLL "cp $scratch6[$ib] BENCHFAM/$pfam1[$ib]/PDB \n";
				print SHLL "cp $scratch7[$ib] BENCHFAM/$pfam1[$ib]/$pfam1[$ib]_test.fa \n";
				print SHLL "cp $scratch8[$ib] BENCHFAM/$pfam1[$ib]/$pfam1[$ib]_ref.template_list \n";
			}

			if ($step[$ib] eq '3_align_libraries')
			{
				print SHLL "mkdir BENCHFAM/$pfam1[$ib]/EVAL \n";
				print SHLL "mkdir BENCHFAM/$pfam1[$ib]/ALN \n";
				print SHLL "mkdir BENCHFAM/$pfam1[$ib]/LIB \n";
				print SHLL "cp $scratch1[$ib] BENCHFAM/$pfam1[$ib]/EVAL \n";
				print SHLL "cp $scratch2[$ib] BENCHFAM/$pfam1[$ib]/ALN \n";
				print SHLL "cp $scratch3[$ib] BENCHFAM/$pfam1[$ib]/LIB \n";
				print SHLL "cp $scratch4[$ib] BENCHFAM/$pfam1[$ib]/$pfam1[$ib]_ref.aln \n";
			}
		}
	}
}

close (TEMP);
close (SHLL);
