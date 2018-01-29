#!/usr/bin/perl
# 
###################################################################################################
# PFAM_extract.pl read the Pfam database Pfam-A.full and extract for each dataset the sequences
# with an existing structure. Then it will generate a fasta file in fasta format for each of 
# them.
# Argument 1 : decide to extract all sequences or only PDB sequences ("FULL" or "PDB")
# Argument 2 : number max of family to extract (anything but a number extracts all) 
# Argument 3 : input name of the PFAM format database to be extracted
####################################################################################################

use strict;
use File::Copy;

my $ia=0;
my $ib;
my $ic;
my $id1;
my $id2;
my $ie;
my $if;
my $ig;
my $ih;

my @seq_all;
my @seq_pdb;
my $name_pdb;
my $name;
my $option=$ARGV[0];
my $limits=$ARGV[1];
my $file=$ARGV[2];
my $mode;

if ($option eq "FULL") { $mode=1;}
if ($option eq "PDB")  { $mode=2;}

open(PFAM,"$file");

while (my $line=<PFAM>)
{
	last if ($ia eq  $limits);
	chomp $line;
	if ($line=~/.*GF\sAC.*/)
	{
		$ia++;
		if ($ib ne $ia)
		{
			@seq_all="";
			@seq_pdb="";
			$id1=0;
			$id2=0;
			close (SET_ALL);	
			close (SET_PDB);	
		}
		$line=~/(\w+)(\s+)AC(\s+)(\w+)/;
		print "FAMILY: $4 \n";
		$ib=$ia;
		if ($mode eq "1")  
		{ 
			open (SET_ALL,">$4_full.aln");
			open (SET_ALL_FA,">$4_full.fa");
		}
		if ($mode eq "2")  { 
			open (SET_PDB,">$4_pdb.aln");
			open (SET_PDB_FA,">$4_pdb.fa");
		}
	}

	if ($mode eq  "1")
	{
	        if (substr($line,0,1) ne '#' && substr($line,0,2) ne '//')
		{
#			print "$line\n";
	                $line=~/(\w+)\/(\w+)-(\w+)(\s+)(\D+)/;
			print SET_ALL ">$1/$2-$3\n";
			print SET_ALL_FA ">$1/$2-$3\n";
			my $buffer=$5;
			$buffer=~tr/\./\-/;
			$buffer=~tr/[a-z]/[A-Z]/;
			print SET_ALL "$buffer\n";
			$buffer=~s/\.||-//g;
			print SET_ALL_FA "$buffer\n";
		}
	}

	if ($mode eq "2")
	{
        	if ($line=~/#=GS/ && $line=~/.*DR\sPDB.*/)
        	{
               		$line=~/(\w+)(\s+)(\w+)\/(\w+)-(\w+)(\s+)/;
			$name_pdb="$3/$4-$5";
			for ($if=0;$if<=$id2;$if++)
			{
				if ($name_pdb eq $seq_pdb[$if]) { $ig++; }
			}
			if ($ig eq 0)
			{
               			push (@seq_pdb,"$3/$4-$5");
				$id2++;
			}
			else { $ig=0;}
        	}
		if (substr ($line,0,1) ne '#')
		{
			$line=~/(\w+)\/(\w+)-(\w+)(\s+)(\D+)/;
			for ($ie=0;$ie<=$id2;$ie++)
			{
				if ($seq_pdb[$ie] eq "$1/$2-$3")
				{	
					print SET_PDB ">$1/$2-$3\n";
					print SET_PDB_FA ">$1/$2-$3\n";
					$line=~/(\w+)\/(\w+)-(\w+)(\s+)(\D+)/;
					my $buffer=$5;
					$buffer=~tr/\./\-/;
					$buffer=~tr/[a-z]/[A-Z]/;
					print SET_PDB "$buffer\n";
					$buffer=~s/\.||-//g;
					print SET_PDB_FA "$buffer\n";
				}
			}
		}					
	}
}

close (SET_ALL);
close (SET_ALL_FA);
close (SET_PDB);
close (SET_PDB_FA);
close(PFAM);
