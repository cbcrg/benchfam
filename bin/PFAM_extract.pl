#! /usr/bin/perl
# 
# PFAM_extract.pl read the Pfam database Pfam-A.full and extract for each dataset the sequences
# with an existing structure. Then it will generate a fasta file in fasta format for each of 
# them.

use strict;
use File::Copy;

my $ia;
my $ib;
my $ic;
my $id;
my $ie;
my $if;

my @seq_all;
my @seq_pdb;
my @seq;

open(PFAM,"Pfam-A.full");

while (my $line=<PFAM>)
{
	chomp $line;
	if ($line=~/.*GF\sAC.*/)
	{
		$ia++;
		if ($ib ne $ia)
		{
			@seq_all="";
			@seq_pdb="";
			$id=0;
			close (SET);	
		}
		$line=~/(\w+)(\s+)AC(\s+)(\w+)/;
		print "$4 \n";
		$ib=$ia;
		open (SET,">$4.fasta");
	}
	if ($line=~/.*DR\sPDB.*/)
	{
                $line=~/(\w+)(\s+)(\w+)\/(\w+)-(\w+)(\s+)/;
		my $name="$3/$4-$5";
		for ($ic=0;$ic<=$id;$ic++)
		{
			if ($name eq @seq_pdb[$ic]) { $ie++;}
		}
		if ($ie eq 0)
		{
			push (@seq_pdb,"$3/$4-$5");
			$id++;
		}
		else {$ie=0;}
	}
	if (substr ($line,0,1) ne '#')
	{
		$line=~/(\w+)\/(\w+)-(\w+)(\s+)(\S+)/;	
		for ($if=0;$if<=$id;$if++)
		{
			if ($seq_pdb[$if] eq "$1/$2-$3")
			{
				print SET ">$1/$2-$3\n";
				$line=~/(\w+)\/(\w+)-(\w+)(\s+)(\w+)/;
				my $buffer=$5;
				$buffer=~s/\.||-//g;
				$buffer=~tr/[a-z]/[A-Z]/;
				print SET "$buffer\n";
			}
		}
	}
}

close(PFAM);
