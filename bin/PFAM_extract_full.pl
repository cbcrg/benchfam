#!/usr/bin/perl
# 
# PFAM_extract.pl read the Pfam database Pfam-A.full and extract for each dataset the sequences
# with an existing structure. Then it will generate a fasta file in fasta format for each of 
# them.

use File::Copy;

$ia=0;
$ib;
$ic;
$id1;
$id2;
$ie;
$if;
$ig;
$ih;

@seq_all;
@seq_pdb;
$name_pdb;
$name;
$option=$ARGV[0];
$limits=$ARGV[1];
$file=$ARGV[2];
$mode;

if ($option eq "FULL") { $mode=1;}
if ($option eq "PDB")  { $mode=2;}

open(PFAM, "< " . $file);

while ($line=<PFAM>)
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
		if ($mode eq  "1") { open (SET_ALL,">$4_full.fasta");}
		if ($mode eq "2")  { open (SET_PDB,">$4_pdb.fasta");}
	}

	if ($mode eq  "1")
	{
		if ($line=~/#=GS/ && $line!~/.*DR\sPDB.*/)
		{
                	$line=~/(\w+)(\s+)(\w+)\/(\w+)-(\w+)(\s+)/;
			push (@seq_all,"$3/$4-$5");
			$id1++;
		}
	        if (substr ($line,0,1) ne '#')
		{
			$line=~/(\w+)\/(\w+)-(\w+)(\s+)(\S+)/;
			for ($ie=0;$ie<=$id1;$ie++)
			{
				if ($seq_all[$ie] eq "$1/$2-$3")
				{
                        $name="$1/$2-$3";	
                        $line=~/(\w+)\/(\w+)-(\w+)(\s+)(\w+)/;
                        $buffer=$5;
                        $buffer=~s/\.||-//g;
                        $buffer=~tr/[a-z]/[A-Z]/;
                        $Bindex = index($buffer, 'B'); 
                        $Zindex = index($buffer, 'Z');
                        $Xindex = index($buffer, 'X');
                        if( $Bindex == -1 && $Zindex==-1 && $Xindex==-1){
                            print SET_ALL ">$1/$2-$3\n";
                            print SET_ALL "$buffer\n";
                        }
			    }
			}
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
			$line=~/(\w+)\/(\w+)-(\w+)(\s+)(\S+)/;
			for ($ie=0;$ie<=$id2;$ie++)
			{
				if ($seq_pdb[$ie] eq "$1/$2-$3")
				{	
					$line=~/(\w+)\/(\w+)-(\w+)(\s+)(\w+)/;
					$buffer=$5;
					$buffer=~s/\.||-//g;
					$buffer=~tr/[a-z]/[A-Z]/;
					$Bindex = index($buffer, 'B'); 
                    $Zindex = index($buffer, 'Z');
                    $Xindex = index($buffer, 'X');
                    if( $Bindex == -1 && $Zindex==-1 && $Xindex==-1){
                       print SET_PDB ">$1/$2-$3\n";
                       print SET_PDB "$buffer\n";
                    }
				}
			}
		}					
	}
}

close(PFAM);

