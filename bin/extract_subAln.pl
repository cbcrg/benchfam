#!/bin/env perl
($subSeq_file, $largeSeq_file)=@ARGV;


if($subSeq_file=~/.*\/(PF.*)\.sp_lib/){  $subSeq_file_name=$1; }
if($largeSeq_file=~/.*\/(.*)\.fa/){  $largeSeq_file_name=$1; }

$out_file=$subSeq_file_name."_".$largeSeq_file_name.".fa";
$erout_file="$subSeq_file_name"."_error.log";
open OUT, ">$out_file";
open EROUT,">>$erout_file";

$i=0;
%hashStru=(); 
open SUBseq, $subSeq_file;
while (<SUBseq>)
{ 
    if($_=~/^#.*/){ last; }
    if($_=~/^(.*)\s\d+.*/){ 
	$hashStru{$1}=$i;
	$i++;
    } 
}
#foreach $key( sort { $hashStru{$a} <=> $hashStru{$b} } (keys %hashStru) ){  print $key."\t".$hashStru{$key}."\n";   }

$/=">";
%hashFinal=(); 
open LARGEseq, $largeSeq_file;
while (<LARGEseq>)
{   
    $entry=$_;
    chop $entry;
    $entry= ">"."$entry";
    $entry=~/>(.+?)\n(\C*)/g;
    $title=$1;$sequence=$2;
    #$sequence=~s/\n//g;
   
    if ($title ne "")
    {   
	if( exists $hashStru{$title} ){  $hashFinal{$title}=$sequence;  }
    }
}
$/="\n";

foreach $key( sort { $hashStru{$a} <=> $hashStru{$b} } (keys %hashStru) ){ print OUT ">".$key."\n".$hashFinal{$key};  }

##### --- NOT NECESSARY! This is to check if all seqs in lib exist in the large scale MSA --- #####
foreach $key( sort { $hashStru{$a} <=> $hashStru{$b} } (keys %hashStru) ){ 
  if( !exists $hashFinal{$key}) { print EROUT "This sequence didn't exist in the large scale MSA: $key \n"; }
}
