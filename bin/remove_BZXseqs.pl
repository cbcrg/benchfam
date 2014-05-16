#!/usr/bin/perl
use POSIX;

($aln_file, $ref_file, $out_file)=@ARGV;
# SOS: the output name is the same as the input name. 
# SOS2: to input file prepei na einai se absolute path
#if($aln_file=~/.*\/(.*\.fa)/){  $out_file_name=$1; } 
#else{  print "Unable to read input file name!!! Most probably doesn't end with .fa or is not given in it's absolute path!! "; exit;  }

#_________ Read Ref MSA file and store ref. seqs _________#

$/=">";
%REFhash=(); 
open REFseq, $ref_file;
while (<REFseq>)
{   
    $entry=$_;
    chop $entry;
    $entry= ">"."$entry";
    $entry=~/>(.+?)\n(\C*)/g;
    $title=$1;$sequence=$2;
    $sequence=~s/\n//g;
   if ($title ne "" && $sequence ne ""){   $REFhash{$title}=$sequence;  }
}
$/="\n"; print "Reading of Ref. file DONE!!\n";
close(REFseq);

#_________ Read MSA file and create column hashes _________#

$/=">";
%ALNhash=(); 
%SEQhash=();
$SeqNum=0;
open ALNseq, $aln_file;
open OUT, ">>", $out_file or die $!; 
while (<ALNseq>)
{   
    $entry=$_; 
    chop $entry;
    $entry= ">"."$entry"; 
    $entry=~/>(.+?)\n(\C*)/g;
    $title=$1;$sequence=$2; $sequence=~s/\n//g;
    if ( ($title ne "" && $sequence ne "") )
    {   
	$Bindex = index($sequence, 'B'); 
        $Zindex = index($sequence, 'Z');
	$Xindex = index($sequence, 'X');
	if( $Bindex == -1 && $Zindex==-1 && $Xindex==-1){
	  print OUT ">".$title."\n".$sequence."\n";
	  $ALNhash{$title}=$sequence;
	}
	else{
	     if (exists $REFhash{$title}){
		 print OUT ">".$title."\n".$sequence."\n";
		 $ALNhash{$title}=$sequence;
	    }
	}
    }
} 
$/="\n"; print "Reading of MSA file DONE!!\n"; 
close(ALNseq);
close(OUT);
