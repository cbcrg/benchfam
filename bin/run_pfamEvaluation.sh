#!/bin/bash
export PATH=/users/cn/mhatzou/tcoffee_source/t_coffee/src:$PATH
unset MAFFT_BINARIES


if [ ! -d "fam_largeScaleData_MSA" ]
then
    mkdir -p fam_largeScaleData_MSA
fi

largeScalePath="/users/cn/mhatzou/Datasets/PFAM_25/fam_largeScaleData_Cedric"   #where largeScale families Data are  
out_dir="/users/cn/mhatzou/projects/pfam_evaluation/fam_largeScaleData_MSA"     #where largeScale MSA will be output  



usage="$(basename $0) -s n [-h] [-t n]  -- program to run run_pfamEvaluation.sh

where:
    -h  show this help text
    -t  sets the number of threads to use
    -s  step 1=Create MSA for families' largeScale datasets, 2=Create sp_library, 3=Extract the MSA of the seqs in the lib from the respective largeScale MSA, 4=Evaluate, using as reference your sp_library. (steps 1,2 can be running at the same time)"

n_threads=1
step=0

while getopts ":t:s:h" Option
# Initial declaration.
# a, b, c, d, e, f, and g are the options (flags) expected.
# The : after option 'e' shows it will have an argument passed with it.
do
	case $Option in
		h) echo "$usage"
			exit
			;;
		t)
			n_threads=$OPTARG
			;;
		s)
			step=$OPTARG
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1
			;;
	esac
done
shift $(($OPTIND - 1))


famNames="PF00004 PF00005 PF00006 PF00008 PF00013 PF00014 PF00016 PF00017 PF00018 PF00019 PF00020 PF00023 PF00025 PF00026 PF00027 PF00028 PF00030 PF00034 PF00036 PF00037 PF00040 PF00041 PF00042 PF00043 PF00044 PF00045 PF00046 PF00047 PF00048 PF00051 PF00056 PF00059 PF00061 PF00062 PF00067 PF00068 PF00069 PF00070 PF00071 PF00072 PF00073 PF00074 PF00075 PF00076 PF00077 PF00079 PF00080 PF00081 PF00082 PF00084 PF00085 PF00087 PF00089 PF00092 PF00096 PF00102 PF00104 PF00105 PF00106 PF00107 PF00109 PF00111 PF00112 PF00113 PF00116 PF00117 PF00121 PF00125 PF00127 PF00128 PF00129 PF00132 PF00134 PF00139 PF00141 PF00144 PF00147 PF00149 PF00150 PF00152 PF00155 PF00156 PF00160 PF00161 PF00162 PF00167 PF00168 PF00169 PF00171 PF00175 PF00179 PF00180 PF00185 PF00186 PF00187 PF00190 PF00191 PF00194 PF00197 PF00199 PF00200 PF00202 PF00205 PF00206 PF00208 PF00210 PF00211 PF00215 PF00217 PF00227 PF00228 PF00229 PF00230 PF00232 PF00233 PF00235 PF00238 PF00240 PF00241 PF00246 PF00248 PF00254 PF00255 PF00258 PF00266 PF00268 PF00270 PF00271 PF00285 PF00288 PF00289 PF00291 PF00293 PF00294 PF00296 PF00300 PF00303 PF00306 PF00307 PF00318 PF00326 PF00328 PF00331 PF00334 PF00337 PF00342 PF00347 PF00348 PF00352 PF00355 PF00370 PF00373 PF00378 PF00383 PF00384 PF00385 PF00389 PF00392 PF00394 PF00400 PF00405 PF00406 PF00407 PF00412 PF00413 PF00415 PF00431 PF00432 PF00435 PF00439 PF00440 PF00441 PF00445 PF00448 PF00452 PF00456 PF00457 PF00459 PF00462 PF00463 PF00465 PF00467 PF00480 PF00483 PF00484 PF00491 PF00496 PF00497 PF00501 PF00502 PF00514 PF00515 PF00531 PF00532 PF00533 PF00534 PF00537 PF00541 PF00545 PF00550 PF00551 PF00554 PF00557 PF00560 PF00561 PF00571 PF00575 PF00576 PF00578 PF00579 PF00581 PF00583 PF00586 PF00587 PF00588 PF00590 PF00591 PF00595 PF00615 PF00621 PF00622 PF00626 PF00627 PF00628 PF00630 PF00644 PF00652 PF00653 PF00657 PF00670 PF00675 PF00680 PF00682 PF00685 PF00696 PF00698 PF00701 PF00702 PF00704 PF00705 PF00710 PF00719 PF00722 PF00724 PF00730 PF00740 PF00753 PF00754 PF00756 PF00775 PF00782 PF00787 PF00795 PF00806 PF00856 PF00857 PF00866 PF00881 PF00884 PF00885 PF00890 PF00891 PF00903 PF00905 PF00929 PF00936 PF00963 PF00969 PF00970 PF00989 PF00993 PF00994 PF01011 PF01023 PF01030 PF01037 PF01041 PF01042 PF01047 PF01048 PF01053 PF01063 PF01094 PF01118 PF01123 PF01138 PF01161 PF01168 PF01177 PF01180 PF01182 PF01187 PF01188 PF01193 PF01230 PF01239 PF01243 PF01248 PF01261 PF01263 PF01266 PF01315 PF01323 PF01327 PF01336 PF01344 PF01353 PF01370 PF01380 PF01381 PF01408 PF01419 PF01421 PF01423 PF01436 PF01451 PF01467 PF01477 PF01479 PF01487 PF01494 PF01497 PF01510 PF01546 PF01547 PF01564 PF01565 PF01568 PF01571 PF01588 PF01590 PF01593 PF01613 PF01627 PF01656 PF01661 PF01728 PF01791 PF01799 PF01833 PF01842 PF01850 PF01965 PF01966 PF01979 PF02085 PF02136 PF02210 PF02310 PF02332 PF02518 PF02525 PF02566 PF02597 PF02627 PF02629 PF02635 PF02729 PF02738 PF02746 PF02769 PF02770 PF02771 PF02774 PF02775 PF02776 PF02777 PF02779 PF02780 PF02782 PF02784 PF02788 PF02797 PF02798 PF02800 PF02801 PF02803 PF02806 PF02807 PF02812 PF02826 PF02832 PF02837 PF02852 PF02861 PF02866 PF02874 PF02875 PF02876 PF02878 PF02881 PF02885 PF02894 PF02907 PF03009 PF03061 PF03070 PF03099 PF03104 PF03129 PF03143 PF03144 PF03167 PF03358 PF03372 PF03422 PF03446 PF03466 PF03471 PF03496 PF03725 PF03727 PF03810 PF03952 PF03989 PF03992 PF04616 PF05191 PF05193 PF05199 PF05221 PF05368 PF05738 PF05932 PF06026 PF06628 PF06696 PF07650 PF07653 PF07654 PF07679 PF07686 PF07687 PF07702 PF07714 PF07715 PF07731 PF07732 PF07859 PF07883 PF07969 PF07980 PF07992 PF08240 PF08241 PF08245 PF08282 PF08447 PF08501 PF08534 PF08541 PF08544 PF08545 PF10535 PF10584 PF12680 PF12681 PF12695 PF12697 PF12796 PF12799 PF12840 PF12847"

if [[ "$step" -eq 1 ]]
then
    #################################################
    #	   Create MSA for large-scale datasets 	    #
    #################################################
    counter=0
    for i in $famNames
    {
        aln_dir=$out_dir/$i
        seq_dir=$largeScalePath

        if [ ! -d $aln_dir ]
        then
                mkdir -p $aln_dir
        fi

        ((counter++))
        if [[ ! -f $aln_dir/mafft.fa || ! -s $aln_dir/mafft.fa ]]
        then
               echo "#!/bin/bash" >xmaf$counter.sh
	       echo "unset MAFFT_BINARIES" >>xmaf$counter.sh
               echo "time mafft --anysymbol --parttree --thread $n_threads --quiet $seq_dir/$i.fasta > $aln_dir/mafft.fa" >> xmaf$counter.sh
               chmod u+x xmaf$counter.sh
               #qsub -cwd -V -o $counter.o -e $counter.e -q mem_6,$MY_NODES ./x$counter.sh
               qsub -cwd -V -l h_vmem=10G -pe smp $n_threads -o x$counter.o -j y ./xmaf$counter.sh
        fi

        ((counter++))
        if [[ ! -f $aln_dir/clustalo.fa || ! -s $aln_dir/clustalo.fa ]]
        then
               echo "#!/bin/bash" >xclu$counter.sh
               echo "time clustalo --threads $n_threads -i $seq_dir/$i.fasta -o $aln_dir/clustalo.fa" >> xclu$counter.sh
               chmod u+x xclu$counter.sh
               #qsub -cwd -V -o $counter.o -e $counter.e -q mem_6,$MY_NODES ./x$counter.sh
               qsub -cwd -V -l h_vmem=10G -pe smp $n_threads -o x$counter.o -j y ./xclu$counter.sh
        fi
     }
fi


if [[ "$step" -eq 2 ]]
then
    #################################################
    #		  Create sp_library 		    #
    #################################################
    if [ ! -d "splib" ]
    then
        mkdir -p splib
    fi
    
    for i in $famNames
    {
      outlib="splib/${i}_sp_lib"

      if [[ ! -f $outlib || ! -s $outlib ]]
      then
	t_coffee -lib /users/cn/cmagis/PROJET/PHYLO3D/DATASET/$i/PFAM25.0/sap.lib /users/cn/cmagis/PROJET/PHYLO3D/DATASET/$i/PFAM25.0/mustang.lib /users/cn/cmagis/PROJET/PHYLO3D/DATASET/$i/PFAM25.0/TMalign.lib -output sp_ascii -outfile $outlib
      fi
    }
fi


if [[ "$step" -eq 3 ]]
then
    #########################################################################################################
    #		  Extract the MSA of the seqs in the lib from the respective large scale MSA		    #
    #########################################################################################################
    
    if [ ! -d "fam_extractedMSA" ]
    then
        mkdir -p fam_extractedMSA
    fi

    for i in $famNames
    {
	# --- input files --- #
	lmsa_clustal="$out_dir/$i/clustalo.fa"
	lmsa_mafft="$out_dir/$i/mafft.fa"
	#lmsa_km_tc="$largeScalePath/km_tc_k200_gap300.fa"
	lib_file="splib/${i}_sp_lib"

	# --- output files --- #
	msa_clustal="fam_extractedMSA/${i}_clustalo.fa"
	msa_mafft="fam_extractedMSA/${i}_mafft.fa"
	#msa_km_tc="${i}_km_tc_k200_gap300.fa"
	errorLog="fam_extractedMSA/${i}_error.log"

	if [ ! -f $lib_file ]  		# check if lib_file exists else go to next one
	then
	    echo "The splib doesn't exist: $lib_file" >> $errorLog
	    continue
	fi
	if [ ! -d $out_dir/$i ] 	# check if large_scale folder exists else print error and go to next one
	then
	    echo "This folder doesn't exist: $aln_dir/" >> $errorLog
	    continue
	fi
   

	## clustalo ##
	if [ ! -f $lmsa_clustal ] 
	then
	    echo "The large scale MSA of clustalo doesn't exist: $lmsa_clustal" >> $errorLog
	    continue
	fi
	if [ ! -f $msa_clustal ] 
	then
	    perl extract_subAln.pl $lib_file $lmsa_clustal
	fi


	## mafft ##
	if [ ! -f $lmsa_mafft ] 
	then
	    echo "The large scale MSA of mafft doesn't exist: $lmsa_mafft" >> $errorLog
	    continue
	fi
	if [ ! -f $msa_mafft ] 
	then
	    perl extract_subAln.pl $lib_file $lmsa_mafft
	fi

	## km_coffee ##
<<COMMENT
	if [ ! -f $lmsa_km_tc ]  # check if large scale exist
	then
	    echo "The large scale MSA of km_coffee doesn't exist: $lmsa_km_tc" >> PF000${i}_error.log
	    continue
	fi
	    if [ ! -f $msa_km_tc ]   # check if you have created the small MSA file before
	then
	    perl extract_subAln.pl $lib_file $lmsa_km_tc
	fi
COMMENT
    }
fi


if [[ "$step" -eq 4 ]]
then
    ##################################################################
    #		Evaluate, using as reference your sp_library	     #
    ##################################################################

    if [ ! -d "evalRES" ]
    then
          mkdir -p evalRES
    fi

    for i in $famNames
    {
	# --- input files --- #
	lib_file="splib/${i}_sp_lib"
	msa_clustal="fam_extractedMSA/${i}_clustalo.fa"
	msa_mafft="fam_extractedMSA/${i}_mafft.fa"
	#msa_km_tc="${i}_km_tc_k200_gap300.fa"
	errorLog="fam_extractedMSA/${i}_error.log"

	# --- output files --- #
	evalRes_file="evalRES/${i}_evaluation.Res"

	if [ ! -s $errorLog ]
	then
	    #rm $errorLog

	    if [ -f $msa_clustal ]
	    then
		t_coffee -other_pg aln_compare -lib $lib_file -al2 $msa_clustal >> $evalRes_file
	    fi

	    if [ -f $msa_mafft ]
	    then
		t_coffee -other_pg aln_compare -lib $lib_file -al2 $msa_mafft >> $evalRes_file
	    fi
    
#	    if [ -f $msa_km_tc ]
#	    then
#	        t_coffee -other_pg aln_compare -lib $lib_file -al2 $msa_km_tc >> PF000${i}_evaluation.Res
#	    fi
	fi

#	t_coffee -other_pg aln_compare -lib PF00004_sp_lib -al2 /users/cn/cmagis/PROJET_PHYLO3D/ANALYSIS/PFAM_PDB/PF00004_set/PFAM25.0_9.03/mustang.aln
#	t_coffee -other_pg aln_compare -lib PF00004_sp_lib -al2 /users/cn/cmagis/PROJET_PHYLO3D/ANALYSIS/PFAM_PDB/PF00004_set/PFAM25.0_9.03/tmalign.aln
#	t_coffee -other_pg aln_compare -lib PF00004_sp_lib -al2 /users/cn/cmagis/PROJET_PHYLO3D/ANALYSIS/PFAM_PDB/PF00004_set/PFAM25.0_9.03/sap.aln
#	t_coffee -other_pg aln_compare -lib PF00004_sp_lib -al2 /users/cn/cmagis/PROJET_PHYLO3D/ANALYSIS/PFAM_PDB/PF00004_set/PFAM25.0_9.03/psicoffee.aln
#	t_coffee -other_pg aln_compare -lib PF00004_sp_lib -al2 /users/cn/cmagis/PROJET_PHYLO3D/ANALYSIS/PFAM_PDB/PF00004_set/PFAM25.0_9.03/msaprobs.aln 
#	t_coffee -other_pg aln_compare -lib PF00004_sp_lib -al2 /users/cn/cmagis/PROJET_PHYLO3D/ANALYSIS/PFAM_PDB/PF00004_set/PFAM25.0_9.03/probcons.aln
#	t_coffee -other_pg aln_compare -lib PF00004_sp_lib -al2 /users/cn/cmagis/PROJET_PHYLO3D/ANALYSIS/PFAM_PDB/PF00004_set/PFAM25.0_9.03/clustalw.aln
#	t_coffee -other_pg aln_compare -lib PF00004_sp_lib -al2 /users/cn/cmagis/PROJET_PHYLO3D/ANALYSIS/PFAM_PDB/PF00004_set/PFAM25.0_9.03/mustang.aln
#	t_coffee -other_pg aln_compare -lib PF00004_sp_lib -al2 /users/cn/cmagis/PROJET_PHYLO3D/ANALYSIS/PFAM_PDB/PF00004_set/PFAM25.0_9.03/mustang.aln

  
    }
fi