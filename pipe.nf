#!/usr/bin/env nextflow

params.in = 'tutorial/sample.fasta'
params.pdb_db='/db/pdb/derived_data_format/blast/latest/pdb_seqres.fa'
params.pdb_dir='/db/pdb/data/structures/divided/pdb'

structuralMethods = 'sap mustang TMalign'
normalMethods = 'tcoffee mcoffee psicoffee 3Dmcoffee sap_proba mustang_proba tmalign_proba clustalw mafft msaprobs muscle prank probcons sate'.split()


fasta = channel( file(params.in))

structural_fasta = list()
structural_template = list()

task {
  input fasta
  output 'modified.fasta': structural_fasta 
  output 'modified.template': structural_template
  echo true

  """ 
  export PDB_DIR=${params.pdb_dir}

  t_coffee $fasta  -mode expresso -blast=LOCAL -pdb_db=${params.pdb_db} -pdb_type d -pdb_min_sim 95 -pdb_min_cov 95 -cache \$PWD
  grep _P_  *_pdb1.template_list > temp.list
  t_coffee -other_pg seq_reformat -in $fasta -action +extract_seq_list temp.list > temp.fasta

  PDB_extract.pl
  
  """
}


task {
  input structural_fasta 
  input structural_template
  echo true

  """

  ln -s ${structural_fasta} sap.fasta
  ln -s ${structural_fasta} mustang.fasta
  ln -s ${structural_fasta} tmalign.fasta
  ln -s ${structural_fasta} tcoffee.fasta
  ln -s ${structural_fasta} mcoffee.fasta
  ln -s ${structural_fasta} psicoffee.fasta
  ln -s ${structural_fasta} 3Dmcoffee.fasta
  ln -s ${structural_fasta} sap_proba.fasta
  ln -s ${structural_fasta} mustang_proba.fasta
  ln -s ${structural_fasta} tmalign_proba.fasta
  ln -s ${structural_fasta} clustalw.fasta
  ln -s ${structural_fasta} mafft.fasta
  ln -s ${structural_fasta} msaprobs.fasta
  ln -s ${structural_fasta} muscle.fasta
  ln -s ${structural_fasta} prank.fasta
  ln -s ${structural_fasta} probcons.fasta
  ln -s ${structural_fasta} sate.fasta

  # ALIGNMENT 
  t_coffee sap.fasta -template_file ${structural_template -method sap_pair -out_lib sap.lib
  t_coffee mustang.fasta -template_file ${structural_template -method mustang_pair -out_lib mustang.lib
  t_coffee tmalign.fasta -template_file structural_template -method TMalign_pair -out_lib tmalign.lib

  t_coffee tcoffee.fasta -out_lib tcoffee.lib
  t_coffee mcoffee.fasta -mode mcoffee -out_lib mcoffee.lib

  t_coffee psicoffee.fasta -mode psicoffee -blast=LOCAL -protein_db=/db/ncbi/201304/blast/db/nr.fa -out_lib psicoffee.lib multi_core=no


  t_coffee sap_proba.fasta -lib tcoffee.lib sap.lib
  t_coffee mustang_proba.fasta -lib tcoffee.lib mustang.lib
  t_coffee tmalign_proba.fasta -lib tcoffee.lib tmalign.lib

  t_coffee 3Dmcoffee.fasta -lib sap.lib mustang.lib tmalign.lib


  clustalw2 clustalw.fasta 
  mafft mafft.fasta > mafft.temp
  t_coffee -other_pg seq_reformat mafft.temp -output clustalw > mafft.aln
  msaprobs msaprobs.fasta -o msaprobs.temp
  t_coffee -other_pg seq_reformat msaprobs.temp -output clustalw > msaprobs.aln
  muscle -in muscle.fasta -out muscle.temp
  t_coffee -other_pg seq_reformat muscle.temp -output clustalw > muscle.aln
  prank -d=prank.fasta -o=prank.temp
  t_coffee -other_pg seq_reformat prank.temp.2.fas -output clustalw > prank.aln
  probcons probcons.fasta > probcons.temp
  t_coffee -other_pg seq_reformat probcons.temp -output clustalw > probcons.aln
  
  export PYTHONPATH="/users/cn/jchang/local/Python/lib/python2.6/site-packages/"
  python /users/cn/jchang/program/satesrc-v2.2.5-2012Oct16/sate-core/run_sate.py --input sate.fasta --datatype=Protein --num-cpus=1 --output-directory=tmp_sate --auto
  t_coffee -other_pg seq_reformat ./tmp_sate/satejob*.marker001.sate.aln -output clustalw > sate.aln

  """  
}
