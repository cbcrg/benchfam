#!/usr/bin/env nextflow

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
  output '*.lib'
  echo true

  """
  for x in ${structuralMethods}; do 
    t_coffee ${structural_fasta} -template_file ${structural_template} -method \${x}_pair -out_lib \${x}.lib
  done
  """  
}
