#!/usr/env nextflow

params.blastDb = "/db/pdb/derived_data_format/blast/latest/pdb_seqres.fa"
params.cpus = 1
params.pfamFullGz = '/db/pfam/latest/Pfam-A.full.gz'
params.limit = '50'
params.dbCache = 'db'

all_methods = ['mafft','clustalo']


// -- create a couple of channel
// 'dataset' emits all the PFAM fasta file to process
// 'famNames' emits the unique PFAM familiy names


db_pdb = file(params.dbCache).resolve('pdb') 
db_full = file(params.dbCache).resolve('full')
expresso_params = "-blast=LOCAL -pdb_db=${params.blastDb}" 

/* 
 * Uncompress the PFAM database extracing only sequences with structures
 */
process extractPdb {
  storeDir db_pdb  

  output: 
  file '*_pdb.fasta' into pdb_files mode flatten

  """
  gzip -c -d ${params.pfamFullGz} | PFAM_extract_full.pl PDB ${params.limit} -
  for x in *.fasta; do [ `grep '>' \$x -c` -lt 10 ] && rm \$x; done
  """
}

/* 
 * Uncompress the PFAM database extracing ALL sequences
 */
process extractFull {
  storeDir db_full

  output: 
  file '*_full.fasta' into full_files mode flatten

  """
  gzip -c -d ${params.pfamFullGz} | PFAM_extract_full.pl FULL ${params.limit} -
  """
}



/* 
 * receive in input the PFXXXX_pdb.fasta
 */
process filter {

    input:
    file fasta from pdb_files

    output:
    set ( fam, 'temp.list', 'temp.fasta', '*.pdb') into temp_struct

    script:
    fam = fasta.baseName.endsWith('_pdb') ? fasta.baseName.replace('_pdb','') : fasta.baseName   

    """
    mkdir OUTPUT
    t_coffee -other_pg seq_reformat -in $fasta -action +trim _seq_%%99_ > data_99.fasta
    t_coffee data_99.fasta -mode expresso -pdb_type d -pdb_min_sim 95 -pdb_min_cov 95 -multi_core=${params.cpus} -cache \$PWD $expresso_params
    mv *tmp *results *dnd *html plot.rms OUTPUT
    grep _P_  *_pdb1.template_list > temp.list
    t_coffee -other_pg seq_reformat -in data_99.fasta -action +extract_seq_list temp.list > temp.fasta
    """
}


process pdb_extract {
    input:
    set ( fam, 'temp.list','temp.fasta','*') from temp_struct

    output:
    set ( fam, 'modified.template','modified.fasta', '*-1.pdb' ) into modified_struct
	set ( fam, 'modified.fasta' ) into seq3d
    """
    PDB_extract.pl
    """

}

process Lib_and_Aln {

    input:
    set ( fam, 'modified.template', 'modified.fasta', '*' ) from modified_struct

    output:
    file '*_irmsd' into irmsd_files
    set (fam, 'sap.lib:mustang.lib:tmalign.lib' ) into lib_files

    """
    cp modified.fasta sap.fasta
    cp modified.fasta mustang.fasta
    cp modified.fasta tmalign.fasta
    cp modified.fasta tcoffee.fasta
    cp modified.fasta psicoffee.fasta
    cp modified.fasta mcoffee.fasta
    cp modified.fasta 3Dmcoffee.fasta
    cp modified.fasta sap_proba.fasta
    cp modified.fasta mustang_proba.fasta
    cp modified.fasta tmalign_proba.fasta
    cp modified.fasta sap_mustang.fasta
    cp modified.fasta mustang_tmalign.fasta
    cp modified.fasta tmalign_sap.fasta
    cp modified.fasta clustalw.fasta
    cp modified.fasta mafft.fasta
    cp modified.fasta msaprobs.fasta
    cp modified.fasta muscle.fasta
    cp modified.fasta prank.fasta
    cp modified.fasta probcons.fasta
    cp modified.fasta sate.fasta

    # Create libraries by combining other methods
    t_coffee sap.fasta -template_file modified.template -method sap_pair -out_lib sap.lib -multi_core=${params.cpus}
    t_coffee mustang.fasta -template_file modified.template -method mustang_pair -out_lib mustang.lib -multi_core=${params.cpus}
    t_coffee tmalign.fasta -template_file modified.template -method TMalign_pair -out_lib tmalign.lib -multi_core=${params.cpus}
    t_coffee tcoffee.fasta -out_lib tcoffee.lib -multi_core=${params.cpus}
    t_coffee mcoffee.fasta -mode mcoffee -out_lib mcoffee.lib -multi_core=${params.cpus}

    # This doesn't need the PDBs
    t_coffee 3Dmcoffee.fasta -lib sap.lib mustang.lib tmalign.lib -multi_core=${params.cpus}
    t_coffee sap_proba.fasta -lib tcoffee.lib sap.lib -multi_core=${params.cpus}
    t_coffee mustang_proba.fasta -lib tcoffee.lib mustang.lib -multi_core=${params.cpus}
    t_coffee tmalign_proba.fasta -lib tcoffee.lib tmalign.lib -multi_core=${params.cpus}
    t_coffee sap_mustang.fasta -lib sap.lib mustang.lib -multi_core=${params.cpus}
    t_coffee mustang_tmalign.fasta -lib mustang.lib tmalign.lib -multi_core=${params.cpus}
    t_coffee tmalign_sap.fasta -lib tmalign.lib sap.lib -multi_core=${params.cpus}

    clustalw clustalw.fasta
    mafft --thread ${params.cpus} mafft.fasta > mafft.temp
    t_coffee -other_pg seq_reformat mafft.temp -output clustalw > mafft.aln

    msaprobs msaprobs.fasta -o msaprobs.temp
    t_coffee -other_pg seq_reformat msaprobs.temp -output clustalw > msaprobs.aln

    muscle -in muscle.fasta -out muscle.temp
    t_coffee -other_pg seq_reformat muscle.temp -output clustalw > muscle.aln

    prank -d=prank.fasta -o=prank.temp
    t_coffee -other_pg seq_reformat prank.temp.2.fas -output clustalw > prank.aln

    probcons probcons.fasta > probcons.temp
    t_coffee -other_pg seq_reformat probcons.temp -output clustalw > probcons.aln

    python \$SATE_HOME/sate-core/run_sate.py --input sate.fasta --datatype=Protein --num-cpus=${params.cpus} --output-directory=tmp_sate --auto
    t_coffee -other_pg seq_reformat ./tmp_sate/satejob*.marker001.sate.aln -output clustalw > sate.aln

    # IRMSD-NIRMSD OF ALL MSAS
    t_coffee -other_pg irmsd sap.aln -template_file modified.template > sap_irmsd
    t_coffee -other_pg irmsd mustang.aln -template_file modified.template > mustang_irmsd
    t_coffee -other_pg irmsd tmalign.aln -template_file modified.template > tmalign_irmsd
    t_coffee -other_pg irmsd tcoffee.aln -template_file modified.template > tcoffee_irmsd
    t_coffee -other_pg irmsd mcoffee.aln -template_file modified.template > mcoffee_irmsd
    t_coffee -other_pg irmsd 3Dmcoffee.aln -template_file modified.template > 3Dmcoffee_irmsd
    t_coffee -other_pg irmsd sap_proba.aln -template_file modified.template > sap_proba_irmsd
    t_coffee -other_pg irmsd mustang_proba.aln -template_file modified.template > mustang_proba_irmsd
    t_coffee -other_pg irmsd tmalign_proba.aln -template_file modified.template > tmalign_proba_irmsd
    t_coffee -other_pg irmsd sap_mustang.aln -template_file modified.template > sap_mustang_irmsd
    t_coffee -other_pg irmsd mustang_tmalign.aln -template_file modified.template > mustang_tmalign_irmsd
    t_coffee -other_pg irmsd tmalign_sap.aln -template_file modified.template > tmalign_sap_irmsd
    t_coffee -other_pg irmsd clustalw.aln -template_file modified.template > clustalw_irmsd
    t_coffee -other_pg irmsd mafft.aln -template_file modified.template > mafft_irmsd
    t_coffee -other_pg irmsd msaprobs.aln -template_file modified.template > msaprobs_irmsd
    t_coffee -other_pg irmsd muscle.aln -template_file modified.template > muscle_irmsd
    t_coffee -other_pg irmsd prank.aln -template_file modified.template > prank_irmsd
    t_coffee -other_pg irmsd probcons.aln -template_file modified.template > probcons_irmsd
    t_coffee -other_pg irmsd sate.aln -template_file modified.template > sate_irmsd
    """
}

/* 
 * - Discards all the fasta files having less than 10 sequences 
 * - Collects all the family names for which there are at least 10 sequences and
 *   sends these names over the channel 'fam_names' 
 * - Sends tuple ( family name, fasta file ) over the channel 'fam_full'
 */

fam_full = Channel.create()
fam_names = Channel.create()

seq3d
	.filter { tuple -> 
		def file = tuple[1]
		int count = 0
		file.chopFasta { count++ } 
		return count >= 10 
	}
	
	.map { fam, file -> fam }
	
	.separate( fam_names, fam_full ) {  fam ->
		def fasta = db_full.resolve("${fam}_full.fasta")
		if( !fasta.exists() ) 
			log.warm "Missing file: $fasta"
		[  fam, [fam, fasta] ] 
	}


/*
 * Apply a MSA step
 * 
 * it received in input the PFXXXX_full.fasta
 */
 
process Large_scale_MSAs {

    input:
    set (fam, file(sequences)) from fam_full
    each method from all_methods

    output:
    set (fam, method, '*.aln') into large_msa

    script:
    alnName = "${fam}_${method}.aln"
    if( method=='mafft')
        """
        mafft --anysymbol --parttree --thread ${params.cpus} --quiet ${sequences} > $alnName
        """

    else if( method=='clustalo' )
        """
        clustalo --threads ${params.cpus} -i ${sequences} -o $alnName
        """

    else
        error "Unknown align method: $method"

}

fam_lib = fam_names
        .phase( lib_files )
        .map { fam, lib ->  lib }

process splib {
    input:
    set ( fam, '*' ) from fam_lib

    output:
    set (fam, '*.sp_lib') into sp_lib

    """
    t_coffee -lib sap.lib mustang.lib tmalign.lib -output sp_lib -outfile ${fam}.sp_lib
    """
}


/* 
 * split the channel in two to handle them separately
 */
 
(sp_lib1, sp_lib2) = sp_lib.split(2)

/* 
 * - Join each lib1 with the large msa for the corresponding family name 
 * - Create a channel named 'lib_and_msa' that will emit tuples like ( familyName, align method, sp_lib file, alignment file ) 
 */ 
lib_and_msa = sp_lib1
				.cross(large_msa)
				.map { lib, aln -> [ lib[0], aln[1], lib[1], aln[2] ] }    


process Extracted_msa {

    input:
    set fam, method, file(splib), file(aln) from lib_and_msa

    output:
    set fam, method, '*.extracted_msa' into extracted_msa


    """
    extract_subAln.pl \$PWD/${splib} \$PWD/${aln}

    #if [ -s ${fam}_error.log ]; then
    # echo There are erros in the log file. Check ${fam}_error.log
    # exit 1
    #fi
    mv ${fam}_${aln.baseName}.fa ${fam}_${aln.baseName}.extracted_msa
    """
}


msa_eval = sp_lib2
        .cross(extracted_msa)
        .map { lib,aln -> [ lib[0], aln[1], lib[1], aln[2] ] }  //  ( familyName, method, sp_lib file, alignment file )

process evaluate {

    input:
    set family, method, file(splib), file(msa) from msa_eval

    output:
    set family, method, '*.Res' into evaluation

    """
   t_coffee -other_pg aln_compare -lib ${splib} -al2 ${msa} >> ${family}_evalution.Res
   """

}

evaluation.subscribe { tuple ->

    println "Family: ${tuple[0]}; Method: ${tuple[1]}; score: ${tuple[2].text.trim().readLines()[1].split(/\s+/)[3]}"
    println '\n\n'

}

