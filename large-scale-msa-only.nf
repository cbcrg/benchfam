#!/usr/env nextflow

/*
 * Copyright (c) 2016, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'BENCHFAM'.
 *
 *   BENCHFAM is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   BENCHFAM is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with BENCHFAM.  If not, see <http://www.gnu.org/licenses/>.
 */


/* 
 * Main pipeline script 
 * 
 * @authors 
 * Cedrik Magis <cedrik.1978@gmail.com>ma
 * Maria Chatzou <mxatzou@gmail.com>
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

params.limit = 'all'
params.blastDb = "/db/pdb/derived_data_format/blast/latest/pdb_seqres.fa"
params.pfamFullGz = '/db/pfam/latest/Pfam-A.full.gz'
params.dbCache = "db_${params.limit}"
params.methods = 'mafft,clustalo,pasta,upp'
params.outdir = 'results'
params.alndir = 'alignments'

params.min_pdb = 10
params.id_max = 95
params.id_min = 95
params.cov_min = 95
params.window=5
params.min_length=0.75
params.max_length=1.50
params.gaps_max=0.05
params.id_filter=0.90

// --validate result directory
resultDir = checkResultDir(params.outdir)

// -- given a comma separated list of methods converts it to a list object 
all_methods = params.methods.split(',').collect { it.trim() }

// -- local paths where are stored sequence files extracted by the Pfam database 
params.db_pdb = "${params.dbCache}/pdb"
params.db_full = "${params.dbCache}/full"
params.pfam_aln = "${params.dbCache}/pfam/*"

// -- the LOCAL BAST database required by T-Coffee
db_blast = file(params.blastDb)
expresso_params = params.blastDb in ['NCBI','EBI'] ? "-blast=${db_blast}" :  "-blast=LOCAL -pdb_db=${db_blast}"

db_pdb = file(params.db_pdb)
db_full = file(params.db_full)

// -- summary 

log.info "B E N C H - F A M     ~   v. 1.4.2"
log.info "=================================="
log.info "blastDb           : ${params.blastDb}"
log.info "pfamFullGz        : ${params.pfamFullGz}"
log.info "dbCache           : ${params.dbCache}"
log.info "db_pdb            : ${params.db_pdb}"
log.info "db_full           : ${params.db_full}"
log.info "pfam_aln          : ${params.pfam_aln}"
log.info "limit             : ${params.limit}"
log.info "methods           : ${params.methods}"
log.info "expresso_params   : ${expresso_params}"
log.info "min_pdb           : ${params.min_pdb}"
log.info "id_min            : ${params.id_min}"
log.info "id_max            : ${params.id_max}"
log.info "cov_min           : ${params.cov_min}"
log.info "window            : ${params.window}"
log.info "min_length        : ${params.min_length}"
log.info "max_length        : ${params.max_length}"
log.info "gaps_max          : ${params.gaps_max}"
log.info "id_filter         : ${params.id_filter}"

/* 
 * Uncompress the PFAM database extracing only sequences with structures
 */
process '1_extractPdb' {
  storeDir db_pdb  

  output: 
  file '*_pdb.fa' into pdb_files mode flatten

  """
  gzip -c -d ${params.pfamFullGz} | PFAM_extract_full.pl PDB ${params.limit} -
  for x in *.fa; do [ `grep '>' \$x -c` -lt ${params.min_pdb} ] && rm \$x; done
  """
}

/* 
 * Uncompress the PFAM database extracting ALL sequences
 */
process '2_extractFull' {
  storeDir db_full

  output: 
  file '*_full.fa' into full_files mode flatten

  """
  gzip -c -d ${params.pfamFullGz} | PFAM_extract_full.pl FULL ${params.limit} -
  """
}


full_files
    .map { file -> [ file.baseName.replace('_full',''), file ] }
    .set { full_files2 }


/* 
 * receive in input the PFXXXX_pdb.fasta
 */
process '3_filter' {
    tag { fasta.name }

    input:
    file fasta from pdb_files

    output:
    set ( fam, 'data.fasta', 'data_pdb1.template_list', '*.pdb') into temp_struct

    script:
    fam = fasta.baseName.endsWith('_pdb') ? fasta.baseName.replace('_pdb','') : fasta.baseName   

    """
    t_coffee -other_pg seq_reformat -in $fasta -action +trim _seq_%%${params.id_max}_ > data.fasta
    t_coffee data.fasta -mode expresso -pdb_type d -pdb_min_sim ${params.id_min} -pdb_min_cov ${params.cov_min} -multi_core=${task.cpus} -cache \$PWD $expresso_params
    """
}


process '4_pdb_extract' {
    tag { fam }
    errorStrategy 'ignore'

    input:
    set ( fam, 'data.fasta', 'data_pdb1.template_list', '*') from temp_struct

    output:
    set ( fam, 'modified.fasta', 'modified.template', '*-*.pdb' ) into modified_struct
    set ( fam, 'super.pml' ) into pml

    """
    PDB_extract.pl data.fasta data_pdb1.template_list ${params.window} ${params.min_length} ${params.max_length} ${params.gaps_max} ${params.id_filter}
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

modified_struct.filter { tuple ->
            def count = tuple[1].countFasta()
            def valid = count >= params.min_pdb
            if( !valid )
                log.info "Discarding family: ${tuple[0]} because 'PDB_extract' returns less than ${params.min_pdb} structures ($count)"
            return valid
        }
        .tap{ modified_struct1 }
        .tap{ modified_struct2 }
        .map { tuple -> tuple[0] }
        .phase( full_files2 )
        .map { f, t ->  [ f, t ]  }
        .set { fam_full }



/*
 * Apply a MSA step
 * 
 * it received in input the PFXXXX_full.fasta
 */
 
process '6_Large_scale_MSAs' {
    tag { "$fam-$method" }
    errorStrategy 'ignore'
    publishDir params.alndir 

    input:
    set (fam, file(sequences)) from fam_full
    each method from all_methods

    output:
    set (fam, method, '*.aln') 

    script:
    alnName = "${fam}_${method}.aln"
    if( method=='mafft')
        """
        unset MAFFT_BINARIES 
        replace_U.pl ${sequences}
        mafft --quiet --anysymbol --parttree --thread ${task.cpus} --quiet ${sequences} > $alnName
        """

    else if( method=='clustalo' )
        """
        clustalo --threads ${task.cpus} -i ${sequences} -o $alnName 
        """
        
    else if( method == 'pasta' ) 
        """
        replace_U.pl $sequences 
	run_pasta.py --num-cpus ${task.cpus} -i $sequences -d Protein -j $sequences -o out
        mv out/${sequences}.marker* $alnName 
        """
  
    else if( method == 'upp' ) 
        """
        replace_U.pl $sequences 
	run_upp.py -s $sequences -m amino --cpu ${task.cpus} -d outdir  -o $alnName
        mv outdir/${alnName}_alignment.fasta $alnName
 	"""
    else if( method == 'mega' ) 
	"""
        replace_U.pl $sequences 
        mega_coffee -i $sequences --cluster_size 2 --cluster_number 5000 -n ${task.cpus} -o $alnName
	"""
    else
        error "Unknown align method: $method"

}



/*
 * Verify that the result dir is empty or create it if do not exist
 */
def checkResultDir( String path ) {
    def result = file(path)
    if( result.exists() && result.isDirectory() && result.isEmpty() )
        return result

    if( result.exists() && !result.isDirectory())
        exit 1, "The specified result path is a file: $result -- please delete it or provide a different result path"

    if( result.exists()  && !result.isEmpty() )
        exit 2, "The specified result path is not empty: $result -- please delete the content or provide a different path"

    if( !result.exists() && !result.mkdirs() )
        exit 3, "Unable to create the result folder: $result -- please write permissions or provide a different path"

    return result
}





