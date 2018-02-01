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
 * Cedrik Magis <cedrik.1978@gmail.com>
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

// -- result directory 
params.outdir = 'results'
resultDir = checkResultDir(params.outdir)

// -- parameters for structure filtering, trimming and aligning
//params.min_pdb = 10
//params.id_max = 95
//params.id_min = 95
//params.cov_min = 95
//params.window=5
//params.min_length=0.75
//params.max_length=1.50
//params.gaps_max=0.05
//params.id_filter=0.90

// -- parameters for BLAST/EXPRESSO

params.blastDb = '/db/ncbi/201604/blast/db/pdbaa.fa'
db_blast = file(params.blastDb)
expresso_params = params.blastDb in ['NCBI','EBI'] ? "-blast=${db_blast}" :  "-blast=LOCAL -pdb_db=${db_blast}"


// -- log summary 

log.info "  B E N C H F A M   ~   version 2.1 "
log.info "===================================="
log.info "blastDb          : ${params.blastDb}"
log.info "expresso_params  : ${expresso_params}"
log.info "min_pdb          : ${params.min_pdb}"
log.info "id_min           : ${params.id_min}"
log.info "id_max           : ${params.id_max}"
log.info "cov_min          : ${params.cov_min}"
log.info "window           : ${params.window}"
log.info "min_length       : ${params.min_length}"
log.info "max_length       : ${params.max_length}"
log.info "gaps_max         : ${params.gaps_max}"
log.info "id_filter        : ${params.id_filter}"

// -- start of the workflow

//proteins = Channel.fromPath('/users/cn/projects/PFAM_29.0/PFAM_more10/*_pdb.fa')
proteins = Channel.fromPath('/users/cn/cmagis/HOME/01_SOFT/BENCHFAM-cbcrg/benchfam/demo/*_pdb.fa')

process '1_filter_PDB' {
    tag { fasta.name }
    
    input:
    file fasta from proteins

    output:
    set ( fam, 'data.fasta', 'data_pdb1.template_list', '*.pdb') into temp_struct

    script:
    fam = fasta.baseName.endsWith('_pdb') ? fasta.baseName.replace('_pdb','') : fasta.baseName

    """
    grep '>' $fasta | wc -l > count
    echo count
    """
    
    """
    t_coffee -other_pg seq_reformat -in $fasta -action +trim _seq_%%${params.id_max}_ > data.fasta
    t_coffee data.fasta -mode expresso -pdb_type d -pdb_min_sim ${params.id_min} -pdb_min_cov ${params.cov_min} -multi_core=${task.cpus} -cache \$PWD $expresso_params
    """

}


process '2_extract_PDB' {
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
 */


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


process '3_align_libraries' {
    tag { fam }

    input:
    set ( fam, 'modified.fasta', 'modified.template', '*' ) from modified_struct1

    output:
    set (fam, '*.aln') into aln_files mode flatten
    set (fam, '*_irmsd') into irmsd_files mode flatten
    set (fam, 'sap.lib:mustang.lib:tmalign.lib' ) into lib_files

    """
    unset MAFFT_BINARIES
    replace_U.pl modified.fasta

    cp modified.fasta 3Dmcoffee.fasta
    cp modified.fasta clustalw.fasta
    cp modified.fasta mafft.fasta
    cp modified.fasta mcoffee.fasta
    cp modified.fasta msaprobs.fasta
    cp modified.fasta muscle.fasta
    cp modified.fasta mustang.fasta
    cp modified.fasta mustang_proba.fasta
    cp modified.fasta mustang_tmalign.fasta
    cp modified.fasta prank.fasta
    cp modified.fasta probcons.fasta
    cp modified.fasta sap.fasta
    cp modified.fasta sap_mustang.fasta
    cp modified.fasta sap_proba.fasta
#    cp modified.fasta sate.fasta
    cp modified.fasta tcoffee.fasta
    cp modified.fasta tmalign.fasta
    cp modified.fasta tmalign_proba.fasta
    cp modified.fasta tmalign_sap.fasta

    # Create libraries to be combine after
    t_coffee sap.fasta -template_file modified.template -method sap_pair -out_lib sap.lib -multi_core=${task.cpus}
    t_coffee mustang.fasta -template_file modified.template -method mustang_pair -out_lib mustang.lib -multi_core=${task.cpus}
    t_coffee tmalign.fasta -template_file modified.template -method TMalign_pair -out_lib tmalign.lib -multi_core=${task.cpus}
    t_coffee tcoffee.fasta -out_lib tcoffee.lib -multi_core=${task.cpus}
    t_coffee mcoffee.fasta -mode mcoffee -out_lib mcoffee.lib -multi_core=${task.cpus}

    # Computes all MSAs which do not require PDB structures
    t_coffee 3Dmcoffee.fasta -lib sap.lib mustang.lib tmalign.lib -multi_core=${task.cpus}
    t_coffee sap_proba.fasta -lib tcoffee.lib sap.lib -multi_core=${task.cpus}
    t_coffee mustang_proba.fasta -lib tcoffee.lib mustang.lib -multi_core=${task.cpus}
    t_coffee tmalign_proba.fasta -lib tcoffee.lib tmalign.lib -multi_core=${task.cpus}
    t_coffee sap_mustang.fasta -lib sap.lib mustang.lib -multi_core=${task.cpus}
    t_coffee mustang_tmalign.fasta -lib mustang.lib tmalign.lib -multi_core=${task.cpus}
    t_coffee tmalign_sap.fasta -lib tmalign.lib sap.lib -multi_core=${task.cpus}

    clustalw clustalw.fasta
    mafft --quiet --thread ${task.cpus} mafft.fasta > mafft.temp
    t_coffee -other_pg seq_reformat mafft.temp -output clustalw > mafft.aln

    msaprobs msaprobs.fasta -o msaprobs.temp
    t_coffee -other_pg seq_reformat msaprobs.temp -output clustalw > msaprobs.aln

    muscle -in muscle.fasta -out muscle.temp
    t_coffee -other_pg seq_reformat muscle.temp -output clustalw > muscle.aln

    prank -d=prank.fasta -o=prank.temp
    t_coffee -other_pg seq_reformat prank.temp.2.fas -output clustalw > prank.aln

    probcons probcons.fasta > probcons.temp
    t_coffee -other_pg seq_reformat probcons.temp -output clustalw > probcons.aln

#    export PYTHONPATH="/users/cn/cmagis/bin/site-packages/"
#    python /users/cn/cmagis/bin/satesrc-v2.2.5-2012Oct16/sate-core/run_sate.py --input sate.fasta --datatype=Protein --num-cpus=${task.cpus} --output-directory=tmp_sate --auto
#    t_coffee -other_pg seq_reformat ./tmp_sate/satejob*.marker001.sate.aln -output clustalw > sate.aln

    # Computes the IRMSD-NIRMSD scores of all  MSAs
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
#    t_coffee -other_pg irmsd sate.aln -template_file modified.template > sate_irmsd

    """
}

/*
 * Save the alignment files
 */

aln_files.collectFile(storeDir: resultDir) { entry ->
    def fam = entry[0]
    def file = entry[1]
    [ "${fam}_${file.name}", file ]
}

/*
 * Save fasta and templates
 */
modified_struct2.subscribe { entry ->
    def fam = entry[0]
    def fasta = entry[1]
    def template = entry[2]
    fasta.copyTo( resultDir.resolve("${fam}_modified.fasta") )
    template.copyTo( resultDir.resolve("${fam}_template.fasta") )
}

/*
 * Save the irmsd files
 */

irmsd_files.collectFile(storeDir: resultDir) { entry ->
    def fam = entry[0]
    def file = entry[1]
    [ "${fam}_${file.name}", file ]
}

pml.subscribe { fam, file -> file.copyTo( resultDir / "${fam}_super.pml" ) }

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

