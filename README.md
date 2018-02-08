# How to build BENCHFAM !

BenchFam creates automatically a protein benchmark dataset that can be used to
evaluate the performance of alignment programs, even for large scale aligners.

This is done by getting all protein families from the PFAM database that belong
that contain at least 10 members with known crystallographic structures. This
treshold can be modify in an early stage when preparing the data to give to
the nextflow workflow. The quality of the stuctures depends on the results of
the BLAST against the PDB database (it may not always be the best).

## 1 - Installing NextFlow and benchfam

### NextFlow
Install Nextflow runtime with this command:

    curl -fsSL get.nextflow.io | bash

#### Prerequisites
 
* Java 7 or later
* Docker engine 1.0 or later (in alternative you can install the required
dependecies as shown in the included [Dockerfile](Dockerfile))


### BenchFam
BenchFam is on GitHub; clone it from the CBCRG git account:

    git clone https://github.com/cbcrg/benchfam.git

All PERL and bash scripts, binaries, etc...are in the $CWD/benchfam/bin/.


## 2 - Getting the data

### Download the PFAM database 
You can either download from the webpage directly or by using wget:

    wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam28.0/Pfam-A.full.gz

### Uncompress the PFAM database 
It is a big single text file (>100G), it takes time !!!

    gunzip Pfam-A.full.gz

### Extract sequences and alignment from PFAM

The PFAM database is a single Stockholm formatted text file. It is therefore
necessary to extract all the families beforehand. To do so, use the perl
script PFAM_extract.pl. It requires three argument:

- "FULL" or "PDB" (full extract all sequences from a family; PDB only 3D)
- number (number max of family to extract; anything except a number extract all)
- the Pfam database file (or any Pfam formatted file)

The script will extract in the current working directory all the families of 
PFAM aligned and unaligned; you have to run it twice (two terminals) to get all
sequences or only PDB ones:

    perl PFAM_extract.pl PDB X Pfam-A.full
    perl PFAM_extract.pl FULL X Pfam-A.full


CAUTION: this step takes ~1 day, because of the size of the database 

### Create a working directory for NextFlow workflow

You must create a folder containing families with a minimum of 10 sequences
that will be analyze by the nextflow pipeline. To do this use the following 
bash script:

    ./PFAM_count.sh

For the moment the script must be located where the files are, but this can be
changed. It is an interactive script that will ask you two parameters: 
- X: the minimum number of sequences per file(10)
- folder where to move the families with Nseq >= X sequences
- folder where to move the families with Nseq < X sequences

## 3 - Run BenchFam with NextFlow

There are still some improvment to do; for now, the path of the files extracted
from Pfam-A.full has to be harcoded inside the pipeline. First, you need to get
BenchFam from github and have NextFlow installed (step 1).

The BenchFam pipeline performs 3 operations:
- 1_filter_PDB
Performs a first layer of filtering of the number of sequences, %id, coverage,
etc...and identifies 3D template from the PDB.
- 2_extract_PDB
Performs an extract of the exact part of the structure file corresponding to the
sequences within the sequence file. It cuts out the 3D structure of each protein
and rename the file with an index according to the number of domains within the
same PDB file (it is important especially when domains are repeats).
- 3_align_libraries
Performs all alignments using a multitude of aligners based on structures or on
sequence (e.g. SAP, TMalign, T-Coffee...)

The command line to run benchfam.nf is:

    nextflow run benchfam.nf -c pfam28_2018.config

If the run is stopped for any reason, it can be resumed:

    nextflow run benchfam.nf -c pfam28_2018.config -resume

## 4 - Organize & Test BenchFam output

### Organize BenchFam from the nextflow log file

The BenchFam pipeline will generate many output files stored in a scratch folder
even when for some datasets fail at some step of the workflow. To get & organize
the data, we use a perl script : BENCHFAM_LOG_extract.pl. The script needs two
input arguments:
- the log file from BenchFam nextFlow
- the full PATH for the scratch directory

Generate the script to organize BenchFam:

    perl BENCHFAM_LOG_extract.pl trace_pfam_28.0_2018 $PATH/PFAM_28.0_scratch

The script generated before is a bash script run_copy_scratch.sh which has to be
run in the SCRATCH folder $PATH/PFAM_28.0_scratch in order to organize BENCHFAM:

    ./run_copy_scratch.sh


### Test the final results

Some steps of the pipeline are not full proof, espcially the PDB_extract.pl due
to discrepancies between the PDB information from BLAST and the information in
the corresponding PDB file (missing regions, mutations, etc...)
In order to verify that there is no problematic case, the best way is to check
in the EVAL folder the IRMSD files. In case of disgreement between PDB sequence
and PFAM sequence, the IRMSD package of T-Coffee will attribute a value of -1.0
to the comparisons involving these structure. 

    grep " -1.00" xxx.irmsd | wc -l

If the count is not NULL, then the dataset must be discarded or manually curated.




