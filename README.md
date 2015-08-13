BENCHFAM
=========

BenchFam creates in an automatic way a benchmark dataset for evaluating the performance of alignment programs. 

This is done by getting all sequences from PFAM database that belong to families that have more than 10 members with known high quality structures. The stuctures as well as their quality are determined by quering PDB database.

Quick start 
-----------

Install Nextflow runtime with this command: 

    curl -fsSL get.nextflow.io | bash

Launch the pipeline executed using the embedded demo dataset entering the following command in your shell terminal: 

    nextflow run cbcrg/benchfam -profile demo -with-docker 


Prerequisites 
-------------

* Java 7 or later 
* Docker engine 1.0 or later (in alternative you can install the required dependecies as shown in the included [Dockerfile](Dockerfile))
 
