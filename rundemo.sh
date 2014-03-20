CFG=''; if [ -f cluster.config ]; then CFG='-c cluster.config'; fi
./nextflow $CFG pfam3d.nf --limit 50 --blast-db tutorial/db_50/blast/pdb_small --db-cache tutorial/db_50/ "$@"

