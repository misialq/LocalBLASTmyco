# LocalBLASTmyco

## Windows installation

1. Download the current version of standalone BLAST from NCBI [FTP server](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
2. Install in the default location
3. If not done automatically, specify BLAST database environment variable (name: BLASTDB, value: path to the database storage location)

## Database preparation

1. Prepare a FASTA file containing all the protein sequences to be contained in the database
2. Using the command line, go to the location where the database should be stored (same as specified in the env variable)
3. Also, copy the FASTA files into this location
4. Execute: `makeblastdb -in filename.fasta -dbtype prot` (can use option `-out` to specify a name for the output files)
