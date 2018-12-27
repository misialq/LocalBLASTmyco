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

## Usage example

We are given a list of sequencing files (_fas_ format; files contain genomic CDS fragments) for which we want to identify 
which proteins are encoded by respective CDS'. We need to perform a BLASTX search (nucleotide-to-protein) against the 
local database of mycobacterial proteins. We prepared the database as described above using UniProt as a source of all 
proteins encoded in the mycobacterial genome (`myco_proteins_uniprot`).

1. Example config file required for this test can be found in: `config/LocalBLASTtest.cfg`.
2. Install and activate the pipenv environment by executing `pipenv shell` (for more information on how to use Pipenv see 
[this link](https://pipenv.readthedocs.io/en/latest/)).
3. Make sure that the BLAST database is located in the folder specified in the config file.
3. In the _local-blast-myco_ location execute `python LocalBLAST.py -cf config/LocalBLASTtest.cfg`.