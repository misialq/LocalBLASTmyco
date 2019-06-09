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

## Arguments

| Argument | Default value | Description |
|---|---|---|
| --verbosity, -v | 2 | Verbosity level (1-3) |
| --use-uniprot, -uni |  | If true, will fetch additional info from UniProt (requires network connection) |
| --clean-up, -c |  | If true, will remove the folder containing temporary BLAST XML files |
| --sequences-loc, -sloc | /sequences | Location of original sequences used for BLAST |
| --results-loc, -rloc | /results | Location where the results will be saved |
| --blastdb, -db | /blast/db/mycoRv_proteins_uniprot.fasta | Name of the BLAST database |
| --input-format, -inp | fas | Format of the original sequences that will be submitted to BLAST (default is FAS) |
| --blast-engine, -be | blastx | BLAST engine to be employed (blastx, blastp, blastn) |
| --sequence-len-threshold, -slt | 100 | Minimum length of a sequence required for BLAST (shorter sequences will be discarded) |
| --max-proc-count, -p | 2 | Number of parallel processes to be used |


## Usage example

We are given a list of sequencing files (_fas_ format; files contain genomic CDS fragments) for which we want to identify 
which proteins are encoded by respective CDS. We need to perform a BLASTX search (nucleotide-to-protein) against the 
local database of mycobacterial proteins. We prepared the database as described above using UniProt as a source of all 
proteins encoded in the mycobacterial genome (`myco_proteins_uniprot`).

1. Install and activate the pipenv environment by executing `pipenv shell` (for more information on how to use Pipenv see 
[this link](https://pipenv.readthedocs.io/en/latest/)).
2. Make sure you know the location of the BLAST database.
3. In the _local-blast-myco_ location execute `python LocalBLAST.py -sloc "~/test/input_sequences" -rloc "~/test/results" -slt 60 -db "~/test/db/myco_proteins_uniprot"
`.

## Docker

1. Install Docker (for platform-specific instructions see [Docker docs](https://docs.docker.com/))
2. Run: ` docker run -v ~/test/input_sequences:/sequences -v ~/test/results:/results misialq/local-blast-myco:v0.2.0 -slt 60` (adjust the mount points using your desired locations)

Required mount points:

| System (adjust accordingly) | Container | Description |
|---|---|---|
| ~/test/input_sequences | /sequences | Location of the input sequences |
| ~/test/results | /results | Location where the results should be saved |
| ~/test/db | /blast/db| Location of additional BLAST databases (only if you want to overwrite the existing ones) |

Databases available within the Docker image (under `/blast/db/`):

| DB name | Description |
|---|---|
| mycoCDC_proteins_uniprot.fasta | All _M. tuberculosis_ CDC1551 proteins |
| mycoRv_proteins_uniprot.fasta | All _M. tuberculosis_ H37Rv proteins |
| myco_proteins_uniprot.fasta | Combination of the two above |