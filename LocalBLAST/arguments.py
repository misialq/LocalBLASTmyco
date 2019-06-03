from argparse import ArgumentParser

# Input arguments to the LocalMultiBLAST


def create_argument_parser():
    arg_parser = ArgumentParser(
        description='See below the list of possible arguments.')
    arg_parser.add_argument(
        '--verbosity',
        '-v',
        type=int,
        default=2,
        choices=range(1, 4),
        help='Verbosity level (1-3, default 2)')
    arg_parser.add_argument(
        '--use-uniprot',
        '-uni',
        action='store_true',
        help='If true, will fetch additional info from UniProt (requires network connection)')
    arg_parser.add_argument(
        '--clean-up',
        '-c',
        action='store_true',
        help='If true, will remove the folder containing temporary BLAST XML files')
    arg_parser.add_argument(
        '--sequences-loc',
        '-sloc',
        type=str,
        default='/sequences',
        help='Location of original sequences used for BLAST')
    arg_parser.add_argument(
        '--results-loc',
        '-rloc',
        type=str,
        default='/results',
        help='Location where the results will be saved')
    arg_parser.add_argument(
        '--blastdb',
        '-db',
        type=str,
        default='/blast/db/mycoRv_proteins_uniprot.fasta',
        help='Name of the BLAST database')
    arg_parser.add_argument(
        '--input-format',
        '-inp',
        type=str,
        default='fas',
        help='Format of the original sequences that will be submitted to BLAST (default is FAS)')
    arg_parser.add_argument(
        '--blast-engine',
        '-be',
        type=str,
        default='blastx',
        choices=('blastp', 'blastx', 'blastn'),
        help='BLAST engine to be employed')
    arg_parser.add_argument(
        '--sequence-len-threshold',
        '-slt',
        type=int,
        default=100,
        choices=range(50, 100),
        help='Minimum length of a sequence required for BLAST')
    arg_parser.add_argument(
        '--max-proc-count',
        '-p',
        type=int,
        default=2,
        choices=range(1, 20),
        help='Number of parallel processes to be used. Defaults to 2.')



    return arg_parser
