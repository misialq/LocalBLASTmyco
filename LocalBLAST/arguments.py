from argparse import ArgumentParser

# Input arguments to the LocalMultiBLAST


def create_argument_parser():
    arg_parser = ArgumentParser(
        description='See below the list of possible arguments.')
    arg_parser.add_argument(
        '-cf',
        '--config-file',
        type=str,
        help='Path to the config file')
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
        '--temp-clean-up',
        '-clean',
        action='store_true',
        help='If true, will remove the folder containing temporary BLAST XML files')

    return arg_parser
