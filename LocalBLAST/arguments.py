from argparse import ArgumentParser

# Input arguments to the LocalMultiBLAST


def create_argument_parser():
    arg_parser = ArgumentParser(
        description='Fill out this description.')
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

    return arg_parser