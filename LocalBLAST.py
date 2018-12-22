import sys

from LocalBLAST.arguments import create_argument_parser
from LocalBLAST.blast import LocalMultiBLAST
from LocalBLAST.helpers import set_logger


def run(args):
    logger = set_logger(args.verbosity, log_to_stdout=True)

    local_blast_myco = LocalMultiBLAST(args.config_file)
    local_blast_myco.run()


if __name__ == '__main__':
    arg_parser = create_argument_parser()
    args = arg_parser.parse_args(sys.argv[1:])
    run(args)
