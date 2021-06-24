"""
Free Divergence Error-Correcting Barcodes

Usage:
  freebarcodes decode       <barcode_files> <fastq_files> [--output-dir=<output_dir>] [--prefixes=<prefixes>] [--max-prefix-err=<max_prefix_err>] [-v | -vv | -vvv]
  freebarcodes generate     <barcode_length> <num_errors> [--output-dir=<output_dir>] [--cont=<prev_bc_fpath>] [--exclude=<exclude_bc_fpath>] [--4sets] [-v | -vv | -vvv]
  freebarcodes prune        <raw_barcodes_file> <num_errors> [--output-dir=<output_dir>] [-v | -vv | -vvv]
  freebarcodes concatenate  <barcode_files> [--output-dir=<output_dir>] [--max_bc=<max_bc>] [-v | -vv | -vvv]

Options:
  -h --help     Show this screen.
  --version     Show version.

Commands:
  decode        Decode barcodes in fastq files. Separate file names with commas.
  generate      Generate new sets of FREE barcodes of given length and number correctable errors
  prune         Prune a previous set of barcodes to a subset of valid FREE barcodes
  concatenate   Concatenate and filter multiple sets of barcodes

"""
import logging
import os
from docopt import docopt
from .__init__ import __version__
from .config import CommandLineArguments
from .decode import decode_fastqs
from .generate import generate_barcodes
from .prune import prune_barcodes
from .concatenate import concatenate_barcodes


def main(**kwargs):
    docopt_args = docopt(__doc__, version=__version__)
    arguments = CommandLineArguments(docopt_args, os.getcwd())

    log = logging.getLogger()
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s   %(message)s", "%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(arguments.log_level)
    log.debug(docopt_args)

    commands = {
        'decode': decode_fastqs,
        'generate': generate_barcodes,
        'prune': prune_barcodes,
        'concatenate': concatenate_barcodes
    }

    commands[arguments.command](arguments)


if __name__ == '__main__':
    main()
