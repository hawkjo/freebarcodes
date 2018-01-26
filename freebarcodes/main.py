"""
Free Divergence Error-Correcting Barcodes

Usage:
  freebarcodes decode       <barcode_files> <fastq_files> [--output-dir=<output_dir>] [--prefixes=<prefixes>] [-v | -vv | -vvv]
  freebarcodes generate     <barcode_length> <num_errors> [--output-dir=<output_dir>] [-v | -vv | -vvv]
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
from freebarcodes.constants import VERSION
from freebarcodes.config import CommandLineArguments
#from freebarcodes.decode import decode_fastqs
from freebarcodes.generate import generate_barcodes 
from freebarcodes.prune import prune_barcodes
from freebarcodes.concatenate import concatenate_barcodes
from docopt import docopt


def main(**kwargs):
    docopt_args = docopt(__doc__, version=VERSION)
    arguments = CommandLineArguments(docopt_args, os.getcwd())

    log = logging.getLogger()
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s   %(message)s", "%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(arguments.log_level)
    log.debug(docopt_args)

    commands = {
#        'decode': decode_fastqs,
        'generate': generate_barcodes,
        'prune': prune_barcodes,
        'concatenate': concatenate_barcodes
    }

    commands[arguments.command](arguments)


if __name__ == '__main__':
    main()
