"""
Free Divergence Error-Correcting Barcodes

Usage:
  freebarcodes.py decode       <barcode_files> <fastq_files> [--output-dir=<output_dir>] [--prefixes=<prefixes>] [-v | -vv | -vvv]
  freebarcodes.py generate     <barcode_length> <num_errors> [--output-dir=<output_dir>] [-v | -vv | -vvv]
  freebarcodes.py prune        <raw_barcodes_file> <num_errors> [--output-dir=<output_dir>] [-v | -vv | -vvv]
  freebarcodes.py concatenate  <barcode_files> [--output-dir=<output_dir>] [-v | -vv | -vvv]

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
#from freebarcodes.generate import generate, concatenate
#from freebarcodes.decode import decode_fastqs
from freebarcodes.prune import prune_barcode_file
from docopt import docopt


def main(**kwargs):
    docopt_args = docopt(__doc__, version=VERSION)
    arguments = CommandLineArguments(docopt_args, os.getcwd())
    return

    log = logging.getLogger()
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s   %(message)s", "%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(arguments.log_level)
    log.debug(docopt_args)

    commands = {'decode': decode,
                'generate': generate,
                'prune': prune,
                'concatenate': concatenate}

    commands[arguments.command].main(arguments)


if __name__ == '__main__':
    main()
