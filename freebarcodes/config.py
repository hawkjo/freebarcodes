import logging
import os


class CommandLineArguments(object):
    """
    Wraps the raw arguments provided by docopt.
    """
    def __init__(self, arguments, current_directory):
        self._arguments = arguments
        self._current_directory = current_directory

    def _comma_delimited_arg(self, key):
        if self._arguments[key]:
            return self._arguments[key].split(',')
        return None

    @property
    def barcode_files(self):
        return self._comma_delimited_arg('<barcode_files>')

    @property
    def command(self):
        # We have to do this weird loop to deal with the way docopt stores the command name
        for possible_command in ('decode',
                                 'generate',
                                 'prune',
                                 'concatenate'):
            if self._arguments.get(possible_command):
                return possible_command

    @property
    def fastq_files(self):
        return self._comma_delimited_arg('<fastq_files>')

    @property
    def log_level(self):
        log_level = {0: logging.ERROR,
                     1: logging.WARN,
                     2: logging.INFO,
                     3: logging.DEBUG}
        # default to silent if the user supplies no verbosity setting
        return log_level.get(self._arguments['-v'], logging.ERROR)

    @property
    def num_errors(self):
        return int(self._arguments['<num_errors>'] or 0)

    @property
    def output_dir(self):
        return self._arguments['<output_dir>'] or None

    @property
    def prefixes(self):
        return self._comma_delimited_arg('<prefixes>')

    @property
    def raw_barcodes_file(self):
        return self._arguments['<raw_barcodes_file>'] or None

