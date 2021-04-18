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
        return []

    @property
    def barcode_files(self):
        return [os.path.expanduser(fp) for fp in self._comma_delimited_arg('<barcode_files>')]

    @property
    def barcode_length(self):
        return int(self._arguments['<barcode_length>'] or 0)

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
        return [os.path.expanduser(fp) for fp in self._comma_delimited_arg('<fastq_files>')]

    @property
    def log_level(self):
        log_level = {0: logging.ERROR,
                     1: logging.WARN,
                     2: logging.INFO,
                     3: logging.DEBUG}
        # default to silent if the user supplies no verbosity setting
        return log_level.get(self._arguments['-v'], logging.ERROR)

    @property
    def max_bc(self):
        return int(self._arguments['--max_bc'] or 0)

    @property
    def max_prefix_err(self):
        return [int(max_err) for max_err in self._comma_delimited_arg('--max-prefix-err')]
    @property
    def num_errors(self):
        return int(self._arguments['<num_errors>'] or 0)

    @property
    def output_dir(self):
        return self._arguments['--output-dir'] or '.'

    @property
    def prefixes(self):
        return self._comma_delimited_arg('--prefixes')

    @property
    def raw_barcodes_file(self):
        return os.path.expanduser(self._arguments['<raw_barcodes_file>']) or None

    @property
    def generate_as_4sets(self):
        return self._arguments['--4sets']

    @property
    def prev_bc_fpath(self):
        return self._arguments['--cont']

    @property
    def exclude_bc_fpath(self):
        return self._arguments['--exclude']
