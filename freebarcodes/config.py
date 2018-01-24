import logging
import os


class CommandLineArguments(object):
    """
    Wraps the raw arguments provided by docopt.
    """
    def __init__(self, arguments, current_directory):
        self._arguments = arguments
        self._current_directory = current_directory

    @property
    def command(self):
        # We have to do this weird loop to deal with the way docopt stores the command name
        for possible_command in ('map',
                                 'init',
                                 'h5',
                                 'align',
                                 'info',
                                 'notebooks'):
            if self._arguments.get(possible_command):
                return possible_command

    @property
    def fastq_directory(self):
        return self._arguments['FASTQ_DIRECTORY']

    @property
    def log_level(self):
        log_level = {0: logging.ERROR,
                     1: logging.WARN,
                     2: logging.INFO,
                     3: logging.DEBUG}
        # default to silent if the user supplies no verbosity setting
        return log_level.get(self._arguments['-v'], logging.ERROR)

    @property
    def output_directory(self):
        return self._arguments['OUTPUT_DIRECTORY']

    @property
    def process_limit(self):
        # 0 indicates unlimited
        return int(self._arguments['--process-limit'] or 0)


class PathInfo(object):
    """ Parses user-provided alignment parameters and provides a default in case no value was given. """
    def __init__(self, image_directory, mapped_reads, perfect_target_name, alternate_fiducial_reads=None,
                 alternate_perfect_reads_filename=None, alternate_good_reads_filename=None):
        self._alternate_fiducial_reads = alternate_fiducial_reads
        self._alternate_good_reads_filename = alternate_good_reads_filename
        self._alternate_perfect_reads_filename = alternate_perfect_reads_filename
        self._image_directory = image_directory
        self._mapped_reads = mapped_reads
        self._perfect_target_name = perfect_target_name

    @property
    def aligning_read_names_filepath(self):
        if self._alternate_fiducial_reads:
            return os.path.join(self._mapped_reads, self._alternate_fiducial_reads)
        return os.path.join(self._mapped_reads, 'phix_read_names.txt')

    @property
    def all_read_names_filepath(self):
        return os.path.join(self._mapped_reads, 'all_read_names.txt')

    @property
    def figure_directory(self):
        return os.path.join(self._image_directory, 'figs')

    @property
    def on_target_read_names(self):
        if self._alternate_good_reads_filename:
            return os.path.join(self._mapped_reads, self._alternate_good_reads_filename)
        if not self._perfect_target_name:
            raise ValueError("This experiment did not have a perfect target set!")
        return os.path.join(self._mapped_reads, 'target_{}_read_names.txt'.format(self._perfect_target_name.lower()))

    @property
    def perfect_read_names(self):
        if self._alternate_perfect_reads_filename:
            return os.path.join(self._mapped_reads, self._alternate_perfect_reads_filename)
        if not self._perfect_target_name:
            raise ValueError("This experiment did not have a perfect target set!")
        return os.path.join(self._mapped_reads, 'perfect_target_{}_read_names.txt'.format(self._perfect_target_name.lower()))

    @property
    def results_directory(self):
        return os.path.join(self._image_directory, 'results')
