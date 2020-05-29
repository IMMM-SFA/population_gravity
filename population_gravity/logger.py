"""
Logger for the popultation_gravity model

@author   Chris R. Vernon
@email:   chris.vernon@pnnl.gov

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

"""

import logging
import sys

from population_gravity.read_config import ReadConfig


class Logger(ReadConfig):
    """Initialize project-wide logger. The logger outputs to both stdout and a file."""

    # output format for log string
    LOG_FORMAT_STRING = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    @property
    def log_format(self):
        """Generate log formatter."""

        return logging.Formatter(self.LOG_FORMAT_STRING)

    @property
    def logger(self):
        """Initialize logger as level info."""

        logger = logging.getLogger()
        logger.setLevel(logging.INFO)

        return logger

    def initialize_logger(self):
        """Initialize logger to stdout and file."""

        # logger console handler
        self.console_handler()

        # logger file handler
        self.file_handler()

    def console_handler(self):
        """Construct console handler."""

        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(self.log_format)
        self.logger.addHandler(console_handler)

    def file_handler(self):
        """Construct file handler."""

        file_handler = logging.FileHandler(self.logfile)
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(self.log_format)
        self.logger.addHandler(file_handler)

    @staticmethod
    def close_logger():
        """Shutdown logger."""

        # Remove logging handlers
        logger = logging.getLogger()

        for handler in logger.handlers[:]:
            handler.close()
            logger.removeHandler(handler)

        logging.shutdown()
