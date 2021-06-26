from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np

import codecs
import os.path

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


if __name__ == '__main__':
    setup(
        name='freebarcodes',
        packages=['freebarcodes'],
        version=get_version("freebarcodes/__init__.py"),
        entry_points={
          'console_scripts': [
              'freebarcodes = freebarcodes.main:main'
          ]
        },
        ext_modules=cythonize('freebarcodes/editmeasures.pyx'),
        include_package_data=True,
        include_dirs=[np.get_include()],
        install_requires=[
            "pathos>=0.2.1",
            "psutil>=5.8.0",
            "h5py>=3.2.1",
            "numpy>=1.19.0",
            "docopt>=0.6.2",
            "biopython>=1.70",
            "cython>=0.29.23",
            ],
        zip_safe=False,
        author='John Hawkins',
        author_email='hawkjo@gmail.com',
        description='FREE Divergence Error-Correcting DNA Barcodes',
        url='https://github.com/hawkjo/freebarcodes',
        download_url='https://github.com/hawkjo/freebarcodes/archive/refs/tags/v3.0.tar.gz',
        keywords=['DNA', 'NGS', 'bioinformatics', 'barcodes'],
        python_requires='>=3.0',
        classifiers=['Development Status :: 3 - Alpha',
                     'Natural Language :: English',
                     'Intended Audience :: Science/Research',
                     'Operating System :: POSIX :: Linux',
                     'Programming Language :: Python :: 3',
                     'Topic :: Scientific/Engineering :: Bio-Informatics',
                     ]
    )
