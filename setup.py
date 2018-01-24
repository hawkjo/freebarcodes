from setuptools import setup
from freebarcodes.constants import VERSION
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np


if __name__ == '__main__':
    setup(
        name='freebarcodes',
        packages=['freebarcodes'],
        version=VERSION,
        entry_points={
          'console_scripts': [
              'freebarcodes = freebarcodes.main:main'
          ]
        },
        ext_modules=cythonize('freebarcodes/editmeasures.pyx'),
        include_package_data=True,
        zip_safe=False,
        description='FREE Divergence Error-Correcting DNA Barcodes',
        url='http://www.finkelsteinlab.org',
        keywords=['DNA', 'NGS', 'bioinformatics', 'barcodes'],
        classifiers=['Development Status :: 3 - Alpha',
                     'Natural Language :: English',
                     'Intended Audience :: Science/Research',
                     'License :: Freely Distributable',
                     'Operating System :: POSIX :: Linux',
                     'Programming Language :: Python :: 2.7',
                     'Topic :: Scientific/Engineering :: Bio-Informatics',
                     ]
    )
