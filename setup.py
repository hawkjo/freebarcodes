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
        include_dirs=[np.get_include()],
        zip_safe=False,
        author='John Hawkins',
        author_email='hawkjo@gmail.com',
        description='FREE Divergence Error-Correcting DNA Barcodes',
        url='https://github.com/hawkjo/freebarcodes',
        download_url='https://github.com/hawkjo/freebarcodes/archive/refs/tags/v2.1.tar.gz',
        keywords=['DNA', 'NGS', 'bioinformatics', 'barcodes'],
        python_requires='>=2.7,!=3.*',
        classifiers=['Development Status :: 3 - Alpha',
                     'Natural Language :: English',
                     'Intended Audience :: Science/Research',
                     'Operating System :: POSIX :: Linux',
                     'Programming Language :: Python :: 2.7',
                     'Topic :: Scientific/Engineering :: Bio-Informatics',
                     ]
    )
