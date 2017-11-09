from setuptools import setup
from nucleaseq.constants import VERSION
from distutils.extension import Extension
#from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy as np


if __name__ == '__main__':
    setup(
        name='nucleaseq',
        packages=['nucleaseq'],
        version=VERSION,
        #entry_points={
          #'console_scripts': [
              #'nucleaseq = nucleaseq.main:main'
          #]
        #},
        ext_modules=cythonize('nucleaseq/editmeasures.pyx'),
        include_package_data=True,
        zip_safe=False,
        description='Processes Nucleaseq data',
        url='http://www.finkelsteinlab.org',
        keywords=['DNA', 'protein', 'NGS', 'bioinformatics', 'crispr'],
        classifiers=['Development Status :: 3 - Alpha',
                     'Natural Language :: English',
                     'Intended Audience :: Science/Research',
                     'License :: Freely Distributable',
                     'Operating System :: POSIX :: Linux',
                     'Programming Language :: Python :: 2.7',
                     'Topic :: Scientific/Engineering :: Bio-Informatics',
                     ]
    )
