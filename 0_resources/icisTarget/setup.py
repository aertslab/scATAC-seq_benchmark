#!/usr/bin/env python
import os
import numpy
from setuptools import setup, Extension, find_packages

from Cython.Distutils import build_ext


setup(
    name='cistargetx',
    version='3.0',
    description='Motif, enhancer, regulatory track and target gene prediction in Drosophila, Human and Mouse',
    author='Bram Van de Sande',
    author_email='bram.vandesande@med.kuleuven.be',
    maintainer='Gert Hulselmans',
    maintainer_email='gert.hulselmans@kuleuven.be',
    url='http://gbiomed.kuleuven.be/english/research/50000622/lcb/tools',
    ext_modules=[
        Extension("cistargetx.common.orderstatistics",
                  ["cistargetx/common/orderstatistics.pyx"],
                  include_dirs=[numpy.get_include()],
                 ),
    ],
    cmdclass={'build_ext': build_ext},
    packages=find_packages(),
    package_data={
        '': ['*.bed', '*.ini', '*.vm'],
        'cistargetx.daemon': ['*.gif', '*.css', '*.sql'],
        'cistargetx.iregulon': ['*.sql'],
        'cistargetx.cistargetdb': ['*.sql'],
        'cistargetx.pwmrankings': ['*.jar'],
        'cistargetx.web': ['*.php', '*.js', '*.css', '.htaccess', '*.html'],
        'cistargetx.web.images': ['*.gif', '*.jpg'],
        'unittests.common': ['*.rr', '*.r', '*.hgnc'],
        'unittests.recoveryanalysis': ['*.hgnc'],
    },
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'ctdb-query = cistargetx.cistargetdb.query:main',
            'ctdb-upload = cistargetx.cistargetdb.upload2database:main',
            'ctx-check-availability = cistargetx.daemon.checkavailability:main',
            'ctx-convert = cistargetx.conversion.convert:main',
            'ctx-convids = cistargetx.miscellaneous.convertfeatureids:main',
            'ctx-daemon = cistargetx.daemon.cistargetxdaemon:main',
            'ctx-db2r = cistargetx.miscellaneous.db2rankings:main',
            'ctx-dbadd = cistargetx.miscellaneous.adddatabases:main',
            'ctx-dbfilter = cistargetx.miscellaneous.filterdatabase:main',
            'ctx-dbrandomize = cistargetx.miscellaneous.randomizedb:main',
            'ctx-dbrm = cistargetx.miscellaneous.rmfeatures:main',
            'ctx-gt = cistargetx.pwmrankings.genetrees:main',
            'ctx-pwm2r = cistargetx.pwmrankings.createrankings:main',
            'ctx-r2db = cistargetx.miscellaneous.rankings2db:main',
            'ctx-r2gmt = cistargetx.miscellaneous.r2gmt:main',
            'ctx-rcc = cistargetx.recoveryanalysis.recoveryanalysis:main',
            'ctx-recombine-rankings = cistargetx.miscellaneous.recombinerankings:main',
            'ctx-track2db = cistargetx.miscellaneous.track2database:main',
            'iregulon-check-availability = cistargetx.iregulon.checkavailability:main',
            'iregulon-daemon = cistargetx.iregulon.daemon:main',
            'iregulon-upload = cistargetx.iregulon.uploadmetatargetome:main',
        ]
    },
    # install_requires=[
    #        'Cython >= 0.14.1',
    #        'numpy >= 1.5.1',
    #        'matplotlib >= 0.99.3',
    #        'scipy >= 0.8.0',
    #        'MySQL-python >= 1.2.2',
    #        'bx-python >= 0.6.0',
    #        'airspeed == 0.1dev-20101001'
    #        ],
    test_suite='unittests'
)
