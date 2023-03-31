import os

from setuptools import find_packages, setup

with open('README.rst', 'r') as f:
    readme = f.read()

with open('CHANGELOG.rst', 'r') as f:
    changes = f.read()

def parse_requirements(filename):
    ''' Load requirements from a pip requirements file '''
    with open(filename, 'r') as fd:
        lines = []
        for line in fd:
            line.strip()
            if line and not line.startswith("#"):
                lines.append(line)
    return lines

requirements = parse_requirements('requirements.txt')


if __name__ == '__main__':
    setup(
        name='ctxcore',
        use_scm_version=True,
        setup_requires=['setuptools_scm'],
        description='Core functions for pycisTarget and the SCENIC tool suite',
        long_description='\n\n'.join([readme, changes]),
        license='GNU General Public License v3',
        url='https://github.com/aertslab/ctxcore',
        author='Bram Van de Sande',
        maintainer='Christopher Flerin',
        maintainer_email='christopher.flerin@kuleuven.be',
        install_requires=requirements,
        keywords=['ctxcore','cisTarget','pycisTarget','SCENIC','pySCENIC'],
        include_package_data=True,
        package_dir={'': 'src'},
        packages=find_packages('src'),
        zip_safe=False,
        classifiers=['Development Status :: 3 - Alpha',
                     'Intended Audience :: Developers',
                     'Programming Language :: Python :: 3.6',
                     'Programming Language :: Python :: 3.7',
                     'Programming Language :: Python :: 3.8']
    )
