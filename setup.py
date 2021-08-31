import os
import glob
import setuptools
from distutils.core import setup

with open("README.md", 'r') as readme:
    long_description = readme.read()

setup(
    name='vivarium-notebooks',
    version='0.0.4',
    packages=[
        'bioscrape_cobra',
    ],
    author='Eran Agmon, William Poole',
    author_email='agmon.eran@gmail.com, wpoole@caltech.edu',
    url='https://github.com/vivarium-collective/vivarium-notebooks',
    license='MIT',
    entry_points={
        'console_scripts': []},
    short_description='Jupyter notebooks and supporting functions for the simulations '
                      'described in the paper "Vivarium: an interface and engine for '
                      'integrative multiscale modeling in computational biology"',
    long_description=long_description,
    long_description_content_type='text/markdown',
    package_data={},
    # package_data={
    #     # If any package contains *.xml files, include them:
    #     "": ["*.xml"],
    # },
    include_package_data=True,
    install_requires=[
        'vivarium-core>=0.3.5',
        'vivarium-bioscrape==0.0.0.7',
        'vivarium-cobra==0.0.18',
        'vivarium-multibody==0.0.13',
        'pytest',
        'tqdm',
        'jupyter'
    ],
)
