import os
import glob
import setuptools
from distutils.core import setup

with open("README.md", 'r') as readme:
    long_description = readme.read()

setup(
    name='vivarium-notebooks',
    version='0.0.1',
    packages=[],
    author='Eran Agmon, William Poole',
    author_email='eagmon@stanford.edu, wpoole@caltech.edu',
    url='https://github.com/vivarium-collective/vivarium-notebooks',
    license='',
    entry_points={
        'console_scripts': []},
    short_description='Includes composite models of vivarium-bioscrape and vivarium-cobra',
    long_description=long_description,
    long_description_content_type='text/markdown',
    package_data={},
    include_package_data=True,
    install_requires=[
        'vivarium-core==0.2.8',
        'vivarium-bioscrape==0.0.0.6',
        'vivarium-cobra==0.0.8',
        'vivarium-multibody==0.0.7',
        'pytest',
        'tqdm',
    ],
)
