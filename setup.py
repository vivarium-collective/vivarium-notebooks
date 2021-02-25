import os
import glob
import setuptools
from distutils.core import setup

with open("README.md", 'r') as readme:
    long_description = readme.read()

setup(
    name='bioscrape-cobra',
    version='0.0.1',
    packages=[
        'bioscrape_cobra',
        'bioscrape_cobra.processes',
        'bioscrape_cobra.composites',
        'bioscrape_cobra.models',
    ],
    author='Eran Agmon, William Poole',
    author_email='',
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
        'vivarium-core',
        'pytest',
    ],
)
