# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='Cherri',
    version='0.1',
    description='Cherri is a tool for classification of RNA RNA interactions.',
    long_description=readme,
    author='Teresa Mueller',
    author_email='muellert@informatik.uni-freiburg.de',
    url='https://github.com/teresa-m/Cherri',
    license=license,
    scripts=['bin/peakhood', 'bin/gtf_add_transcript_biotype_info.py',
             'bin/bed_generate_unique_ids.py'],
    packages=find_packages(exclude=('tests', 'docs'))
)
