# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='cherri',
    version='0.1',
    description='Cherri is a tool for classification of RNA-RNA interactions.',
    long_description=readme,
    author='Teresa Mueller',
    author_email='muellert@informatik.uni-freiburg.de',
    url='https://github.com/teresa-m/Cherri',
    license=license,
    scripts=['bin/cherri', 'bin/find_occupied_regions.py',
             'bin/find_trusted_RRI.py', 'bin/generate_pos_neg_with_context.py',
             'bin/get_features.py'],
    packages=find_packages(exclude=('tests', 'docs')),
    package_data={'': ['IntaRNA_param/*']},
    include_package_data=True,
    #packages=['rrieval'],
    # package_data={'rrieval': ['lib.py','IntaRNA_param/IntaRNA_param.txt']},
    zip_safe=False
)
