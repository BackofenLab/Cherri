# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='RRIevaluator',
    version='1.0.0',
    description='RRIevaluator is a tool for classification of RNA RNA interactions.',
    long_description=readme,
    author='Teresa Mueller',
    author_email='muellert@informatik.uni-freiburg.de',
    url='https://github.com/teresa-m/RNA_RNA_binding_evaluation',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)
