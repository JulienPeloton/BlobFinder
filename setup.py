#!/usr/bin/env python
# Copyright (C) 2016 Peloton
"""Distutils based setup script for BlobFinder

For the easiest installation just type the command (you'll probably need
root privileges for that):

    python setup.py install

This will install the library in the default location. For instructions on
how to customize the install procedure read the output of:

    python setup.py --help install

In addition, there are some other commands:

    python setup.py clean -> will clean all trash (*.pyc and stuff)
    python setup.py test  -> will run the complete test suite
    python setup.py bench -> will run the complete benchmark suite
    python setup.py audit -> will run pyflakes checker on source code

To get a full list of avaiable commands, read the output of:

    python setup.py --help-commands
"""
from setuptools import find_packages
from setuptools import setup

if __name__ == "__main__":
    setup(name='BlobFinder',
        version='0.1.0',
        author='Julien Peloton',
        author_email='j.peloton@sussex.ac.uk',
        url='https://github.com/JulienPeloton/BlobFinder',
        install_requires=[''],
        packages=['BlobFinder'],
        description='Find point sources in sky maps using a match filtering algorithm',
        classifiers=[
            "Programming Language :: Python :: 2",
            'Programming Language :: Python :: 2.7'],)
