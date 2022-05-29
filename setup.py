#!/usr/bin/env python3

from setuptools import setup, find_packages
from os import path

PKG_NAME = "OffTarget"
MOD_NAME = "OffTarget"

DESCRIPTION = """ 
Off-Target is a short-string (<50bp) search query tool which maps queries to hits with high mismatches or Insertion Deletions (InDels) with a default setting of 1 mismatch per 3 nucleotides. 
"""

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md")) as fh, open(path.join(here, "requirements.txt")) as req:
    install_requires = [pkg.strip() for pkg in req]

__version__ = ""

exec(open("{}/_version.py".format(MOD_NAME)).read())

setup(
    name=PKG_NAME,
    version=__version__,
    author="M Irfan",
    author_email="bioirfanatics@gmail.com",
    description=f"{PKG_NAME}",
    long_description=DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/chewlabSB2/OffTarget",
    packages=find_packages(exclude=["*.test", "*.test.*", "test.*", "test"]),
    package_dir={'OffTarget': 'OffTarget'}, 
    entry_points={
        "console_scripts": [
            "OffTarget={}.query:main".format(MOD_NAME),
        ],
    },
    
    install_requires=install_requires,
    include_package_data=True,
    python_requires=">=3.5",
)
