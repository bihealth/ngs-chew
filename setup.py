#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Installation driver (and development utility entry point) for ngs-chew
"""

import os

from setuptools import find_packages, setup

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


def parse_requirements(path):
    """Parse ``requirements.txt`` at ``path``."""
    requirements = []
    with open(path, "rt") as reqs_f:
        for line in reqs_f:
            line = line.strip()
            if line.startswith("-r"):
                fname = line.split()[1]
                inner_path = os.path.join(os.path.dirname(path), fname)
                requirements += parse_requirements(inner_path)
            elif line != "" and not line.startswith("#"):
                requirements.append(line)
    return requirements


with open("README.md") as readme_file:
    readme = readme_file.read()

with open("CHANGELOG.md") as history_file:
    history = history_file.read()

test_requirements = parse_requirements("requirements/test.txt")
install_requirements = parse_requirements("requirements/base.txt")
serve_requirements = parse_requirements("requirements/serve.txt")


def bash_scripts(names):
    """Return generator with bash scripts of the given ``names``"""
    return (os.path.join("scripts", name) for name in names)


package_root = os.path.abspath(os.path.dirname(__file__))
version = {}
with open(os.path.join(package_root, "chew/_version.py")) as fp:
    exec(fp.read(), version)
version = version["__version__"]

setup(
    author="Manuel Holtgrewe",
    author_email="manuel.holtgrewe@bih-charite.de",
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    description="NGS Chew",
    entry_points={"console_scripts": (("ngs-chew=chew.cli:cli",),)},
    install_requires=install_requirements,
    tests_require=test_requirements,
    extras_require={"serve": serve_requirements},
    license="MIT license",
    long_description=readme + "\n\n" + history,
    long_description_content_type="text/markdown",
    include_package_data=True,
    keywords="bioinformatics",
    name="ngs-chew",
    packages=find_packages(),
    package_dir={"chew": "chew"},
    test_suite="tests",
    url="https://github.com/bihealth/ngs-chew",
    version=version,
    zip_safe=False,
)
