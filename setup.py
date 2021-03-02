#!/usr/bin/env python

import pathlib
from setuptools import setup, find_packages

# setup.py directory also contains README.md
setup_dir = pathlib.Path(__file__).parent
README = (setup_dir / "README.md").read_text()


setup(
    name="metaclassifier",
    version="1.0.0",
    description="MetaClassifier is an integrated pipeline for classifying and quantifying DNA metabarcoding data into taxonomy groups",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/ewafula/MetaClassifier",
    author='Eric Wafula',
    author_email="ekw10@gmail.com",
    license="GNU General Public License v3 (GPLv3)",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
    ],
    packages=find_packages(exclude=["venv", "docs", "test", "bin", "db"]),
    include_package_data=True,
    install_requires=["pandas"],
    entry_points={
        "console_scripts": [
            "metaclassifier=metaclassifier.__main__:main",
            "process_reads=metaclassifier.process_reads.py:main",
            "classify_reads.py=metaclassifier.classify_reads.py:main"
        ],
    },
    zip_safe=True
)
