# MetaClassifier
MetaClassifier is an integrated pipeline for identifying the floral composition of honey using DNA metabarcoding to determine the plants that honey bees visit. MetaClassifier utilizes a database of marker sequences and their corresponding taxonomy information to classify high-throughput metabarcoding sample sequencing reads data into taxonomic groups and quantify taxon abundance. MetaClassifier can also be employed in other studies that utilize barcoding, metabarcoding, and metagenomics techniques to characterize richness, abundance, relatedness, and interactions in ecological communities.

In addition to this README file, you can consult the MetaClassifier [manual](docs/MetaClassifier.md) for more detailed information.

## Installation
MetaClassifier requires dependencies and external tools that need to be installed and available on the environment the pipeline can be used. Not a requirement if insatalling using Bioconda.

### Dependecies

* [Python](https://www.python.org/) (version >=3.6)
* [pandas](https://pandas.pydata.org/) (version >=1.2.2)

### External tools
* [PEAR](https://cme.h-its.org/exelixis/web/software/pear/) for merging overlapping paired-end (PE) reads
* [seqtk](https://github.com/lh3/seqtk/) for converting FASTQ to FASTA sequence format
* [VSEARCH](https://github.com/torognes/vsearch/) for searching searching high-throughtput sequence read data against marker database

### Python package installation
MetaClassifier is available through [pypi] (https://pypi.python.org/pypi/metaclassifier). To install, type:
```
pip install metaclassifier
```
### Repository from GitHub
```
git clone https://github.com/ewafula/MetaClassifier.git
cd MetaClassifier/
python setup.py install
```
### Bioconda package
This requires a working [Conda](https://docs.conda.io/en/latest/miniconda.html#) installation.
```
conda install -c bioconda metaclassifier
```
We recommend installing MetaClassifier in a new separate environment from the base for all dependencies to be properly resolved by conda. To install, type:
```
conda create -n "metaclassifier" -c bioconda metaclassifier=1.0.0
```



