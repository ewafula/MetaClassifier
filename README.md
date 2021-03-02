# MetaClassifier
## Overview
MetaClassifier is an integrated pipeline for identifying the floral composition of honey using DNA metabarcoding to determine the plants that honey bees visit. MetaClassifier utilizes a database of marker sequences and their corresponding taxonomy lineage information to classify high-throughput metabarcoding sample sequencing reads data into taxonomic groups and quantify taxon abundance. MetaClassifier can also be employed in other studies that utilize barcoding, metabarcoding, and metagenomics techniques to characterize richness, abundance, relatedness, and interactions in ecological communities.

In addition to this README file, you can consult the MetaClassifier [manual](https://github.com/ewafula/MetaClassifier/blob/main/docs/MetaClassifier.md) for more detailed information.

## Installation
MetaClassifier requires dependencies and external tools that need to be installed and available on the environment the pipeline can be used. Not a requirement if insatalling using Bioconda.

### Dependecies

* [Python](https://www.python.org/) (version >=3.7)
* [pandas](https://pandas.pydata.org/) (version >=1.2.2)

### External tools
* [PEAR](https://cme.h-its.org/exelixis/web/software/pear/) for merging overlapping paired-end (PE) reads
* [seqtk](https://github.com/lh3/seqtk/) for converting FASTQ to FASTA sequence format
* [VSEARCH](https://github.com/torognes/vsearch/) for searching searching high-throughtput sequence read data against marker database

### Python package installation
MetaClassifier is available through [pypi](https://pypi.python.org/pypi/metaclassifier). To install, type:
```
pip install metaclassifier
```
### Repository from GitHub
```
Either "git clone https://github.com/ewafula/MetaClassifier.git" or download and "unzip MetaClassifier-main.zip"
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
## Marker reference databases
[MetaCurator](https://github.com/RTRichar/MetaCurator) reference databases with taxonomy lineage information reformated to work with MetaClassifier. 
Detailed step by step tutorial workflow for creating reference marker database is descrribe on the GitHub [MetaCurator database repository](https://github.com/RTRichar/MetabarcodeDBsV2/blob/master/Workflow.md)
* [MetabarcodeDBs](http://bigdata.bx.psu.edu/MetaClassifier_databases/)

## Basic usage
```
metaclassifier [options] <SAMPLE_FILE> <DB_DIR> <CONFIG_FILE>
```

Please consult the MetaClassifier [manual](https://github.com/ewafula/MetaClassifier/blob/main/docs/MetaClassifier.md) for a detailed description and usage of all options.

## Citation
If you use MetaClassifier please cite the following paper the describes the methodology:

**Characterizing the floral resources of a North American metropolis using a honey bee foraging assay.**  
_Douglas B. Sponsler_, _Don Shump_,  _Rodney T. Richardson_,  _Christina M. Grozinger_.  
Ecosphere 11, no. 4 (2020): e03102.   
DOI: [https://doi.org/10.1002/ecs2.3102](https://doi.org/10.1002/ecs2.3102)

## License
MetaClassifier is distributed under the GNU GPL v3.0 For more information, see [license](https://github.com/ewafula/MetaClassifier/blob/main/LICENSE).
