# MetaClassifier
## Overview
MetaClassifier is an integrated pipeline for identifying the floral composition of honey using DNA metabarcoding to determine the plants that honey bees visit. MetaClassifier utilizes a database of marker sequences and their corresponding taxonomy lineage information to classify high-throughput metabarcoding sample sequencing reads data into taxonomic groups and quantify taxon abundance. MetaClassifier can also be employed in other studies that utilize barcoding, metabarcoding, and metagenomics techniques to characterize richness, abundance, relatedness, and interactions in ecological communities.

In addition to this README file, you can consult the MetaClassifier [manual](docs/MetaClassifier.md) for more detailed information.

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

## Using MetaClassifier
- **MetaClassifier pipeline package**: 
  - Display all usage options: 
    ```
    python3 -m metaclassifier -h
    ```
  - Basic usage with defaults for optional arguments:
    ```
    python3 -m metaclassifier [options] <SAMPLE_FILE> <DB_DIR> <CONFIG_FILE>
    ```
- **MetaClassifier pipeline wrapper script**: 
  - Display all usage options: 
    ```
    python3 MetaClassifier/metaclassifier.py -h
    ```
  - Basic usage with defaults for optional arguments:
    ```
    python3 MetaClassifier/metaclassifier.py [options] <SAMPLE_FILE> <DB_DIR> <CONFIG_FILE>
    ```
- **MetaClassifier read processing script**: 
  - Display all usage options: 
    ```
    python3 MetaClassifier/metaclassifier/process_reads.py -h
    ```
  - Basic usage with defaults for optional arguments:
    ```
    python3 MetaClassifier/metaclassifier/process_reads.py [options] <SAMPLE_FILE>
    ```
- **MetaClassifier read classification script**: 
  - Display all usage options: 
    ```
    python3 MetaClassifier/metaclassifier/classify_reads.py -h
    ```
  - Basic usage with defaults for optional arguments:
    ```
    python3 MetaClassifier/metaclassifier/classify_reads.py [options] <FASTA_DIR> <DB_DIR> <CONFIG_FILE>
    ```
Where:
* `<SAMPLE_FILE>` is the input tab-separated file specifying 1) sample names, 2) file names for forward paired-end reads, and 3) file names for reverse paired-end reads. The full file path is required if the files are not in the current directory, and the second file is not required for single-end frangments. [Example sample input file is available here](test/mv sample_input.tsv).
* `<DB_DIR>` is the input marker database directory with sequence fasta files and corresponding taxonomy lineage files for each marker. Marker data files should be named with marker nomenclature followed `.fa` for the fasta sequence file and `.tax` for the taxonomy lineage file. (i.e., `ITS1.fa` and `ITS1.tax` for the first internal transcribed spacer in eukaryotic ribosomal RNA subunit genes). The marker taxonomy lineage files should formated as a tab-separated file with 1) the NCBI taxon ID,  2) order name, 3) family name, 4) genus name, and 5) species name as shown in the MetaClassifier's [MetaCurator reformated reference database](http://bigdata.bx.psu.edu/MetaClassifier_databases/).
* `<CONFIG_FILE>` is the input tab-separated file specifying marker name `(i.e., ITS1)` and its corresponding VSEARCH's usearch_global function search settings of minimum query coverage `(i.e., 0.8 for 80)` and minimum sequence identity `(i.e., 0.95 for 95%)` for each search marker. [Example sample input file is available here](test/sample_config.tsv).
* `<FASTA_DIR>` is the directory containing samples read data fasta files produced by either the MetaClassifier's `process_reads.py` script or from external sources for classification by the MetaClassifier's `classify_reads.py` script. 

## MetaClassifier output
`python ../metaclassifier.py -m -r ../bin/pear -s ../bin/seqtk -a ../bin/vsearch sample_input.tsv ../db/MetabarcodeDBsV2 sample_config.tsv`

run log:
```
===========================================================
MetaClassifier version 1.0.0 (release date: 05 March 2021)
===========================================================

17:22:40 - Starting MetaClassifier...

17:22:40 - - Merging paired-end (PE) sample read dataset(s)...

17:23:36 - - Converting merged paired-end (PE) sample read dataset(s) from FASTQ to FASTA...

17:23:36 - - Searching for markers in the sample read dataset(s)...

17:23:36 - - - Searching sample sample1...

17:24:51 - - - - Computing genus taxonomy class read proportions for ITS1 marker...

17:27:27 - - - - Computing genus taxonomy class read proportions for ITS2 marker...

17:28:59 - - - - Computing genus taxonomy class read proportions for trnL marker...

17:29:00 - - - - Computing and rescaling median proportions for genus taxonomy class across all markers...

17:29:00 - - - Searching sample sample2...

17:30:17 - - - - Computing genus taxonomy class read proportions for ITS1 marker...

17:32:54 - - - - Computing genus taxonomy class read proportions for ITS2 marker...

17:34:29 - - - - Computing genus taxonomy class read proportions for trnL marker...

17:34:29 - - - - Computing and rescaling median proportions for genus taxonomy class across all markers...

17:34:29 - Completed MetaClassifier...
Total elapsed time 708
```


## Citation
If you use MetaClassifier please cite the following paper the describes the methodology:

**Characterizing the floral resources of a North American metropolis using a honey bee foraging assay.**  
_Douglas B. Sponsler_, _Don Shump_,  _Rodney T. Richardson_,  _Christina M. Grozinger_.  
Ecosphere 11, no. 4 (2020): e03102.   
DOI: [https://doi.org/10.1002/ecs2.3102](https://doi.org/10.1002/ecs2.3102)

## License
MetaClassifier is distributed under the GNU GPL v3.0 For more information, see [license](LICENSE).

