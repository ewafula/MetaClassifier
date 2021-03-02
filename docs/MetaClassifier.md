# MetaClassifier
## Overview
MetaClassifier is an integrated pipeline for identifying the floral composition of honey using DNA metabarcoding to determine the plants that honey bees visit. MetaClassifier utilizes a database of marker sequences and their corresponding taxonomy lineage information to classify high-throughput metabarcoding sample sequencing reads data into taxonomic groups and quantify taxon abundance. MetaClassifier can also be employed in other studies that utilize barcoding, metabarcoding, and metagenomics techniques to characterize richness, abundance, relatedness, and interactions in ecological communities.

In addition to this README file, you can consult the MetaClassifier [manual](../docs/MetaClassifier.md) for more detailed information.

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
conda create -n "metaclassifier" -c bioconda metaclassifier=1.0.1
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
* `<SAMPLE_FILE>` is the input tab-separated file specifying 1) sample names, 2) file names for forward paired-end reads, and 3) file names for reverse paired-end reads. The full file path is required if the files are not in the current directory, and the second file is not required for single-end frangments. [Example sample input file is available here](../test/sample_input.tsv).
* `<DB_DIR>` is the input marker database directory with sequence fasta files and corresponding taxonomy lineage files for each marker. Marker data files should be named with marker nomenclature followed `.fa` for the fasta sequence file and `.tax` for the taxonomy lineage file. (i.e., `ITS1.fa` and `ITS1.tax` for the first internal transcribed spacer in eukaryotic ribosomal RNA subunit genes). The marker taxonomy lineage files should formated as a tab-separated file with 1) the NCBI taxon ID,  2) order name, 3) family name, 4) genus name, and 5) species name as shown in the MetaClassifier's [MetaCurator reformated reference database](http://bigdata.bx.psu.edu/MetaClassifier_databases/).
* `<CONFIG_FILE>` is the input tab-separated file specifying marker name `(i.e., ITS1)` and its corresponding VSEARCH's usearch_global function search settings of minimum query coverage `(i.e., 0.8 for 80)` and minimum sequence identity `(i.e., 0.95 for 95%)` for each search marker. [Example sample input file is available here](../test/sample_config.tsv).
* `<FASTA_DIR>` is the directory containing sample read data FASTA produced by either the MetaClassifier's `process_reads.py` script or from an external source for classification by the MetaClassifier's `classify_reads.py` script. 

## MetaClassifier output
As an example, this section uses the a small test dataset subsampled from the sample read data that was utilized in the [Sponsler et al., 2020](#citation) study to show how to perform an complete analysis using wrapper script. The test datasets and [external tools](#external-tools) Linux binaries are located in the [test](../test) and [bin](../bin) sub-directory of MetaClassifier installation respectively. 

**Analysis**:
- Get into the [test](../test) directory of MetaClassifier:
  - `cd MetaClassifier/test/`
- Unpack the `sample1` and `sample2` test paired-end (PE) read files:
  - `gzip -d *.gz`  
- Download and unpack the [MetaCurator reference database](http://bigdata.bx.psu.edu/MetaClassifier_databases/) in the [db](../db) directory of MetaClassifier:
  - `wget -qO- http://bigdata.bx.psu.edu/MetaClassifier_databases/MetabarcodeDBsV2.tar.gz | tar -xvz -C ../db/`
- Change persmissions to make [external tools](#external-tools) Linux binaries in the [bin](../bin) directory executable:
  - `chmod u+x ../bin/*`
- Execute the `metaclassifier.py` wrapper script with defaults for optional arguments and allowing PE reads merging:
  - `python ../metaclassifier.py -m -r ../bin/pear -s ../bin/seqtk -a ../bin/vsearch sample_input.tsv ../db/MetabarcodeDBsV2 sample_config.tsv`

**Run log**:
```
===========================================================
MetaClassifier version 1.0.1 (release date: 05 March 2021)
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

**Outputs**:

All of the output will be in the `sample_input/` output directory, named after the sample input file name if an alternative output name is not provided by the `-o` option. If running the complete pipeline using the `metaclassifier.py` wrapper script as shown above, the contents of the output directory is as follows:
```
sample_input/
├── merge_dir/
|           ├── sample1_pear.log
|           ├── sample1.assembled.fastq
|           ├── sample1.discarded.fastq
|           ├── sample1.unassembled.forward.fastq
|           ├── sample1.unassembled.reverse.fastq
|           ├── sample2_pear.log
|           ├── sample2.assembled.fastq
|           ├── sample2.discarded.fastq
|           ├── sample2.unassembled.forward.fastq
|           └── sample2.unassembled.reverse.fastq
├── fasta_dir/
|           ├── sample1.fasta
|           └── sample2.fasta
├── sample1/
|         ├── ITS1_vsearch.log
|         ├── ITS1_blast6out.tsv
|         ├── ITS2_vsearch.log 
|         ├── ITS2_blast6out.tsv
|         ├── trnL_vsearch.log
|         ├── trnL_blast6out.tsv
|         └── sample1_rescaled_propotions.tsv
├── sample2/
|         ├── ITS1_vsearch.log
|         ├── ITS1_blast6out.tsv
|         ├── ITS2_vsearch.log 
|         ├── ITS2_blast6out.tsv
|         ├── trnL_vsearch.log
|         ├── trnL_blast6out.tsv
|         └── sample2_rescaled_propotions.tsv
└── sample.tsv
```
Where:  
* **`merge_dir/`**: this output sub-directory contains results of the sample paired-end (PE) read merging process, including `*_pear.log`: the PEAR merger run logs, `*.assembled.fastq`: the merged overlapping PE reads, `*.discarded.fastq`: discarded erroneous reads, `*.unassembled.forward.fastq`: unmerged forward fragments of PE reads, and `*.unassembled.reverse.fastq`: unmerged reverse fragments of PE reads. 
* **`fasta_dir/`**: this output sub-directory contains the merged overlapping paired-end reads sample FASTA converted from FASTQ format. This out directory be used as an input of the `classify_reads.py` script for reruns to optimize classification parameters.
* **`sample1&2/`** this output sub-directories contains results of the each sample's paired-end reads classification process, including `*_vsearch.log`: the reference marker VSEARCH run logs, `*_blast6out.tsv`: the reference marker VSEARCH search results, and `*_rescaled_propotions.tsv`: the identified sample taxonomic composition and taxon read abundance propotions. 
* **`sample.tsv`** this a recreated sample input file for merged sample paired-end reads.

## Parameters

### metaclassifier.py
```
usage: metaclassifier.py [-h] [-o OUTPUT_DIR] [-f {paired,single}] [-m]
                         [-c {order,family,genus,species}] [-p MIN_PROPORTION]
                         [-i MAX_MARKERS] [-r PEAR_MERGER]
                         [-s SEQTK_CONVERTER] [-a VSEARCH_ALIGNER]
                         [-t THREADS] [-v]
                         SAMPLE_FILE DB_DIR CONFIG_FILE
                         
The metaclassifier.py uses DNA metabarcoding sequence reads and a database of marker sequences, including
corresponding taxonomy classes to identify and quantify the floral composition of honey

positional arguments:
  SAMPLE_FILE           Input tab-delimited file specifying sample names, file names for forward paired-end
                        reads, and file names for reverse paired-end (file path if not in working directory)
                        The second file not required for single-end frangments
                        
  DB_DIR                Input marker database directory with sequence fasta and corresponding taxonomy lineage
                        files for each marker
                        
  CONFIG_FILE           Input tab-delimited file specifying marker name, and its corresponding VSEARCH's
                        usearch_global function minimum query coverage (i.e. 0.8 for 80%) and minimun sequence
                        identity (i.e. 0.95 for 95%) for each search marker (provide the file path if not in
                        if the VSEARCH settings configuration is not in working directory)
                        

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Specify output directory name, otherwise it will automatically be created using the
                        input sample table file name
                        
  -f {paired,single}, --frag_type {paired,single}
                        Specify the sequence fragment type in the input sample file, available options are:
                        paired: single-end read fragments (default)
                        single: paired-end read fragments
                        
  -m, --merge           Merge overlapping paired-end reads (default: False)
                        
  -c {order,family,genus,species}, --tax_class {order,family,genus,species}
                        Taxonomy class for quantify taxon level marker read abundance (default: genus)
                        
  -p MIN_PROPORTION, --min_proportion MIN_PROPORTION
                        Minimum taxon read proportion allowed to retain a sample taxon, allowed proportion,
                        ranges from 0.00 to 0.01 (default = 0.00)
                        
  -i MAX_MARKERS, --max_markers MAX_MARKERS
                        Maximum missing markers allowed to retain a sample taxon (default = 0)
                        
  -r PEAR_MERGER, --pear_merger PEAR_MERGER
                        Path to PEAR, the paired-end read merger if not in environmental variables (ENV)
                        (default: read from ENV)
                        
  -s SEQTK_CONVERTER, --seqtk_converter SEQTK_CONVERTER
                        Path to seqtk, the sequence processing tool if not in environmental variables (ENV)
                        (default: read from ENV)
                        
  -a VSEARCH_ALIGNER, --vsearch_aligner VSEARCH_ALIGNER
                        Path to VSEARCH, the sequence analysis tool if not in environmental variables (ENV)
                        (default: read from ENV)
                        
  -t THREADS, --threads THREADS
                        Specify the number of threads to use (default: 2)
                        
  -v, --version         Print the current metaclassifier.py version and exit
```

### process_reads.py
```
usage: process_reads.py [-h] [-o OUTPUT_DIR] [-f {paired,single}] [-m]
                        [-r PEAR_MERGER] [-s SEQTK_CONVERTER] [-t THREADS]
                        [-v]
                        SAMPLE_FILE

The process_reads.py script optionally merges overlapping paired-end (PE) reads when the DNA fragment is
shorter than two times the read length, and converts fastq to fasta format

positional arguments:
  SAMPLE_FILE           Input tab-delimited file specifying sample names, file names for forward paired-end
                        reads, and file names for reverse paired-end (file path if not in working directory)
                        The second file not required for single-end frangments
                        

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Specify output directory name, otherwise it will automatically be created using the
                        input sample table file name
                        
  -f {paired,single}, --frag_type {paired,single}
                        Specify the sequence fragment type in the input sample file, available options are:
                        paired: single-end read fragments (default)
                        single: paired-end read fragments
                        
  -m, --merge           Merge overlapping paired-end reads (default: False)
                        
  -r PEAR_MERGER, --pear_merger PEAR_MERGER
                        Path to PEAR, the paired-end read merger if not in environmental variables (ENV)
                        (default: read from ENV)
                        
  -s SEQTK_CONVERTER, --seqtk_converter SEQTK_CONVERTER
                        Path to seqtk, the sequence processing tool if not in environmental variables (ENV)
                        (default: read from ENV)
                        
  -t THREADS, --threads THREADS
                        Specify the number of threads to use (default: 2)
                        
  -v, --version         Print the current process_reads.py version and exit
```

### classify_reads.py
```
usage: classify_reads.py [-h] [-a VSEARCH_ALIGNER]
                         [-c {order,family,genus,species}] [-p MIN_PROPORTION]
                         [-i MAX_MARKERS] [-t THREADS] [-v]
                         FASTA_DIR DB_DIR CONFIG_FILE
 
The classify_reads.py script searches for marker sequences in mult-fasta short reads sample datasets,
identifies sample taxonomic composition, and quantifies taxon read abundance

positional arguments:
  FASTA_DIR             Input sample read data fasta file directory - either merged paired-end PE or
                        single-end fasta files
                        
  DB_DIR                Input marker database directory with sequence fasta and corresponding taxonomy lineage
                        files for each marker
                        
  CONFIG_FILE           Input tab-delimited file specifying marker name, and its corresponding VSEARCH's
                        usearch_global function minimum query coverage (i.e. 0.8 for 80%) and minimun sequence
                        identity (i.e. 0.95 for 95%) for each search marker (provide the file path if not in
                        if the VSEARCH settings configuration is not in working directory)
                        

optional arguments:
  -h, --help            show this help message and exit
  -a VSEARCH_ALIGNER, --vsearch_aligner VSEARCH_ALIGNER
                        Path to VSEARCH, the sequence analysis tool if not in environmental variables (ENV)
                        (default: read from ENV)
                        
  -c {order,family,genus,species}, --tax_class {order,family,genus,species}
                        Taxonomy class for quantify taxon level marker read abundance (default: genus)
                        
  -p MIN_PROPORTION, --min_proportion MIN_PROPORTION
                        Minimum taxon read proportion allowed to retain a sample taxon, allowed proportion,
                        ranges from 0.00 to 0.01 (default = 0.00)
                        
  -i MAX_MARKERS, --max_markers MAX_MARKERS
                        Maximum missing markers allowed to retain a sample taxon (default = 0)
                        
  -t THREADS, --threads THREADS
                        Specify the number of threads to use (default: 2)
                        
  -v, --version         Print the current classify_reads.py version and exit
```

## Citation
If you use MetaClassifier please cite the following paper the describes the methodology:

**Characterizing the floral resources of a North American metropolis using a honey bee foraging assay.**  
_Douglas B. Sponsler_, _Don Shump_,  _Rodney T. Richardson_,  _Christina M. Grozinger_.  
Ecosphere 11, no. 4 (2020): e03102.   
DOI: [https://doi.org/10.1002/ecs2.3102](https://doi.org/10.1002/ecs2.3102)

## License
MetaClassifier is distributed under the GNU GPL v3.0 For more information, see [license](../LICENSE).

