#!/usr/bin/env python3

"""
metaclassifier.__main___
~~~~~~~~~~~~~~~~~~~~~~~~

Function quantify the taxonomic floral composition of honey samples read datasets

A wrapper module the "process_reads" and "classify_reads" modules "for merging paired-end (PE)
and converting FASTQ to FASTA" and "for classifying and quantifying sample reads into taxonomy
groups" respectively.
"""

__author__ = ('Eric Wafula (ewafula@gmail.com)')
__version__ = '1.0.0'
__date__ = '05 March 2021'


import os
import sys
import time
import argparse
from . import process_reads, classify_reads


OUTPUT_DIR = ""
FRAGMENT_CHOICES = ['paired', 'single']
TAX_CLASS_CHOICES = ["order", "family", "genus", "species"]
MIN_PROPORTION = 0.00
MAX_MARKERS = 0
PEAR = "pear"
SEQTK = "seqtk"
VSEARCH = "vsearch"
THREADS = "2"


def read_parameters():
    p = argparse.ArgumentParser(description=("The metaclassifier.py uses DNA metabarcoding sequence reads and"        
                                             " a database of marker sequences, including\ncorresponding taxonomy"
                                             " classes to identify and quantify the floral composition of honey"),
                                formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('SAMPLE_FILE', type=str, default=None,
                    help="Input tab-delimited file specifying sample names, file names for forward paired-end\n"
                         "reads, and file names for reverse paired-end (file path if not in working directory)\n"
                         "The second file not required for single-end frangments\n\n")
    p.add_argument('DB_DIR', type=str, default=None,
                    help="Input marker database directory with sequence fasta and corresponding taxonomy lineage\n"
                         "files for each marker\n\n")
    p.add_argument('CONFIG_FILE', type=str, default=None,
                    help="Input tab-delimited file specifying marker name, and its corresponding VSEARCH's\n"
                         "usearch_global function minimum query coverage (i.e. 0.8 for 80%%) and minimun sequence\n"
                         "identity (i.e. 0.95 for 95%%) for each search marker (provide the file path if not in\n"
                         "if the VSEARCH settings configuration is not in working directory)\n\n")
    p.add_argument('-o', '--output_dir', type=str, default=None,
                    help="Specify output directory name, otherwise it will automatically be created using the\n"
                         "input sample table file name\n\n")   
    p.add_argument('-f', '--frag_type', type=str, default='paired', choices=FRAGMENT_CHOICES,
                    help="Specify the sequence fragment type in the input sample file, available options are:\n"
                         "paired: single-end read fragments (default)\n"
                         "single: paired-end read fragments\n\n")
    p.add_argument('-m', '--merge', action='store_true',
                    help="Merge overlapping paired-end reads (default: False)\n\n")
    p.add_argument('-c', '--tax_class', type=str, default='genus', choices=TAX_CLASS_CHOICES,
                    help="Taxonomy class for quantify taxon level marker read abundance (default: genus)\n\n")
    p.add_argument('-p', '--min_proportion', type=float, default=MIN_PROPORTION,
                    help="Minimum taxon read proportion allowed to retain a sample taxon, allowed proportion,\n"
                         "ranges from 0.00 to 0.01 (default = 0.00)\n\n")
    p.add_argument('-i', '--max_markers', type=int, default=MAX_MARKERS,
                    help="Maximum missing markers allowed to retain a sample taxon (default = 0)\n\n")
    p.add_argument('-r', '--pear_merger', type=str, default=PEAR,
                    help="Path to PEAR, the paired-end read merger if not in environmental variables (ENV)\n"
                         "(default: read from ENV)\n\n")
    p.add_argument('-s', '--seqtk_converter', type=str, default=SEQTK,
                    help="Path to seqtk, the sequence processing tool if not in environmental variables (ENV)\n"
                         "(default: read from ENV)\n\n")
    p.add_argument('-a', '--vsearch_aligner', type=str, default=VSEARCH,
                    help="Path to VSEARCH, the sequence analysis tool if not in environmental variables (ENV)\n"
                         "(default: read from ENV)\n\n")                 
    p.add_argument('-t', '--threads', type=str, default=THREADS,
                    help="Specify the number of threads to use (default: 2)\n\n")
    p.add_argument('-v', '--version', action='version',
                   version=f"metaclassifier.py version {__version__} ({__date__})",
                   help="Print the current metaclassifier.py version and exit\n\n")
    return p.parse_args()

def run_metaclassifier(sample_input, db_dir, config_file, output_dir, frag_type, merge, merger,  
                    converter, aligner, tax_class, min_proportion, max_markers, threads):
    """Merges overlapping paired-end (PE) sample FASTQ reads

    Parameters
    ----------
    sample_input : str
        Input tab-delimited sample file specifying sample name
        and PE FASTQ files
    db_dir : str 
        Input marker database directory
    config_file : str
        Input tab-delimited file specifying marker name, and its
        corresponding VSEARCH search parameters
    output_dir : str
        output directory
    frag_type : str
        The sequence fragment type in the input sample file
        (either "paired" or "single")
    merge : boolean
        Either True or False
    merger : str
        The PEAR, the PE read merger executable
    converter : str
        The seqtk, the sequence processing tool executable
    aligner : str
        VSEARCH, the sequence analysis tool executable
    tax_class : str
        Taxonomy class for quantify taxon level marker read abundance
    min_proportion : float
        Minimum taxon read proportion allowed to retain a sample taxon
    max_markers = int
        Maximum missing markers allowed to retain a sample taxon    
    threads : str 
        The number of threads for PEAR, the PE read merger

    Returns
    -------
        None
    """

    print(f"\n===========================================================")
    print(f"MetaClassifier version {__version__} (release date: {__date__})")
    print(f"===========================================================\n")
    print(f"{get_localtime()} - Starting MetaClassifier...\n")

    if frag_type == "single" and merge == True:
        raise Exception(f"single-end read frangments cannot be merged!")

    if frag_type == "paired" and merge == False:
        raise Exception(f"Overlapping paired-end read frangments need to be merged!")

    if min_proportion >= 0.01:
        raise Exception(f"The minimum taxon read proportion allowed to retain a sample taxon cannot exceed 0.1%!")

    if merge == True:
        process_reads.merge_pairs(sample_input, output_dir, merger, threads)
        sample_input = f"{output_dir}/sample.tsv"
        process_reads.get_fasta(sample_input, output_dir, converter)
    else:
        process_reads.get_fasta(sample_input, output_dir, converter)

    fasta_dir = f"{output_dir}/fasta_dir"
    classify_reads.search_markers(fasta_dir, db_dir, config_file, aligner, tax_class, min_proportion, max_markers, threads)

    print(f"{get_localtime()} - Completed MetaClassifier...")


def main():
    t0 = time.time()
    
    args = read_parameters()
    if not args.output_dir:
        OUTPUT_DIR = os.path.splitext(os.path.basename(args.SAMPLE_FILE))[0]
        os.mkdir(OUTPUT_DIR,  mode=0o775)
    else:
        OUTPUT_DIR = args.output_dir
    
    run_metaclassifier(args.SAMPLE_FILE, args.DB_DIR, args.CONFIG_FILE, OUTPUT_DIR, args.frag_type, args.merge, args.pear_merger, 
                    args.seqtk_converter, args.vsearch_aligner, args.tax_class, args.min_proportion, args.max_markers, args.threads)

    t1 = time.time()
    print(f"Total elapsed time {int(t1 - t0)}\n")
    sys.exit(0)


def get_localtime():
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    return current_time


if __name__ == '__main__':
    main()
    