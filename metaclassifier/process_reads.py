#!/usr/bin/env python3

"""
metaclassifier.process_reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Functions for merging paired-end (PE) and converting FASTQ to FASTA

A collection of functions for merging overlapping PE reads when the DNA fragment is shorter
than two times the read length, and converting the resulting merged FASTQ to FASTA format.
"""

__author__ = ('Eric Wafula (ewafula@gmail.com)')
__version__ = '1.0.0'
__date__ = '05 March 2021'


import os
import sys
import time
import subprocess
import argparse


OUTPUT_DIR = ""
FRAGMENT_CHOICES = ['paired', 'single']
PEAR = "pear"
SEQTK = "seqtk"
THREADS = "2"
  

def read_parameters():
    p = argparse.ArgumentParser(description=("The process_reads.py script optionally merges overlapping"        
                                             " paired-end (PE) reads when the DNA fragment is\nshorter than"
                                             " two times the read length, and converts fastq to fasta format"),
                                formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('SAMPLE_FILE', type=str, default=None,
                    help="Input tab-delimited file specifying sample names, file names for forward paired-end\n"
                         "reads, and file names for reverse paired-end (file path if not in working directory)\n"
                         "The second file not required for single-end frangments\n\n")
    p.add_argument('-o', '--output_dir', type=str, default=None,
                    help="Specify output directory name, otherwise it will automatically be created using the\n"
                         "input sample table file name\n\n")   
    p.add_argument('-f', '--frag_type', type=str, default='paired', choices=FRAGMENT_CHOICES,
                    help="Specify the sequence fragment type in the input sample file, available options are:\n"
                         "paired: single-end read fragments (default)\n"
                         "single: paired-end read fragments\n\n")
    p.add_argument('-m', '--merge', action='store_true',
                    help="Merge overlapping paired-end reads (default: False)\n\n")
    p.add_argument('-r', '--pear_merger', type=str, default=PEAR,
                    help="Path to PEAR, the paired-end read merger if not in environmental variables (ENV)\n"
                         "(default: read from ENV)\n\n")
    p.add_argument('-s', '--seqtk_converter', type=str, default=SEQTK,
                    help="Path to seqtk, the sequence processing tool if not in environmental variables (ENV)\n"
                         "(default: read from ENV)\n\n")                
    p.add_argument('-t', '--threads', type=str, default=THREADS,
                    help="Specify the number of threads to use (default: 2)\n\n")
    p.add_argument('-v', '--version', action='version',
                   version=f"process_reads.py version {__version__} ({__date__})",
                   help="Print the current process_reads.py version and exit\n\n")
    return p.parse_args()


def merge_pairs(sample_input, output_dir, merger, threads):
    """Merges overlapping paired-end (PE) sample FASTQ reads

    Parameters
    ----------
    sample_input : str
        Input tab-delimited sample file specifying sample name
        and PE FASTQ files
    output_dir : str
        output directory
    merger : str
        The PEAR, the PE read merger executable
    threads : str 
        The number of threads for PEAR, the PE read merger

    Returns
    -------
        None
    """

    print(f"{get_localtime()} - - Merging paired-end (PE) sample read dataset(s)...\n")

    # get paired-end sample read file and merge with PEAR
    merge_dir = f"{output_dir}/merge_dir"
    os.mkdir(merge_dir,  mode=0o775)
    sample_output = f"{output_dir}/sample.tsv"
    with open(sample_input, 'r') as input:
        with open(sample_output, 'w') as output:
            for sample in input:  
                fields = sample.strip().split("\t")
                if len(fields) != 3:
                    raise Exception(f"The fields (columns) in input sample file, {sample_input} are not the required three: sample name, forward reads, and reverse reads")
                if os.path.exists(fields[1]) and os.path.exists(fields[2]):
                    merged_fastq = f"{merge_dir}/{fields[0]}"
                    log_file = f"{merge_dir}/{fields[0]}_pear.log"
                    with open(log_file, 'w') as log:
                        subprocess.run([merger, '-f', fields[1], '-r', fields[2], '-o', merged_fastq, '-j', threads], stdout=log, check=True)
                    output.write(f"{fields[0]}\t{merge_dir}/{fields[0]}.assembled.fastq\n")
                else:
                    raise Exception(f"file {fields[1]} and {fields[2]} for sample {fields[0]} do not exist")


def get_fasta(sample_input, output_dir, converter):
    """Converts merged paired-end (PE) reads from FASTQ to FASTA 

    Parameters
    ----------
    sample_input : str
        Input tab-delimited sample file specifying sample name
        and PE FASTQ files
    output_dir : str
        output directory
    converter : str
        The seqtk, the sequence processing tool executable

    Returns
    -------
        None
    """
    print(f"{get_localtime()} - - Converting merged paired-end (PE) sample read dataset(s) from FASTQ to FASTA...\n")

    # convert merged paired-end fastq sample reads files to fasta format
    fasta_dir = f"{output_dir}/fasta_dir"
    os.mkdir(fasta_dir,  mode=0o775)
    with open(sample_input, 'r') as input:
        for sample in input: 
            fields = sample.strip().split("\t")
            if os.path.exists(fields[1]):
                fasta_file = f"{fasta_dir}/{fields[0]}.fasta"
                with open(fasta_file, 'w') as fasta:
                    subprocess.run([converter, 'seq', '-A', fields[1]], stdout=fasta, check=True)


def main():
    t0 = time.time()
    args = read_parameters()

    print(f"{get_localtime()} - Starting read processing...\n")
    if args.frag_type == "single" and args.merge == True:
        raise Exception(f"single-end read fragments cannot be merged!")
    if args.frag_type == "paired" and args.merge == False:
        raise Exception(f"Overlapping paired-end read fragments need to be merged!") 

    if not args.output_dir:
        OUTPUT_DIR = os.path.splitext(os.path.basename(args.SAMPLE_FILE))[0]
        os.mkdir(OUTPUT_DIR,  mode=0o775)
    else:
        OUTPUT_DIR = args.output_dir

    if args.merge == True:
        merge_pairs(args.SAMPLE_FILE, OUTPUT_DIR, args.pear_merger, args.threads)
        sample_input = f"{OUTPUT_DIR}/sample.tsv"
        get_fasta(sample_input, OUTPUT_DIR, args.seqtk_converter)
    else:
        get_fasta(args.SAMPLE_FILE, OUTPUT_DIR, args.seqtk_converter) 

    t1 = time.time()
    print(f"{get_localtime()} - Completed read processing...")
    print(f"Total elapsed time {int(t1 - t0)}\n")
    sys.exit(0)


def get_localtime():
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    return current_time


if __name__ == '__main__':
    main()