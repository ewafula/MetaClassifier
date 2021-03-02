#!/usr/bin/env python3

"""
metaclassifier.classify_reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Functions for classifying and quantifying sample reads into taxonomy groups

A collection of functions to searching for marker sequences in mult-fasta short reads sample
datasets, identifying sample taxonomic composition, and quantifying taxon read abundance.
"""

__author__ = ('Eric Wafula (ewafula@gmail.com)')
__version__ = '1.0.0'
__date__ = '05 March 2021'


import os
import sys
import time
import subprocess
import argparse
import pandas as pd
from functools import reduce


TAX_CLASS_CHOICES = ["order", "family", "genus", "species"]
MIN_PROPORTION = 0.00
MAX_MARKERS = 0
VSEARCH = "vsearch"
THREADS = "2"
  

def read_parameters():
     p = argparse.ArgumentParser(description=("The classify_reads.py script searches for marker sequences in mult-fasta"       
                                             " short reads sample datasets,\nidentifies sample taxonomic composition,"
                                             " and quantifies taxon read abundance"),
                                formatter_class=argparse.RawTextHelpFormatter)
     p.add_argument('FASTA_DIR', type=str, default=None,
                    help="Input sample read data fasta file directory - either merged paired-end PE or\n"
                         "single-end fasta files\n\n")
     p.add_argument('DB_DIR', type=str, default=None,
                    help="Input marker database directory with sequence fasta and corresponding taxonomy lineage\n"
                         "files for each marker\n\n")
     p.add_argument('CONFIG_FILE', type=str, default=None,
                    help="Input tab-delimited file specifying marker name, and its corresponding VSEARCH's\n"
                         "usearch_global function minimum query coverage (i.e. 0.8 for 80%%) and minimun sequence\n"
                         "identity (i.e. 0.95 for 95%%) for each search marker (provide the file path if not in\n"
                         "if the VSEARCH settings configuration is not in working directory)\n\n")                      
     p.add_argument('-a', '--vsearch_aligner', type=str, default=VSEARCH,
                    help="Path to VSEARCH, the sequence analysis tool if not in environmental variables (ENV)\n"
                         "(default: read from ENV)\n\n") 
     p.add_argument('-c', '--tax_class', type=str, default='genus', choices=TAX_CLASS_CHOICES,
                    help="Taxonomy class for quantify taxon level marker read abundance (default: genus)\n\n")
     p.add_argument('-p', '--min_proportion', type=float, default=MIN_PROPORTION,
                    help="Minimum taxon read proportion allowed to retain a sample taxon, allowed proportion,\n"
                         "ranges from 0.00 to 0.01 (default = 0.00)\n\n")
     p.add_argument('-i', '--max_markers', type=int, default=MAX_MARKERS,
                    help="Maximum missing markers allowed to retain a sample taxon (default = 0)\n\n")               
     p.add_argument('-t', '--threads', type=str, default=THREADS,
                    help="Specify the number of threads to use (default: 2)\n\n")
     p.add_argument('-v', '--version', action='version',
                   version=f"classify_reads.py version {__version__} ({__date__})",
                   help="Print the current classify_reads.py version and exit\n\n")
     return p.parse_args()


def search_markers(fasta_dir, db_dir, config_file, aligner, tax_class, min_proportion, max_markers, threads):
     """Searches for marker sequences in mult-fasta short reads sample datasets

     Parameters
     ----------
     fasta_dir : str
          Input sample read data fasta file directory
     db_dir : str 
          Input marker database directory
     config_file : str
          Input tab-delimited file specifying marker name, and its
          corresponding VSEARCH search parameters
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

     print(f"{get_localtime()} - - Searching for markers in the sample read dataset(s)...\n")

     # get sample name and create sample out directory
     for file in os.listdir(fasta_dir):
          sample_fasta = f"{fasta_dir}/{file}"
          sample = os.path.splitext(file)[0]
          print(f"{get_localtime()} - - - Searching sample {sample}...\n")
          fasta_dir_path_components = os.path.abspath(fasta_dir).split(os.sep)
          output_dir = "/".join(fasta_dir_path_components[:-1])
          sample_dir = f"{output_dir}/{sample}"
          os.mkdir(sample_dir,  mode=0o775)
          

          # search for marker with vsearch aligner
          marker_dfs = []  
          with open(config_file, 'r') as config:
               for line in config:
                    params = line.strip().split("\t")
                    marker = params[0]
                    seq_identity = params[1]
                    query_coverage = params[2]
                    marker_fasta = f"{db_dir}/{params[0]}.fa"
                    marker_search = f"{sample_dir}/{params[0]}_blast6out.tsv"
                    log_file = f"{sample_dir}/{params[0]}_vsearch.log"
                    subprocess.run([aligner, '--usearch_global', sample_fasta, '--db', marker_fasta, '--blast6out', marker_search, '--id', seq_identity, '--query_cov', query_coverage, '--top_hits_only', '--log', log_file, '--quiet', '--threads', threads], capture_output=True, check=True)
                   
                    # get taxon read proportions
                    marker_lineage = f"{db_dir}/{marker}.tax"
                    marker_dfs.append(get_read_proportion(marker_search, marker_lineage, tax_class, min_proportion, marker))
                     
               # get rescaled median taxon read proportions and write final sample results (markers and taxa) to file
               sample_df = get_rescaled_median_proportion(marker_dfs, tax_class, max_markers)
               sample_results = f"{sample_dir}/{sample}_rescaled_propotions.tsv"
               sample_df.to_csv(sample_results, sep="\t", index=False, float_format="%.6f", encoding="utf-8")



def get_read_proportion(align_resluts, marker_lineage, tax_class, min_proportion, marker):
     """Identifies sample taxonomic composition and quantifies class read abundance

     Parameters
     ----------
     align_resluts : str
          VSEARCH blast6out sequence search results
     marker_lineage : str 
          the marker lineage table from the input marker database directory
     tax_class : str
          Taxonomy class for quantify taxon level marker read abundance
     min_proportion : float
          Minimum taxon read proportion allowed to retain a sample taxon
     marker : str 
          The sample marker name to search
     
     Returns
     -------
          list
               a list of pandas dataframe objects with marker taxon class read
               abundance proportions
     """

     print(f"{get_localtime()} - - - - Computing {tax_class} taxonomy class read proportions for {marker} marker...\n")
     
     # load marker vsearch results and taxonomy lineage tables
     blast6out_names = ["query", "target", "id", "alnlen", "mism", "opens", "qlo", "qhi", "tlo", "thi", 
                         "evalue", "bits"]
     blast6out = pd.read_csv(align_resluts, names=blast6out_names, sep="\t")
     lineage_names = ["taxon_id", "order", "family", "genus", "species"]
     lineage = pd.read_table(marker_lineage, names=lineage_names, sep="\t") 
     merged_tables = pd.merge(blast6out, lineage.drop_duplicates(subset=['taxon_id'], keep='first'),
                         how="left", left_on="target", right_on="taxon_id")

     # select taxonomy class, quantify taxon read proportions,
     # and filter out taxonomic classes with low proportions
     merged_tables["reads"] = 1
     merged_tables = merged_tables[[tax_class, "reads"]]
     read_counts = merged_tables.groupby(tax_class).count().reset_index()
     read_counts[marker] = read_counts["reads"] / merged_tables["reads"].count()
     read_counts = read_counts[[tax_class, marker]]
     read_counts = read_counts[read_counts[marker] >= min_proportion]
     return read_counts


def get_rescaled_median_proportion(marker_dfs, tax_class, max_markers):
     """Computing and rescaling median proportions for the taxonomy class across all markers

     Parameters
     ----------
     marker_dfs : list
          a list of pandas dataframe objects with marker taxon class read
          abundance proportions
     tax_class : str
          Taxonomy class for quantify taxon level marker read abundance
     max_markers : int
          Maximum missing markers allowed to retain a sample taxon
     
     Returns
     -------
          pandas dataframe object
               a pandas dataframes object with marker taxon class read
               abundance proportions, median proportions, and rescaled
               median proportions
     """

     print(f"{get_localtime()} - - - - Computing and rescaling median proportions for {tax_class} taxonomy class across all markers...\n")

     # compute rescaled median taxon reads proportions
     sample_df = reduce(lambda x, y: pd.merge(x, y, how="outer", on=[tax_class]), marker_dfs)
     max_markers += 1
     sample_df = sample_df[sample_df.isnull().sum(axis=1) < max_markers].fillna(0)
     sample_df["median"] = sample_df.iloc[:,1:].median(axis=1)
     sample_df["proportion"] = sample_df["median"] / sample_df["median"].sum()
     return sample_df


def main():
     t0 = time.time()
     args = read_parameters()

     print(f"{get_localtime()} - Starting read classification...\n")
     if args.min_proportion >= 0.01:
          raise Exception(f"The minimum taxon read proportion allowed to retain a sample taxon cannot exceed 0.1%!")
     search_markers(args.FASTA_DIR, args.DB_DIR, args.CONFIG_FILE, args.vsearch_aligner, args.tax_class, args.min_proportion, args.max_markers, args.threads)
 
     t1 = time.time()
     print(f"{get_localtime()} - Completed read classification...")
     print(f"Total elapsed time {int(t1 - t0)}\n")
     sys.exit(0)


def get_localtime():
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    return current_time


if __name__ == '__main__':
    main()