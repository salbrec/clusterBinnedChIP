"""
Utils used to preprocess data or to compute distances.
A feq helper functions used by the main script "cluNbinChIP.py"

date:	January 2023
author:	Steffen Albrecht
"""

import os
import pandas as pd

def get_ID_pair(ID1, ID2):
    key = '_'.join(sorted([ID1.replace(' ','_'),
        ID2.replace(' ','_')]))
    return key

def jaccard_dist_sets(set1, set2):
    inter = len(set1 & set2)
    jaccard_index = inter / (len(set1) + len(set2) - inter)
    return 1.0 - jaccard_index

def jaccard_dist_peaks(file1, file2):
    bedtools = 'bedtools jaccard -a %s -b %s '%(file1, file2)
    bedtools += "| awk '{print $3}' | tail -n1"
    return 1.0 - float(os.popen(bedtools).read().strip())

def get_peaks(file_path, chroms_used):
    peaks = {}
    for chrom in chroms_used:
        peaks[chrom] = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if len(line) < 3:
                raise RuntimeError('The given bed-file has less than 3 columns: '+file_path)
            chrom = line[0]
            if not chrom in chroms_used:
                continue
            peaks[chrom].append((chrom, int(line[1]), int(line[2])))
    return peaks

def get_bins_from_file(bed_file, chr_sizes_file, size):
    bin_size = int(size.replace('bp','').replace('kb','000'))
    chrom_sizes = pd.read_csv(chr_sizes_file, sep='\t', names=['chr', 'size'])
    chroms_used = set(chrom_sizes['chr']) - set(['chrX', 'chrY'])
    chrom_sizes = dict(zip(chrom_sizes['chr'], chrom_sizes['size']))

    peaks = get_peaks(bed_file, chroms_used)

    bins_set = set()
    position_map = {}
    bin_ID = 0
    for chrom in chroms_used:
        chrom_size = chrom_sizes[chrom]
        bin_end = bin_size
        while bin_end < chrom_size:
            bin_start = bin_end - bin_size
            if chrom in peaks:
                for peak in peaks[chrom]:
                    if peak[2] > bin_start:
                        if peak[1] < bin_end:
                            bins_set.add(bin_ID)
                        else:				# early exit: since peaks are sorted the loop can be broken
                            break			#             when the start of the bin is "behind" the end of the current peak
            output_end = bin_end if bin_end < chrom_size else chrom_size
            position_map[bin_ID] = (chrom, bin_start, output_end)
            # update the running values
            bin_end += bin_size
            bin_ID += 1
    return bins_set, bin_ID, position_map

def get_path_to_filename(fn):
    return fn[:-fn[::-1].find('/')]

