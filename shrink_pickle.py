#!/usr/bin/env python
import numpy as np
import pandas as pd
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('pickle_file')
args = parser.parse_args()

with open(args.pickle_file, 'rb') as f:
    pwaln = pickle.load(f)

chroms = []
chroms_map = {}
def create_ptr_to_chromosome(chrom_name):
    # Checks whether the given chromosome is already in `chroms`. If so, then
    # its index is returned. Otherwise, it is added to that list.
    if chrom_name in chroms_map:
        return chroms_map[chrom_name]
    else:
        idx = len(chroms)
        chroms_map[chrom_name] = idx
        chroms.append(chrom_name)
        return idx

# this doesn't need to be in the alignment pipeline, not important for binarized stuff
for k,v in pwaln.items():
    for k2,df in v.items():
        df["ref_chrom"] = df["ref_chrom"].apply(create_ptr_to_chromosome)
        df["qry_chrom"] = df["qry_chrom"].apply(create_ptr_to_chromosome)
        pwaln[k][k2] = df.astype({
            "ref_chrom": np.uint16,
            "ref_start": np.uint32,
            "ref_end": np.uint32,
            "qry_chrom": np.uint16,
            "qry_start": np.uint32,
            "qry_end": np.uint32})

outfile = args.pickle_file.replace('.pkl', '.shrunk.pkl')
with open(outfile, "wb") as f:
    data = {
        "chromosomes": chroms,
        "pwaln": pwaln,
    }
    pickle.dump(data,f)
