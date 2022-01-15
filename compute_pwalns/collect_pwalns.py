#!/usr/bin/env python

# This script collects a set of given pairwise alignment tables (.tbl format)
# and saves them in binary format readable for the c++ ipp module
#
# Format:
#   num_chomosomes            [uint16]
#   {
#     chrom_name              [null-terminated string]
#   } num_chromosomes times
#   num_sp1                   [uint8]
#   {
#     sp1_name                [null-terminated string]
#     num_sp2                 [uint8]
#     {
#       sp2_name              [null-terminated string]
#       num_ref_chrom_entries [uint32]
#       {
#         num_pwaln_entries   [uint32]
#         {
#           ref_start         [uint32]
#           ref_end           [uint32]
#           qry_start         [uint32]
#           qry_end           [uint32]
#           ref_chrom         [uint16]
#           qry_chrom         [uint16]
#         } num_pwaln_entries times
#       } num_ref_chrom_entries times
#     } num_sp2 times
#   } num_sp1 times

import argparse
import os
import sys
import tqdm
import pandas as pd

# parse arguments
parser = argparse.ArgumentParser()
args, unknownargs = parser.parse_known_args()
parser.add_argument('tbl_dir')
parser.add_argument('species_list')
parser.add_argument('outfile')
args = parser.parse_args()

# read pairwise alignment tbl files
cols = ['ref_chrom','ref_start','ref_end','qry_chrom','qry_start','qry_end']
species_list = args.species_list.split(',')
def get_tbl_path(ref,qry):
  return os.path.join(args.tbl_dir, '{}.{}.tbl'.format(ref, qry))
pwalns = {ref: {qry: pd.read_csv(get_tbl_path(ref,qry), sep='\t', header=None, names=cols) for qry in species_list if not ref == qry} for ref in species_list}

# change chromosome names to indices
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
for ref,v in pwalns.items():
  for qry,df in v.items():
    df['ref_chrom'].apply(create_ptr_to_chromosome)
    df['qry_chrom'].apply(create_ptr_to_chromosome)

# Compute the total number of rows to complete.
print("Sorting and removing duplicates from the pwalns")
total_pwalns = 0
for _, v in pwalns.items():
    total_pwalns += len(v)
    
pbar = tqdm.tqdm(total=total_pwalns, leave=False)
total_rows = 0
for _, v in pwalns.items():
  for _, df in v.items():
    # Sort the df and compute how many rows there are with the same
    # ref_chrom value.
    df.sort_values(by=["ref_chrom", "ref_start", "qry_chrom", "qry_start"],
                   ignore_index=True,
                   inplace=True)
    df.drop_duplicates(ignore_index=True, inplace=True)
    total_rows += df.shape[0]
    pbar.update()
pbar.close()

# Write the binary output file.
print("Writing the output file")
pbar = tqdm.tqdm(total=total_rows, leave=False)
with open(args.outfile, "wb") as out:
  def write_int(i, length):
    out.write(int(i).to_bytes(length, byteorder=sys.byteorder))
    
  def write_str(s):
    out.write(s.encode())
    out.write(b'\x00')
    
  write_int(len(chroms), 2)
  for chrom in chroms:
    write_str(chrom)

  write_int(len(pwalns), 1)
  for sp1, v in pwalns.items():
    write_str(sp1)
    write_int(len(v), 1)
    for sp2, df in v.items():
      ref_chrom_counts = df.ref_chrom.value_counts(sort=False)

      write_str(sp2)
      write_int(len(ref_chrom_counts), 4)

      last_ref_chrom = None
      for row in df.itertuples(index=False):
        if row.ref_chrom != last_ref_chrom:
          # This is a new ref_chrom -> write a new "header"
          write_int(ref_chrom_counts[row.ref_chrom], 4)
          last_ref_chrom = row.ref_chrom
        write_int(row.ref_start, 4)
        write_int(row.ref_end, 4)
        write_int(row.qry_start, 4)
        write_int(row.qry_end, 4)
        write_int(row.ref_chrom, 2)
        write_int(row.qry_chrom, 2)
        pbar.update()
pbar.close()
