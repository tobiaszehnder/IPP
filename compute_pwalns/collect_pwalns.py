#!/usr/bin/env python
#
# This script collects a set of given pairwise alignment tables (.tbl format),
# converts chromosome names to indices,
# and saves them in binary format readable for the c++ ipp module
import argparse
import os
import sys
import tqdm
import pandas as pd

# Whenever you change the output format in an incompatible way, be sure to
# increase the FORMAT_VERSION.
FORMAT_VERSION = 2
# Format:
#   version                   [uint8]
#   endianness_magic          [uint16]
#   num_chomosomes            [uint32]
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
#         ref_chrom           [uint32]
#         num_pwaln_entries   [uint32]
#         {
#           ref_start         [uint32]
#           qry_start         [uint32]
#           qry_chrom         [uint32]
#           length_and_strand [uint16]
#         } num_pwaln_entries times
#       } num_ref_chrom_entries times
#     } num_sp2 times
#   } num_sp1 times

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('tbl_dir')
parser.add_argument('species_list')
parser.add_argument('outfile')
args = parser.parse_args()

# read pairwise alignment tbl files
cols = ['ref_chrom','ref_start','ref_end','qry_chrom','qry_start','qry_end']
species_list = args.species_list.split(',')
def get_tbl_path(ref,qry):
  return os.path.join(args.tbl_dir, '{}.{}.tbl'.format(ref, qry))
print("Reading pairwise alignment tables")
total_files = len(species_list) * (len(species_list) - 1)
pbar = tqdm.tqdm(total=total_files, leave=False)
pwalns = {}
for ref in species_list:
  if not ref in pwalns.keys():
    pwalns[ref] = {}
  for qry in species_list:
    if ref == qry:
      continue
    pwalns[ref][qry] = pd.read_csv(get_tbl_path(ref,qry), sep='\t', header=None, names=cols)
    pbar.update()

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
pbar = tqdm.tqdm(total=total_files, leave=False)
print("Converting chromosome names to indices")
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
    assert i >= 0 and i < 2**(length*8)
    out.write(int(i).to_bytes(length, byteorder=sys.byteorder))
    
  def write_str(s):
    out.write(s.encode())
    out.write(b'\x00')

  # Write the version of the output format that is used. That enables the
  # consumer to verify that it does not read any outdated pwalns file.
  write_int(FORMAT_VERSION, 1)

  # Write a magic number to ensure that the endianness of the system that
  # produced the pwalns file is the same as the endianness of the system that
  # consumes it.
  write_int(0xAFFE, 2)
  
  write_int(len(chroms), 4)
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
          write_int(row.ref_chrom, 4)
          write_int(ref_chrom_counts[row.ref_chrom], 4)
          last_ref_chrom = row.ref_chrom

        # MSB of length_and_strand is 1 if the qry is on the negative strand.
        assert row.ref_end >= row.ref_start
        length = row.ref_end - row.ref_start + 1
        assert length < 2**15, "The MSB must not be set by the length"
        length_and_strand = length \
                | ((1<<15) if row.qry_start > row.qry_end else 0)

        write_int(row.ref_start, 4)
        write_int(row.qry_start, 4)
        write_int(row.qry_chrom, 4)
        write_int(length_and_strand, 2)
        pbar.update()
pbar.close()
