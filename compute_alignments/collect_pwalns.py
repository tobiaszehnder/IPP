#!/usr/bin/env python
#
# This script collects a set of given pairwise alignment chain files,
# converts chromosome names to indices,
# and saves them in binary format readable for the c++ ipp module
import argparse
import os
import sys

import numpy as np
import pandas as pd
import tqdm


# Whenever you change the output format in an incompatible way, be sure to
# increase the FORMAT_VERSION.
FORMAT_VERSION = 4
# Format:
#   version                   [uint8]
#   endianness_magic          [uint16]
#   num_sp1                   [uint8]
#   {
#     sp1_name                [null-terminated string]
#     sp1_genome_size         [uint64]
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
#   num_chomosomes            [uint32]
#   {
#     chrom_name              [null-terminated string]
#   } num_chromosomes times


def create_ptr_to_chromosome(chrom_name):
    # Checks whether the given chromosome is already in `chroms`. If so, then
    # its index is returned. Otherwise, it is added to that list.
    global chroms
    global chroms_map
    if chrom_name in chroms_map:
        return chroms_map[chrom_name]
    else:
        idx = len(chroms)
        chroms_map[chrom_name] = idx
        chroms.append(chrom_name)
        return idx


def read_alignment_blocks_from_chain(chain_file):
    # read chain file and return data frame with alignment blocks

    # initialize progress bar
    pbar = tqdm.tqdm(leave=False)
    pbar.set_description("Parsing %s" % chain_file)

    # read all lines in chain file
    with open(chain_file, "r") as f:
        lines = f.readlines()

    # predefine number of alignment blocks, reserve memory and update progress bar
    n_alignment_blocks = sum(
        [1 for line in lines if not line.startswith(("#", "chain", "\n"))]
    )
    alignments = np.empty(shape=(n_alignment_blocks, 5), dtype=np.uint32)
    pbar.total = n_alignment_blocks
    pbar.refresh()

    # store coordinates of alignment blocks from chain file to data frame
    i = 0
    for line in lines:
        if (line.startswith("#")) or (line == "\n"):
            # Comment line
            continue
        if line.startswith("chain"):
            # Header line. Format:
            # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
            line = line.strip().split()
            ref_chrom, qry_chrom, qry_strand = [line[i] for i in (2, 7, 9)]
            ref_chrom_size, ref_start, ref_end, qry_chrom_size, qry_start, qry_end = [
                int(line[i]) for i in (3, 5, 6, 8, 10, 11)
            ]
            # convert chromosome names to indices and save in chrom_map
            ref_chrom_idx = create_ptr_to_chromosome(ref_chrom)
            qry_chrom_idx = create_ptr_to_chromosome(qry_chrom)
        else:
            # read alignment block
            line = [int(x) for x in line.strip().split()]
            block_width = line[0]
            assert block_width >= 0
            assert block_width < 2**15, "The MSB must not be set by the length"
            length_and_strand = block_width | ((1 << 15) if qry_strand == "-" else 0)
            # The chain file stores coordinates on the reverse strand
            # from the right end, i.e. chrom_size - x.
            # The end coordinate is inclusive. an alignment block of
            # size 3 at the end of a chain [0,5) should have the
            # coordinates [2,4] (i.e. the positions 2,3,4).
            alignments[i] = aln_block = [
                ref_chrom_idx,
                ref_start,
                qry_chrom_idx,
                qry_start if qry_strand == "+" else (qry_chrom_size - qry_start - 1),
                length_and_strand,
            ]
            i += 1
            pbar.update()
            if len(line) == 3:
                # The last line has only length 1
                # Length 3 means that there are more alignment blocks
                # Move the ref and qry start by the sum of the alignment
                # block and the respective gap size
                ref_gap = line[1]
                qry_gap = line[2]
                ref_start += block_width + ref_gap
                qry_start += block_width + qry_gap

    cols = ["ref_chrom", "ref_start", "qry_chrom", "qry_start", "length_and_strand"]

    df = pd.DataFrame(alignments, columns=cols)

    # sort alignment blocks and remove duplicates
    df.sort_values(
        by=["ref_chrom", "ref_start", "qry_chrom", "qry_start"],
        ignore_index=True,
        inplace=True,
    )
    df.drop_duplicates(ignore_index=True, inplace=True)

    pbar.close()

    return df


def main():
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("chain_dir")
    parser.add_argument("assembly_dir")
    parser.add_argument("species_list")
    parser.add_argument("outfile")
    parser.add_argument(
        "-l",
        "--leave_progressbar",
        action="store_true",
        help="do not remove completed progressbar",
    )
    args = parser.parse_args()
    species_list = args.species_list.split(",")

    # initialize chromosome list and index dict
    global chroms
    chroms = []
    global chroms_map
    chroms_map = {}

    # define functions
    def get_chain_path(ref, qry):
        return os.path.join(args.chain_dir, "{}.{}.all.pre.chain".format(ref, qry))

    # Write the binary output file.
    print("Collecting alignment blocks from all species pairs")
    print("Writing output to %s" % args.outfile)
    n_species = len(species_list)
    n_species_pairs = n_species * (n_species - 1)
    pbar = tqdm.tqdm(total=n_species_pairs, leave=args.leave_progressbar)
    with open(args.outfile, "wb") as out:

        def write_int(i, length):
            assert i >= 0 and i < 2 ** (length * 8)
            out.write(int(i).to_bytes(length, byteorder=sys.byteorder))

        def write_str(s):
            out.write(s.encode())
            out.write(b"\x00")

        # Write the version of the output format that is used. That enables the
        # consumer to verify that it does not read any outdated pwalns file.
        write_int(FORMAT_VERSION, 1)

        # Write a magic number to ensure that the endianness of the system that
        # produced the pwalns file is the same as the endianness of the system that
        # consumes it.
        write_int(0xAFFE, 2)

        write_int(n_species, 1)
        for sp1 in species_list:
            write_str(sp1)

            # read and write genome size
            with open(os.path.join(args.assembly_dir, sp1 + ".sizes"), "r") as f:
                genome_size_sp1 = sum(
                    [int(line.strip().split()[1]) for line in f.readlines()]
                )
            write_int(genome_size_sp1, 8)

            write_int(n_species - 1, 1)
            for sp2 in species_list:
                if sp1 == sp2:
                    continue

                # write sp2 name
                write_str(sp2)

                # read chain file and save alignment blocks in df
                chain_file = get_chain_path(sp1, sp2)
                df = read_alignment_blocks_from_chain(chain_file)

                # write number of ref chrom entries
                ref_chrom_counts = df.ref_chrom.value_counts(sort=False)
                write_int(len(ref_chrom_counts), 4)

                # loop through df rows and write to file
                pbar2 = tqdm.tqdm(total=df.shape[0], leave=False)
                pbar2.set_description("Writing alignment blocks (%s, %s)" % (sp1, sp2))
                last_ref_chrom = None
                for row in df.itertuples(index=False):
                    if row.ref_chrom != last_ref_chrom:
                        # This is a new ref_chrom -> write a new "header"
                        write_int(row.ref_chrom, 4)
                        write_int(ref_chrom_counts[row.ref_chrom], 4)
                        last_ref_chrom = row.ref_chrom
                    write_int(row.ref_start, 4)
                    write_int(row.qry_start, 4)
                    write_int(row.qry_chrom, 4)
                    write_int(row.length_and_strand, 2)
                    pbar2.update()

                pbar2.close()

                pbar.update()

        # write chromosomes
        write_int(len(chroms), 4)
        for chrom in chroms:
            write_str(chrom)

    pbar.close()
    return


if __name__ == "__main__":
    main()
