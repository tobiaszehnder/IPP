#!/usr/bin/env python

import argparse
import ipp
import numpy as np
import os
import pandas as pd
import tabulate
import tqdm
import pyranges as pr
from functions import *

LOG_LEVEL_QUIET = 0
LOG_LEVEL_LOG = 1
LOG_LEVEL_DEBUG = 2

log_level = LOG_LEVEL_LOG

def is_log():
    return log_level >= LOG_LEVEL_LOG

def is_debug():
    return log_level >= LOG_LEVEL_DEBUG

def log(*args, **kwargs):
    if is_log():
        print(*args, **kwargs)
def debug(*args, **kwargs):
    if is_debug():
        print(*args, **kwargs)

def debug_shortest_path(shortest_path, simple):
    # Prints the given shortest path.
    # The "out anchors" are the ref coordinates of the anchors of the next
    # shortest path entry.
    # The "in anchors"  are the qry coordinates of the anchors of the current
    # shortest path entry.
    # If the up and down anchors are the same, then only one of them is printed.
    # If simple is given, then the coords are shifted towards zero to make them
    # more readable.
    if not is_debug() or not len(shortest_path):
        return

    headers = ['species', 'score', 'chrom', 'in anchors (strand)', 'loc',
               'out anchors (strand)']
    data = []
    for i in range(len(shortest_path)):
        e = shortest_path[i]

        def anchors_to_ints(anchors):
            # Converts the ("upStart:upEnd", "downStart:downEnd") anchor tuple
            # to an array [upStart, upEnd, downStart, downEnd].
            return [int(i) for anchor in anchors for i in anchor.split(":")]

        not_first = i > 0
        not_last = i < (len(shortest_path) - 1)
        in_anchors = anchors_to_ints(e[4]) if not_first else []
        out_anchors = anchors_to_ints(shortest_path[i+1][3]) if not_last else []

        min_coord = min(in_anchors + out_anchors) if simple else 0
        in_anchors = [i - min_coord for i in in_anchors]
        out_anchors = [i - min_coord for i in out_anchors]

        def anchors_str(anchors):
            # Returns the given anchors pair.
            if not len(anchors):
                # First in or last out
                return ""

            up_str = "{}-{}".format(anchors[0], anchors[1])
            down_str = "{}-{}".format(anchors[2], anchors[3])
            return "{} ({})".format(
                    up_str if up_str == down_str else ", ".join([up_str, down_str]),
                    "+" if anchors[0] < anchors[1] else "-")

        in_anchors_str = anchors_str(in_anchors)
        out_anchors_str = anchors_str(out_anchors)

        chrom, loc = e[2].split(":")
        loc = int(loc) - min_coord

        data.append([e[0], e[1], chrom, in_anchors_str, loc, out_anchors_str])

    debug(tabulate.tabulate(data, headers=headers, stralign="center"))
    debug()


def main():
    parser = argparse.ArgumentParser(description='Independent Point Projections (IPP).\nA method for projecting genomic point coordinates between genomes with large evolutionary distances.')
    parser.add_argument('regions_file', help='Bed file containing genomic coordinates. regions with width > 1 will be centered.')
    parser.add_argument('ref', help='reference species')
    parser.add_argument('qry', help='query species')
    parser.add_argument('path_pwaln')
    parser.add_argument('--out_dir', default=os.getcwd(), help='directory for output files')
    parser.add_argument('--half_life_distance', type=int, default=10000, help='distance to closest anchor point at which projection score is 0.5')
    parser.add_argument('--n_cores', type=int, default=1, help='number of CPUs')
    parser.add_argument('--quiet', action="store_true", help='do not produce any log output')
    parser.add_argument('-v', '--verbose', action="store_true", help='produce additional debugging output')
    parser.add_argument('-s', '--simple_coords', action="store_true", help='make coord numbers in debug output as small as possible')
    parser.add_argument('--data_dir', default='/project/ipp-data')
    args = parser.parse_args()

    global log_level
    if args.quiet:
        log_level = LOG_LEVEL_QUIET
    elif args.verbose:
        log_level = LOG_LEVEL_DEBUG


    # define variables and create output directory
    assembly_dir = args.data_dir + '/assembly/'
    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)
    regions_file_basename = os.path.splitext(os.path.basename(args.regions_file))[0]

    #input("about to init ipp")
    log("Loading pwaln file")
    myIpp = ipp.Ipp()
    myIpp.load_pwalns(args.path_pwaln)
    myIpp.load_genome_sizes(assembly_dir);
    myIpp.set_half_life_distance(args.half_life_distance)

    #input('Press enter to start')
    log('Projecting regions from %s to %s' %(args.ref, args.qry))

    # Read the regions file and enqueue one projection job per line.
    ref_coords = []
    ids = {}
    with open(args.regions_file) as regions_file:
        for i,line in enumerate(regions_file.readlines()):
            cols = line.split('\t')
            # define the coordinate as the center point between start and end coordinates
            name = cols[3]
            refChrom = cols[0]
            refLoc = int(np.mean([int(cols[1]), int(cols[2])]))
            ref_coords.append((refChrom, refLoc))
            # add the name of the region to a dict with refChrom:refLoc as the key for later translation
            ids['{}:{}'.format(refChrom,refLoc)] = (i, name)

    pbar = None
    if not is_debug():
        pbar = tqdm.tqdm(total=len(ref_coords), leave=False)

    results = pd.DataFrame(
        columns=['id', 'coords_ref', 'coords_direct', 'coords_multi',
                 'score_direct', 'score_multi',
                 'ref_anchor_direct_left', 'ref_anchor_direct_right', 'ref_anchor_multi_left', 'ref_anchor_multi_right',
                 'qry_anchor_direct_left', 'qry_anchor_direct_right', 'qry_anchor_multi_left', 'qry_anchor_multi_right',
                 'bridging_species'])
    unmapped = []
    
    def on_job_done_callback(ref_coord,
                             direct_score,
                             direct_coords,
                             direct_ref_anchors,
                             direct_qry_anchors,
                             multi_shortest_path):
        # ref_coord:           string           "ref_chrom:ref_loc"
        # direct_score:        float            [0-1]
        # direct_coords:       string           "qry_chrom:qry_loc"
        # direct_ref_anchors:  (string, string) ("upRefStart:upRefEnd", "downRefStart:downRefEnd")
        # direct_qry_anchors:  (string, string) ("upQryStart:upQryEnd", "downQryStart:downQryEnd")
        # multi_shortest_path: [ref, intermediate1, intermediate2, ..., qry]
        #   Each entry:
        #     species      string           "spX"
        #     score        float            [0-1]
        #     coords       string           "chrom:loc"
        #     ref_anchors  (string, string) ("upRefStart:upRefEnd", "downRefStart:downRefEnd")
        #     qry_anchors  (string, string) ("upQryStart:upQryEnd", "downQryStart:downQryEnd")
        nonlocal results, unmapped, ids
        debug()
        debug("({})".format(len(results)))
        debug_shortest_path(multi_shortest_path, args.simple_coords)
  
        multi_last_entry = multi_shortest_path[-1]

        # handle unmapped region (multi_last_entry is not args.qry)
        if not multi_last_entry[0] == args.qry:
            print(ref_coord)
            print(coord_name_dict[ref_coord])
            unmapped += coord_name_dict[ref_coord]
            return

        assert multi_last_entry[0] == args.qry
        multi_score = multi_last_entry[1]
        multi_coords = multi_last_entry[2]

        # Ref anchors of the first species in the path (first non-reference species)
        multi_ref_anchors = multi_shortest_path[1][3]
  
        # Qry anchors of the last species in the path.
        multi_qry_anchors = multi_last_entry[4]
  
        multi_bridging_species = \
            ','.join([spe[0] for spe in multi_shortest_path[1:-1]])
  
        values = [ids[ref_coord][1], ref_coord, direct_coords, multi_coords,
                  direct_score, multi_score,
                  direct_ref_anchors[0], direct_ref_anchors[1], multi_ref_anchors[0], multi_ref_anchors[1],
                  direct_qry_anchors[0], direct_qry_anchors[1], multi_qry_anchors[0], multi_qry_anchors[1],
                  multi_bridging_species]
        idx = [ids[ref_coord][0]]
        results = results.append(
                pd.DataFrame([values], columns=results.columns, index=idx),
                ignore_index=False)
  
        if pbar:
            pbar.update()
  
    # Start the projection
    myIpp.project_coords(args.ref,
                         args.qry,
                         ref_coords,
                         args.n_cores,
                         on_job_done_callback)
    if pbar:
        pbar.close()
  
    debug()
    debug(results.to_string())

    # write results table to file
    print(results)
    columns_out=['coords_ref', 'coords_direct', 'coords_multi',
                 'score_direct', 'score_multi', 'bridging_species']
    results = results.sort_index().set_index('id').loc[:,columns_out]
    outfile_table = regions_file_basename + '.proj'
    results.to_csv(os.path.join(args.out_dir, outfile_table), sep='\t', header=True, float_format='%.3f')

    # write list of IDs of unmapped regions to file
    outfile_unmapped = regions_file_basename + '.unmapped'
    with open(outfile_unmapped, 'w') as f:
        f.write('\n'.join(unmapped) + '\n')

    # classify projections according to conservation of sequence (DC/IC/NC) and function (+/-)
    thresh = .95
    maxgap = 500
    target_regions = None
    if hasattr(args, 'target_bedfile'):
        target_regions = pr.read_bed(args.target_bedfile)
    results['conservation'] = classify_conservation(results, target_regions, thresh, maxgap)
      
    # write results to bed files (bed9 format for both reference and target species)
    # name column (4):    string    "id_qryChrom:qryLoc_class"    e.g. "peak_75_chr3:3880_IC+"
    outfile_bed_ref = '{}.{}.bed'.format(os.path.join(args.out_dir, regions_file_basename), args.ref)
    outfile_bed_qry = '{}.{}.bed'.format(os.path.join(args.out_dir, regions_file_basename), args.qry)
    bed_ref = results.apply(lambda row: format_row_table_to_bed(row, 'ref'), axis=1)
    bed_qry = results.apply(lambda row: format_row_table_to_bed(row, 'qry'), axis=1)
    bed_ref.to_csv(outfile_bed_ref, sep='\t', index=False, header=None, float_format='%.3f')
    bed_qry.to_csv(outfile_bed_qry, sep='\t', index=False, header=None, float_format='%.3f')
    
if __name__ == '__main__':
    main()

# vim: tabstop=4 shiftwidth=4 expandtab
