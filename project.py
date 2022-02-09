#!/usr/bin/env python

import argparse
import ipp
import numpy as np
import os
import pandas as pd
import tabulate
import tqdm
import pyranges as pr
import sys

LOG_LEVEL_QUIET = 0
LOG_LEVEL_LOG = 1
LOG_LEVEL_DEBUG = 2

log_level = LOG_LEVEL_LOG

def is_log():
    return log_level >= LOG_LEVEL_LOG

def is_debug():
    return log_level >= LOG_LEVEL_DEBUG

pbar = None
def do_print(*args, **kwargs):
    global pbar
    if pbar:
        pbar.clear()
    print(*args, **kwargs)
    if pbar:
        pbar.refresh()

def log(*args, **kwargs):
    if is_log():
        do_print(*args, **kwargs)
def debug(*args, **kwargs):
    if is_debug():
        do_print(*args, **kwargs)

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
        spe = shortest_path[i]

        not_first = i > 0
        not_last = i < (len(shortest_path) - 1)
        nxt_spe = shortest_path[i+1] if not_last else None
        in_anchors = [spe.up_anchor.qry_start, spe.up_anchor.qry_end,
                       spe.down_anchor.qry_start, spe.down_anchor.qry_end] if not_first else []
        out_anchors = [nxt_spe.up_anchor.ref_start, nxt_spe.up_anchor.ref_end,
                       nxt_spe.down_anchor.ref_start, nxt_spe.down_anchor.ref_end] if not_last else []

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

        data.append([spe.species,
                     spe.score,
                     spe.coords.chrom,
                     in_anchors_str,
                     spe.coords.loc - min_coord,
                     out_anchors_str])

    debug(tabulate.tabulate(data, headers=headers, stralign="center"))
    debug()


def main():
    parser = argparse.ArgumentParser(description='Independent Point Projections (IPP).\nA method for projecting genomic point coordinates between genomes with large evolutionary distances.')
    parser.add_argument('regions_file', help='Bed file containing genomic coordinates. regions with width > 1 will be centered.')
    parser.add_argument('ref', help='reference species')
    parser.add_argument('qry', help='query species')
    parser.add_argument('path_pwaln')
    parser.add_argument('-o', '--out_dir', default=os.getcwd(), help='directory for output files')
    parser.add_argument('-n', '--n_cores', type=int, default=1, help='number of CPUs')
    parser.add_argument('-t', '--target_bedfile', default=None, help='functional regions in target species to check for overlap with projections for classification')
    parser.add_argument('-d', '--half_life_distance', type=int, default=10000, help='distance to closest anchor point at which projection score is 0.5')
    parser.add_argument('-q', '--quiet', action="store_true", help='do not produce any log output')
    parser.add_argument('-v', '--verbose', action="store_true", help='produce additional debugging output')
    parser.add_argument('-s', '--simple_coords', action="store_true", help='make coord numbers in debug output as small as possible')
    parser.add_argument('-a', '--include_anchors', action='store_true', help='include anchors in results table')
    args = parser.parse_args()

    # check if files exist
    with open(args.regions_file):
        pass
    if args.target_bedfile is not None:
        with open(args.target_bedfile):
            pass
    
    global log_level
    if args.quiet:
        log_level = LOG_LEVEL_QUIET
    elif args.verbose:
        log_level = LOG_LEVEL_DEBUG

    # define variables and create output directory
    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    #input("about to init ipp")
    log("Loading pairwise alignments")
    myIpp = ipp.Ipp()
    myIpp.load_pwalns(args.path_pwaln)
    myIpp.set_half_life_distance(args.half_life_distance)

    #input('Press enter to start')
    log('Reading regions from %s' %(args.regions_file))

    # Read the regions file and enqueue one projection job per line.
    ref_coords = []
    coord_names = {}
    with open(args.regions_file) as regions_file:
        for i,line in enumerate(regions_file.readlines()):
            cols = line.strip().split('\t')
            # define the coordinate as the center point between start and end coordinates
            name = cols[3]
            refChrom = cols[0]
            refLoc = int(np.mean([int(cols[1]), int(cols[2])]))
            coords = ipp.Coords((refChrom, refLoc)) 
            ref_coords.append(coords)
            # add the name of the region to a dict with refChrom:refLoc as the key for later translation
            coord_names[coords] = (i, name)

    global pbar
    pbar = tqdm.tqdm(total=len(ref_coords), leave=False)

    results = []
    unmapped_regions = []
    
    def on_job_done_callback(ref_coord,
                             direct_score,
                             direct_coords,
                             direct_up_anchor,
                             direct_down_anchor,
                             multi_shortest_path):
        # ref_coord:           ipp.Coords              (ref_chrom, ref_loc)
        # direct_score:        float                   [0-1]
        # direct_coords:       ipp.Coords              (qry_chrom, qry_loc)
        # direct_up_anchor:    ipp.Anchor              (up.ref_start, up.ref_end, up.qry_start, up.qry_end)
        # direct_down_anchor:  ipp.Anchor              (down.ref_start, down.ref_end, down.qry_start, down.qry_end)
        # multi_shortest_path: [ipp.ShortestPathEntry] [ref, intermediate1, intermediate2, ..., qry]
        #     species      string           "spX"
        #     score        float            [0-1]
        #     coords       ipp.Coords       (chrom, loc)
        #     up_anchor    ipp.Anchor       (up.ref_start, up.ref_end, up.qry_start, up.qry_end)
        #     down_anchor  ipp.Anchor       (down.ref_start, down.ref_end, down.qry_start, down.qry_end)
        nonlocal coord_names, results, unmapped_regions
        pbar.update()

        coord_idx, coord_name = coord_names[ref_coord]

        # handle unmapped region
        debug()
        if not len(multi_shortest_path):
            debug("no mapping found for {} ({}:{})".format(
                  coord_name, ref_coord.chrom, ref_coord.loc))
            debug()
            unmapped_regions.append(coord_name)
            return

        debug("({})".format(coord_name))
        debug_shortest_path(multi_shortest_path, args.simple_coords)
  
        direct_refs = tuple(map(str, (direct_up_anchor.ref_start, direct_up_anchor.ref_end,
                                      direct_down_anchor.ref_start, direct_down_anchor.ref_end))) \
            if direct_up_anchor else ("", "", "", "")
        direct_qrys = tuple(map(str, (direct_up_anchor.qry_start, direct_up_anchor.qry_end,
                                      direct_down_anchor.qry_start, direct_down_anchor.qry_end))) \
            if direct_up_anchor else ("", "", "", "")

        # The first intermediate species on the path.
        multi_first_entry = multi_shortest_path[1]

        # The qry species on the shortest path.
        multi_last_entry = multi_shortest_path[-1]
        assert multi_last_entry.species == args.qry
        multi_score = multi_last_entry.score
        multi_coords = multi_last_entry.coords
  
        multi_bridging_species = \
            ','.join([spe.species for spe in multi_shortest_path[1:-1]])
        results.append([
            coord_idx,
            coord_name, ref_coord, direct_coords, multi_coords,
            direct_score, multi_score,
            multi_bridging_species,
            *direct_refs,
            multi_first_entry.up_anchor.ref_start, multi_first_entry.up_anchor.ref_end,
            multi_first_entry.down_anchor.ref_start, multi_first_entry.down_anchor.ref_end,
            *direct_qrys,
            multi_last_entry.up_anchor.qry_start, multi_last_entry.up_anchor.qry_end,
            multi_last_entry.down_anchor.qry_start, multi_last_entry.down_anchor.qry_end])
  
    # Start the projection
    log('Projecting regions from %s to %s' %(args.ref, args.qry))
    myIpp.project_coords(args.ref,
                         args.qry,
                         ref_coords,
                         args.n_cores,
                         on_job_done_callback)
    pbar.close()

    # create data frame from results dict
    anchor_cols = ['ref_anchor_direct_left_start', 'ref_anchor_direct_left_end', 'ref_anchor_direct_right_start', 'ref_anchor_direct_right_end',
                   'ref_anchor_multi_left_start', 'ref_anchor_multi_left_end', 'ref_anchor_multi_right_start', 'ref_anchor_multi_right_end',
                   'qry_anchor_direct_left_start', 'qry_anchor_direct_left_end', 'qry_anchor_direct_right_start', 'qry_anchor_direct_right_end',
                   'qry_anchor_multi_left_start', 'qry_anchor_multi_left_end', 'qry_anchor_multi_right_start', 'qry_anchor_multi_right_end']
    results_df = pd.DataFrame(results, columns=[
        'sort_idx',
        'id', 'coords_ref', 'coords_direct', 'coords_multi',
        'score_direct', 'score_multi', 'bridging_species',
        *anchor_cols]) \
            .sort_values(by=['sort_idx']) \
            .drop(columns=['sort_idx']) \
            .set_index('id')

    # exclude anchor columns if flag to keep them was not set
    if not args.include_anchors:
        columns=['coords_ref', 'coords_direct', 'coords_multi',
                 'score_direct', 'score_multi', 'bridging_species']
        results_df = results_df.loc[:, columns]

    # Trim scores to third decimal (floor) such that 0.9999 becomes 0.999 instead of 1.000
    # A score of 1.0 is reserved for regions overlapping alignments throughout the path
    results_df[['score_direct','score_multi']] = results_df[['score_direct','score_multi']].apply(lambda x: np.floor(x*1000)/1000)
        
    if is_debug():
        debug(results_df.to_string())

    # write list of coord names of unmapped regions to file
    regions_file_basename = os.path.splitext(os.path.basename(args.regions_file))[0]
    outfile_unmapped = os.path.join(args.out_dir, '{}.{}-{}.unmapped'.format(regions_file_basename, args.ref, args.qry))
    log('Writing unmapped regions to:\n\t%s' %outfile_unmapped)
    with open(outfile_unmapped, 'w') as f:
        f.write('\n'.join(unmapped_regions) + '\n')

    def classify_conservation(df_projections, target_regions=pr.PyRanges(), thresh=.95, maxgap=500):
        ### function for determining the conservation of sequence (DC/IC/NC) and function (+/-)  
        # determine sequence conservation
        sequence_conservation = df_projections.apply(lambda x: 'DC' if x['score_direct'] == 1 else 'IC' if x['score_multi'] >= thresh else 'NC', axis=1)
        # create PyRanges object of projections
        df_multi = pd.DataFrame({'Name':df_projections.index,
                                 'Chromosome': [x.chrom for x in df_projections['coords_multi']],
                                 'Start': [x.loc for x in df_projections['coords_multi']]})
        df_multi['End'] = df_multi['Start'] + 1
        pr_multi = pr.PyRanges(df_multi)  
        # determine functional conservation by checking for overlap with target regions
        functional_conservation = [''] * df_multi.shape[0]
        if target_regions is not None:
            # expand target regions by maxgap
            target_regions.Start -= maxgap
            target_regions.End += maxgap
            # save names of projected regions that overlap with a target region
            names_FC = pr_multi.overlap(target_regions).Name.values
            functional_conservation = pd.Series('-', index=df_projections.index)
            functional_conservation[names_FC] = '+'
        return sequence_conservation, functional_conservation
        
    # classify projections according to conservation of sequence (DC/IC/NC) and function (+/-)
    log('Classifying projections')
    thresh = .95
    maxgap = 500
    target_regions = None
    if args.target_bedfile is not None:
        target_regions = pr.read_bed(args.target_bedfile)
    sequence_conservation, functional_conservation = classify_conservation(results_df, target_regions, thresh, maxgap)
    results_df.insert(5, 'sequence_conservation', sequence_conservation)
    results_df.insert(6, 'functional_conservation', functional_conservation)

    # write results table to file
    outfile_table = os.path.join(args.out_dir, '{}.{}-{}.proj'.format(regions_file_basename, args.ref, args.qry))
    if results_df.shape[0] == 0:
        log('No regions from the input bed files could be projected\nDone')
        return
    log('Writing table with successful projections to:\n\t%s' %outfile_table)
    results_df.to_csv(outfile_table, sep='\t', header=True, float_format='%.3f')
    
    def convert_df_to_bed(df, which):
        # function to convert the projections from data frame to bed format
        # which must be 'ref' or 'qry'
        colors = {'DC':'127,201,127', 'IC':'253,180,98', 'NC':'141,153,174', 'DC+':'255,0,0', 'IC+':'255,0,0', 'DC-':'30,30,30', 'IC-':'30,30,30', 'NC+':'141,153,174', 'NC-':'141,153,174'}
        coords_col = 'coords_ref'
        opposite_coords_col = 'coords_multi'
        if which == 'qry':
            coords_col, opposite_coords_col = opposite_coords_col, coords_col
        bed = pd.DataFrame(df[coords_col].tolist(), index=df.index).rename(columns={0:'chrom', 1:'start'})
        bed['start'] = bed['start'].astype(int)
        bed['end'] = bed['start'] + 1
        bed['name'] = df.apply(lambda x: '_'.join([x.name, str(x[opposite_coords_col]), x['sequence_conservation']+x['functional_conservation']]), axis=1)
        bed['score'] = df['score_multi']
        bed['strand'] = '.'
        bed['thickStart'] = bed['start']
        bed['thickEnd'] = bed['end']
        bed['itemRgb'] = (df['sequence_conservation'] + df['functional_conservation']).apply(lambda x: colors[x])
        return bed
    
    # write results to bed files (bed9 format for both reference and target species)
    # name column (4):    string    "id_qryChrom:qryLoc_class"    e.g. "peak_75_chr3:3880_IC+"
    outfile_bed_ref = '{}.{}.bed'.format(os.path.join(args.out_dir, regions_file_basename), args.ref)
    outfile_bed_qry = '{}.{}.bed'.format(os.path.join(args.out_dir, regions_file_basename), args.qry)
    log('Writing output bed files:\n\t%s\n\t%s' %(outfile_bed_ref, outfile_bed_qry))
    bed_ref = convert_df_to_bed(results_df, 'ref')
    bed_qry = convert_df_to_bed(results_df, 'qry')
    bed_ref.to_csv(outfile_bed_ref, sep='\t', index=False, header=None, float_format='%.3f')
    bed_qry.to_csv(outfile_bed_qry, sep='\t', index=False, header=None, float_format='%.3f')

    log('Done')
    return
    
if __name__ == '__main__':
    main()

# vim: tabstop=4 shiftwidth=4 expandtab
