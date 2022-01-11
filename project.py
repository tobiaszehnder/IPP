#!/usr/bin/env python
import argparse
import ipp
import numpy as np
import os
import pandas as pd
import tqdm

def main():
    parser = argparse.ArgumentParser(description='Independent Point Projections (IPP).\nA method for projecting genomic point coordinates between genomes with large evolutionary distances.')
    parser.add_argument('regions_file', help='Bed file containing genomic coordinates. regions with width > 1 will be centered.')
    parser.add_argument('ref', help='reference species')
    parser.add_argument('qry', help='query species')
    parser.add_argument('path_pwaln')
    parser.add_argument('--out_dir', default=os.getcwd(), help='directory for output files')
    parser.add_argument('--half_life_distance', type=int, default=10000, help='distance to closest anchor point at which projection score is 0.5')
    parser.add_argument('--n_cores', type=int, default=1, help='number of CPUs')
    parser.add_argument('--quiet', help='do not produce any verbose output')
    parser.add_argument('--data_dir', default='/project/ipp-data')
    args = parser.parse_args()

    # define paths
    assembly_dir = args.data_dir + '/assembly/'

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    #input("about to init ipp")
    myIpp = ipp.Ipp()
    myIpp.load_pwalns(args.path_pwaln)
    myIpp.load_genome_sizes(assembly_dir);
    myIpp.set_half_life_distance(args.half_life_distance)

    #input('Press enter to start')
    print('Projecting regions from %s to %s' %(args.ref, args.qry))

    # Read the regions file and enqueue one projection job per line.
    ref_coords = []
    with open(args.regions_file) as regions_file:
        for line in regions_file.readlines():
            cols = line.split('\t')
            # define the coordinate as the center point between start and end coordinates
            refChrom = cols[0]
            refLoc = int(np.mean([int(cols[1]), int(cols[2])]))
            ref_coords.append((refChrom, refLoc))

    pbar = tqdm.tqdm(total=len(ref_coords), leave=False)
    results = pd.DataFrame(
        columns=['coords_ref', 'coords_direct', 'coords_multi',
                 'score_direct', 'score_multi',
                 'ref_anchor_direct_left', 'ref_anchor_direct_right', 'ref_anchor_multi_left', 'ref_anchor_multi_right',
                 'qry_anchor_direct_left', 'qry_anchor_direct_right', 'qry_anchor_multi_left', 'qry_anchor_multi_right',
                 'bridging_species'])
    def on_job_done_callback(ref_coord,
                             direct_score,
                             direct_coords,
                             direct_ref_anchors,
                             direct_qry_anchors,
                             multi_shortest_path):
        assert multi_shortest_path[-1][0] == args.qry
        multi_score = multi_shortest_path[-1][1]
        multi_coords = multi_shortest_path[-1][2]
  
        # Ref anchors of the first species in the path (first non-reference species)
        multi_ref_anchors = multi_shortest_path[1][3]
  
        # Qry anchors of the last species in the path.
        multi_qry_anchors = multi_shortest_path[-1][4]
  
        multi_bridging_species = \
            ','.join([spe[0] for spe in multi_shortest_path[1:-1]])
  
        values = [ref_coord, direct_coords, multi_coords,
                  direct_score, multi_score,
                  direct_ref_anchors[0], direct_ref_anchors[1], multi_ref_anchors[0], multi_ref_anchors[1],
                  direct_qry_anchors[0], direct_qry_anchors[1], multi_qry_anchors[0], multi_qry_anchors[1],
                  multi_bridging_species]
        nonlocal results
        results = results.append(
                pd.DataFrame([values], columns=results.columns),
                ignore_index=True)
  
        pbar.update()
  
    # Start the projection
    myIpp.project_coords(args.ref,
                         args.qry,
                         ref_coords,
                         args.n_cores,
                         on_job_done_callback)
    pbar.close()
  
    print(results)
  
    ### TO DO:
    ### name columns here and not in project_coord function
    ### use ids from bed column 4 as index names (this should be part of the Coord class anyways, implement with that)
    
    outfile = os.path.splitext(os.path.basename(args.regions_file))[0] + '.proj'
    results.to_csv(os.path.join(args.out_dir, outfile), sep='\t', header=True)
  
if __name__ == '__main__':
    main()

# vim: tabstop=4 shiftwidth=4 expandtab
