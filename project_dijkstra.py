#! /usr/bin/env python

import numpy as np, pandas as pd, sys, pickle
# import multiprocessing
# from joblib import Parallel, delayed
from functions_genomic_coordinate_projection import *

def project(coord, ref, qry, species, pwaln, genome_size, scaling_factor):
     # direct projection
     score_direct_projection, coords_direct_projection, ref_anchors_direct, qry_anchors_direct = project_genomic_location(ref, qry, coord, 1.0, pwaln, genome_size, scaling_factor)
     ref_anchors_direct_left, ref_anchors_direct_right = ref_anchors_direct.split(',')
     qry_anchors_direct_left, qry_anchors_direct_right = qry_anchors_direct.split(',')
          
     # multi-species projection
     shortest_path_to_qry, shortest_path, orange = get_shortest_path(ref, qry, coord, species, pwaln, genome_size, scaling_factor, verbose=False)
     score_multi_projection, coords_multi_projection = shortest_path_to_qry.loc[qry,['score','coords']]
     ref_anchors_multi_left, ref_anchors_multi_right = shortest_path_to_qry['ref_anchors'][1].split(',') # ref anchors of the first species in the path (first non-reference species)
     qry_anchors_multi_left, qry_anchors_multi_right = shortest_path_to_qry['qry_anchors'][-1].split(',') # qry anchors of the last species in the path
     bridging_species = ','.join(shortest_path_to_qry.index.values[1:-1])
     values = [coord, coords_direct_projection, coords_multi_projection,
               score_direct_projection, score_multi_projection,
               ref_anchors_direct_left, ref_anchors_direct_right, ref_anchors_multi_left, ref_anchors_multi_right,
               qry_anchors_direct_left, qry_anchors_direct_right, qry_anchors_multi_left, qry_anchors_multi_right,
               bridging_species]
     columns = ['coords_ref', 'coords_direct', 'coords_multi', 'score_direct', 'score_multi', 'ref_anchor_direct_left', 'ref_anchor_direct_right', 'ref_anchor_multi_left', 'ref_anchor_multi_right', 'qry_anchor_direct_left', 'qry_anchor_direct_right', 'qry_anchor_multi_left', 'qry_anchor_multi_right', 'bridging_species']
     df = pd.DataFrame(values, index=columns).T
     return df

def main():
     if len(sys.argv) < 8:
          print('Usage: python project_dijkstra.py reference_species query_species coord<chrN:xxxxxx> id half_life_distance path_pwaln_pkl tmp_dir [optional: -quiet]')
          sys.exit(0)
     _, ref, qry, coord, coord_id, half_life_distance, path_pwaln_pkl, tmp_dir = sys.argv[:8]
     quiet = False
     if len(sys.argv) == 9 and sys.argv[8] == '-quiet':
          quiet = True

     # define paths
     data_dir = '/project/MDL_ChIPseq/data/genome'
     assembly_dir = data_dir + '/assembly/'
     
     # read pwaln and genome size files
     with open(path_pwaln_pkl, 'rb') as pkl: # you can also pass the ce.pkl path to path_pwaln_pkl, so this works for both pairwise alignments and C(N)Es.
          pwaln = pickle.load(pkl)
     species = np.array(list(pwaln.keys())) # The list of species is determined based on the keys of the supplied pairwise aln (pwaln) dict
     genome_size = {s : read_genome_size(assembly_dir + s + '.sizes') for s in species}

     # determine scaling factor based on desired distance_half_life (at which distance to an anchor in the reference species is the score supposed to be 0.5)
     scaling_factor = get_scaling_factor(genome_size[ref], int(half_life_distance))
     
     # project coordinates
     try:
          df = project(coord, ref, qry, species, pwaln, genome_size, scaling_factor)
     except KeyError: # if not even dijkstra finds a projection (mostly due to chromosome border regions)
          if not quiet:
               print('Unable to map %s: %s' %(coord_id, coord))
          return
     df.index = [coord_id]
     if not os.path.exists(tmp_dir):
          os.mkdir(tmp_dir)
     df.to_csv('%s/%s.proj.tmp' %(tmp_dir,coord_id), sep='\t', header=False)
     
     return

if __name__ == '__main__':
     main()
