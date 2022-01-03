#! /usr/bin/env python

import numpy as np, pandas as pd, sys, pickle, argparse, multiprocessing as mp, ctypes, traceback, tqdm
# from joblib import Parallel, delayed
from functions import *

# set global variables
pwaln = {} # assigned in load_pwaln() for main process and init_worker() for worker process
chroms = [] # assigned in load_pwaln()

def project_coord_new(coord, ref, qry, genome_size, scaling_factor):
  return

def project_coord(coord, ref, qry, genome_size, scaling_factor):
  # direct projection
  score_direct_projection, coords_direct_projection, ref_anchors_direct, qry_anchors_direct = project_genomic_location(ref, qry, coord, 1.0, pwaln, genome_size, scaling_factor)
  ref_anchors_direct_left, ref_anchors_direct_right = ref_anchors_direct.split(',')
  qry_anchors_direct_left, qry_anchors_direct_right = qry_anchors_direct.split(',')

  # multi-species projection
  species = pwaln.keys()
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

def load_pwaln(path_pwaln_pkl):
  # load pwaln pickle and copy to shared memory
  global pwaln
  global pwaln_shared
  global chroms
  global columns
  
  with open(path_pwaln_pkl, "rb") as f:
    data = pickle.load(f)
    pwaln = data["pwaln"]
    chroms = data["chromosomes"]
    
  pwaln_shared = {}
  for sp1, v in pwaln.items():
      pwaln_shared[sp1] = {}
      for sp2, df in v.items():
          columns = df.columns
          # Remember the columns. They are required to be the same in all data
          # frames.

          # reserve the shared memory
          mp_arr = mp.RawArray(ctypes.c_uint32, int(df.size))

          # copy the data from df into the shared memory.
          mp_arr_np = np.frombuffer(mp_arr, dtype=np.uint32).reshape(df.shape)
          np.copyto(mp_arr_np, df.to_numpy(dtype=np.uint32))
          pwaln_shared[sp1][sp2] = mp_arr

def init_worker(pwaln_shared):
    # Create a wrapper pwaln object that accesses its data from the shared
    # memory.
    global pwaln
    for sp1, v in pwaln_shared.items():
        pwaln[sp1] = {}
        for sp2, mp_arr in v.items():
            mp_arr_np = np.frombuffer(mp_arr, dtype=np.uint32) \
                .reshape([int(len(mp_arr)/len(columns)), len(columns)])
            pwaln[sp1][sp2] = pd.DataFrame(mp_arr_np, columns=columns)
    
def main():
  parser = argparse.ArgumentParser(description='Independent Point Projections (IPP).\nA method for projecting genomic point coordinates between genomes with large evolutionary distances.')
  parser.add_argument('regions_file', help='Bed file containing genomic coordinates. regions with width > 1 will be centered.')
  parser.add_argument('ref', help='reference species')
  parser.add_argument('qry', help='query species')
  parser.add_argument('path_pwaln_pkl')
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

  # load pwaln pickle and copy to shared memory
  print('Loading pairwise alignments from pickle')
  load_pwaln(args.path_pwaln_pkl)
  
  # The list of species is determined based on the keys of the supplied pwaln dict.
  genome_size = {s: read_genome_size(assembly_dir + s + '.sizes') for s in pwaln.keys()}
  
  # determine scaling factor based on desired distance_half_life (at which distance to an anchor in the reference species is the score supposed to be 0.5)
  scaling_factor = get_scaling_factor(
              genome_size[args.ref], 
              int(args.half_life_distance))

#   input('Press enter to start')
  print('Projecting regions from %s to %s' %(args.ref, args.qry))
  # create a process pool
  pool = mp.Pool(processes=args.n_cores,
                 initializer=init_worker,
                 initargs=[pwaln_shared])
  
  # callback function to append job result to collection of results
  global results
  results = pd.DataFrame()
  def handle_job_complete(result):
    global results, pbar
    results = results.append(result)
    pbar.update()
    
  # error callback function in case project() fails
  def handle_job_error(e):
    global pbar
    pbar.update()
    try:
      raise e
    except KeyError:
      # if not even dijkstra finds a projection (mostly due to chromosome border regions)
      if not args.quiet:
        print('Unable to map some region') #%s: %s' %(coord_id, coord))
        print(e)
        traceback.print_exc()
        
  # Read the regions file and enqueue one projection job per line.
  nregions = sum(1 for _ in open(args.regions_file))
  global pbar
  pbar = tqdm.tqdm(total=nregions)
  with open(args.regions_file) as regions_file:
    for line in regions_file.readlines():
      cols = line.split('\t')
#       coord = Coord(chroms.index(cols[0]), int(cols[1]), cols[3])
      # define the coordinate as the center point between start and end coordinates
      coord = '{}:{}'.format(chroms.index(cols[0]), int(np.mean([int(cols[1]), int(cols[2])])))
      job_args = (coord,
                  args.ref,
                  args.qry,
                  genome_size,
                  scaling_factor)
      pool.apply_async(project_coord,
                       job_args,
                       callback=handle_job_complete,
                       error_callback=handle_job_error)
      
  # Wait for the jobs to complete.
  pool.close()
  pool.join()
  
  ### TO DO:
  ### name columns here and not in project_coord function
  ### use ids from bed column 4 as index names (this should be part of the Coord class anyways, implement with that)
  
  # reformat '0:123456' to 'chr1:123456'
#   results.columns = np.arange(len(results.columns))
  coord_cols = ['coords_ref', 'coords_direct', 'coords_multi']
  results.loc[:,coord_cols] = results.loc[:,coord_cols].apply(lambda x: ['{}:{}'.format(chroms[int(y.split(':')[0])], y.split(':')[1]) for y in x])
  outfile = os.path.splitext(os.path.basename(args.regions_file))[0] + '.proj'
  results.to_csv(os.path.join(args.out_dir, outfile), sep='\t', header=True)
  
if __name__ == '__main__':
  main()

# vim: tabstop=2 shiftwidth=2 expandtab
