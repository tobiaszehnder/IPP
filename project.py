#! /usr/bin/env python

import numpy as np, pandas as pd, sys, pickle, argparse, multiprocessing as mp, ctypes, traceback, tqdm, liftover, pyranges as pr
# from joblib import Parallel, delayed
from functions import *

# set global variables
pwaln = {} # assigned in load_pwaln() for main process and init_worker() for worker process
chroms = [] # assigned in load_pwaln()

def project_coord(coord, ref, qry, genome_size, scaling_factor):
  # function for getting direct and multi projections for a given genomic coordinate
  global chroms
  coords_ref = '{}:{}'.format(coord.chrom_idx, coord.loc)
  
  # direct projection
  score_direct, coords_direct = project_genomic_location(
    ref, qry, coords_ref, 1.0, pwaln, genome_size, scaling_factor
  )

  # multi-species projection
  species = pwaln.keys()
#   shortest_path_to_qry, shortest_path, orange = ...
  shortest_path_to_qry = get_shortest_path(
    ref, qry, coords_ref, species, pwaln, genome_size, scaling_factor, verbose=False
  )
  if qry not in shortest_path_to_qry.index:
    return pd.DataFrame(index=[coord.id])
  score_multi, coords_multi = shortest_path_to_qry.loc[qry,['score','coords']]
  bridging_species = ','.join(shortest_path_to_qry.index.values[1:-1]) # everything except first and last
  
  # combine direct and multi projections and return as a dataframe
  columns = ['coords_ref', 'coords_direct', 'coords_multi', 'score_direct', 'score_multi', 'bridging_species']
  values = [coords_ref, coords_direct, coords_multi, score_direct, score_multi, bridging_species]
  df = pd.DataFrame([values], index=[coord.id], columns=columns)
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
  parser.add_argument('regions_file', help='bed file containing genomic coordinates. regions with width > 1 will be centered.')
  parser.add_argument('ref', help='reference species')
  parser.add_argument('qry', help='query species')
  parser.add_argument('pwaln_pickle_file', help='pairwise alignments collection')
  parser.add_argument('chain_file', help='reference to query chain file for liftover')
  parser.add_argument('--target_bedfile', default=argparse.SUPPRESS, help='path to a bedfile with equivalent regions in the target species used for the classification of functional conservation.')
  parser.add_argument('--out_dir', default=os.getcwd(), help='directory for output files')
  parser.add_argument('--half_life_distance', type=int, default=10000, help='distance to closest anchor point at which projection score is 0.5')
  parser.add_argument('--n_cores', type=int, default=1, help='number of CPUs')
  parser.add_argument('--quiet', help='do not produce any verbose output')
  parser.add_argument('--data_dir', default='/project/ipp-data')
  args = parser.parse_args()

  # define variables and create output directory
  assembly_dir = args.data_dir + '/assembly/'
  if not os.path.exists(args.out_dir):
    os.mkdir(args.out_dir)
  regions_file_basename = os.path.splitext(os.path.basename(args.regions_file))[0]
  columns = ['coords_ref', 'coords_direct', 'coords_multi', 'score_direct', 'score_multi', 'bridging_species']

  # load pairwise alignment collection pickle and copy to shared memory
  print('Loading pairwise alignments from pickle')
  load_pwaln(args.pwaln_pickle_file)  

  # load chromosome sizes
  genome_size = {s: read_genome_size(assembly_dir + s + '.sizes') for s in pwaln.keys()}
  
  # determine scaling factor based on desired distance_half_life
  # (at which distance to an anchor in the reference species is the score supposed to be 0.5)
  scaling_factor = get_scaling_factor(
              genome_size[args.ref], 
              int(args.half_life_distance))
  
  # read regions file
  with open(args.regions_file, 'r') as regions_file:
    lines = [line.strip().split('\t') for line in regions_file.readlines()]
    coords = np.array([Coord(line[0], chroms.index(line[0]), np.mean([int(line[1]), int(line[2])]), line[3]) for line in lines])
    
  # liftover alignable regions
  print('Liftover alignable regions')
  converter = liftover.ChainFile(args.chain_file, args.ref, args.qry)
  liftover_coords = [converter[x.chrom][x.loc] for x in coords]
  idx_successful_liftover = [i for i,x in enumerate(liftover_coords) if x]
  idx_ipp = [i for i,x in enumerate(liftover_coords) if not x] # indices for IPP projections
  results_liftover = pd.DataFrame([['{}:{}'.format(coords[i].chrom, coords[i].loc),
                                    '{}:{}'.format(*liftover_coords[i][0][:2]),
                                    '{}:{}'.format(*liftover_coords[i][0][:2]),
                                    1.,
                                    1.,
                                    '']
                                   for i in idx_successful_liftover], columns=columns)
  results_liftover.index = [x.id for x in coords[idx_successful_liftover]]

#   input('Press enter to start')
  print('Projecting non-alignable regions from %s to %s' %(args.ref, args.qry))
  # create a process pool
  pool = mp.Pool(processes=args.n_cores,
                 initializer=init_worker,
                 initargs=[pwaln_shared])
  
  # callback function to append job result to collection of results
  global results, unmapped
  results = pd.DataFrame()
  unmapped = []
  def handle_job_complete(result):
    global results, unmapped, pbar
    # if result is empty add it's ID to the list of unmapped regions
    if result.shape[1] == 0:
      unmapped += result.index.tolist()
    else:
      results = results.append(result)
    pbar.update()
    
  # error callback function in case project_coord() fails
  def handle_job_error(e):
    global pbar
    pbar.update()
    try:
      raise # e
    except KeyError:
      if not args.quiet:
        print(e)
        traceback.print_exc()
        
  # Enqueue one projection job per coordinate that wasn't projected by liftover
  global pbar
  pbar = tqdm.tqdm(total=len(idx_ipp))
  for i in idx_ipp:
    coord = coords[i]
#     coord = '{}:{}'.format(coords[i].chrom_idx, coords[i].loc)
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
  
  # reformat chromosome names: '0:123456' to 'chr1:123456'
  coord_cols = ['coords_ref', 'coords_direct', 'coords_multi']
  results.loc[:,coord_cols] = results.loc[:,coord_cols].apply(lambda x: [reformat_coordinate(y, chroms) if ':' in y else y for y in x])
  ids_ipped_regions = np.setdiff1d([x.id for x in coords[idx_ipp]], unmapped)
  results.index = ids_ipped_regions
  
  # add liftover projections to ipp projections table
  ids_all_mapped_regions = np.setdiff1d([x.id for x in coords], unmapped)
  results = pd.concat([results, results_liftover]).loc[ids_all_mapped_regions,:]
  
  # write projections table to file
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
  
  # write projections to file (bed9 format for both reference and target species)
  # the name columns contains id, coordinates of the opposite species and the class
  outfile_bed_ref = '{}.{}.bed'.format(regions_file_basename, args.ref)
  outfile_bed_qry = '{}.{}.bed'.format(regions_file_basename, args.qry)
  bed_ref = results.apply(lambda row: format_row_table_to_bed(row, 'ref'), axis=1)
  bed_qry = results.apply(lambda row: format_row_table_to_bed(row, 'qry'), axis=1)
  bed_ref.to_csv(outfile_bed_ref, sep='\t', index=False, header=None, float_format='%.3f')
  bed_qry.to_csv(outfile_bed_qry, sep='\t', index=False, header=None, float_format='%.3f')
  
if __name__ == '__main__':
  main()

# vim: tabstop=2 shiftwidth=2 expandtab
