import numpy as np, pandas as pd, sys, os, glob, heapq, heapq_max, re, itertools, copy

class Coord:
  def __init__(self, chrom, loc, id_):
    self.chrom = chrom
    self.loc = loc
    self.id = id_

def flatten_list(l):
  return [item for sublist in l for item in sublist]

def read_genome_size(filename):
  return float(pd.read_csv(filename, sep='\t', usecols=[1]).values.sum())
  
def longest_sorted_subsequence(seq):
  # wrapper for longest_increasingly_sorted_subsequence(), applying it for increasing and decreasing and returning the longer resulting subsequence.
  # 0.1 - 1 millisec for any sequence of length 20
  seq = list(seq)
  seq_neg = list(-np.array(seq))
  longest_sub = [np.array(longest_increasingly_sorted_subsequence(seq)), -np.array(longest_increasingly_sorted_subsequence(seq_neg))] # increasing and decreasing
  return longest_sub[np.argmax([len(x) for x in longest_sub])]

def longest_increasingly_sorted_subsequence(seq):
    # really fast stackoverflow implementation (50-100 microsec for any sequence of length 20)
    seq = list(seq)
    if not seq:
        return seq

    M = [None] * len(seq)    # offset by 1 (j -> j-1)
    P = [None] * len(seq)

    # Since we have at least one element in our list, we can start by 
    # knowing that the there's at least an increasing subsequence of length one:
    # the first element.
    L = 1
    M[0] = 0

    # Looping over the sequence starting from the second element
    for i in range(1, len(seq)):
        # Binary search: we want the largest j <= L
        #  such that seq[M[j]] < seq[i] (default j = 0),
        #  hence we want the lower bound at the end of the search process.
        lower = 0
        upper = L

        # Since the binary search will not look at the upper bound value,
        # we'll have to check that manually
        if seq[M[upper-1]] < seq[i]:
            j = upper

        else:
            # actual binary search loop
            while upper - lower > 1:
                mid = (upper + lower) // 2
                if seq[M[mid-1]] < seq[i]:
                    lower = mid
                else:
                    upper = mid

            j = lower    # this will also set the default value to 0

        P[i] = M[j-1]

        if j == L or seq[i] < seq[M[j]]:
            M[j] = i
            L = max(L, j+1)

    # Building the result: [seq[M[L-1]], seq[P[M[L-1]]], seq[P[P[M[L-1]]]], ...]
    result = []
    pos = M[L-1]
    for _ in range(L):
        result.append(seq[pos])
        pos = P[pos]

    return result[::-1]    # reversing

def get_anchors(df, chrom, x, return_top_n=False): 
  # return_top_n: do not just return the two closest anchors, but the list of topn anchors  
  anchor_cols = ['ref_chrom','ref_start','ref_end','ref_coord','qry_chrom','qry_start','qry_end','qry_coord','qry_strand']
  ov_aln = df.loc[(df.ref_chrom == chrom) & (df.ref_start < x) & (df.ref_end > x),].reset_index(drop=True) # x lies on an alignment block. return x itself as both anchors
  ov_aln['ref_coord'] = (ov_aln.ref_start + ov_aln.ref_end) / 2 # add the center of the alignment as the `ref_coord` column in order to check for collinearity with up- and downstream anchors later
  ov_aln['qry_coord'] = (ov_aln.qry_start + ov_aln.qry_end) / 2 # add the center of the alignment as the `qry_coord` column in order to check for collinearity with up- and downstream anchors later
  
  # first define anchors upstream, downstream and ov_aln, then do major_chrom / collinearity test, then either return overlapping anchor or closest anchors.
  # take orientation into account for the anchor definition. if start > end, then the aln is to the '-' strand.
  # for speed reasons only select the first 100. the rest just takes longer to compute min / max and most likely will (and should) not be an anchor anyways.
  anchors_upstream = df.loc[abs(df.loc[(df.ref_chrom == chrom) & (df.ref_end < x),].ref_end - x).sort_values().index,['ref_chrom','ref_start','ref_end','qry_chrom','qry_start','qry_end']].iloc[:100,:].reset_index(drop=True) # keeping index slows column creation later
  anchors_downstream = df.loc[abs(df.loc[(df.ref_chrom == chrom) & (df.ref_start > x),].ref_start - x).sort_values().index,['ref_chrom','ref_start','ref_end','qry_chrom','qry_start','qry_end']].iloc[:100,:].reset_index(drop=True)
  anchors_upstream.insert(3, 'ref_coord', anchors_upstream.ref_end)
  anchors_downstream.insert(3, 'ref_coord', anchors_downstream.ref_start)
  # set the corresponding start / end coordinate that is closer to the projected coordinate as the qry_coord. (choosing the SAME LABELED COORDINATE as the reference, i.e. START in case of downsteam and END in case of upstream)
  anchors_upstream['qry_coord'] = anchors_upstream.qry_end
  anchors_downstream['qry_coord'] = anchors_downstream.qry_start
  if min(anchors_upstream.shape[0], anchors_downstream.shape[0]) < 1: # require minimum of 1 anchor on each side. later, the minimum total number of collinear anchors will be set to 5 (but one side is allowed to have as little as 1 anchor).
    return pd.DataFrame(columns=anchor_cols)

  # test collinearity of anchors: take top 20 in each direction (top 10 produced many locally collinear pwalns that were still non-collinear outliers in the global view of the GRB.)
  # note: using ungapped chain blocks might require n to be even larger.
  minn = 5
  topn = 20
  
  # MAJOR CHROMOSOME: retain anchors that point to the majority chromosome in top n of both up- and downstream anchors
  major_chrom = pd.concat([anchors_upstream[:topn], ov_aln, anchors_downstream[:topn]], axis=0, sort=False).qry_chrom.value_counts().idxmax()
  ov_aln[ov_aln.qry_chrom == major_chrom]
  anchors_upstream = anchors_upstream[anchors_upstream.qry_chrom == major_chrom]
  anchors_downstream = anchors_downstream[anchors_downstream.qry_chrom == major_chrom]
  
  # COLLINEARITY: remove pwalns pointing to outliers by getting the longest sorted subsequence of the top n of both up- and downstream anchors.
  closest_anchors = pd.concat([anchors_upstream[:topn][::-1], ov_aln, anchors_downstream[:topn]], axis=0, sort=False).reset_index(drop=True) # reset_index necessary, otherwise working with duplicate indices messing things up
  idx_collinear = closest_anchors.index[np.intersect1d(closest_anchors.qry_coord.values.astype(int), longest_sorted_subsequence(closest_anchors.qry_coord.values.astype(int)), return_indices=True)[1]] # this step takes 2 sec
  closest_anchors = closest_anchors.loc[idx_collinear,].dropna(axis=1, how='all') # drop columns if it only contains NaNs (see explanation below)
  # set minimum number of collinear anchors to 5 (for species pairs with very large evol. distances setting a lower boundary for the number of collinear anchors will help reduce false positives)
  if closest_anchors.shape[0] < minn:
    return pd.DataFrame(columns=anchor_cols)
  
  # return the top anchors if the flag was passed
  if return_top_n:
    return closest_anchors[['ref_chrom','ref_start','ref_end','qry_chrom','qry_start','qry_end']]
  
  # check if the original ov_aln is still present (or ever was) in the filtered closest_anchors (that include a potential ov_aln)
  # if not, it was an outlier alignment and was filtered out
  # if yes, narrow it to the actual position of x and its relative position in the qry such that the returned anchors have distance = 0 to x
  ov_aln = closest_anchors.loc[(closest_anchors.ref_start < x) & (closest_anchors.ref_end > x),].reset_index(drop=True)
  if ov_aln.shape[0] == 1:
    x_relative = x - ov_aln.ref_start[0]
    strand = '+' if ov_aln.qry_start[0] < ov_aln.qry_end[0] else '-'
    
    vals_up = ov_aln.copy().rename(index={0:'upstream'})
    vals_up.ref_start = x-1
    vals_up.ref_end = x
    vals_up.ref_coord = x
    vals_up.qry_start = ov_aln.qry_start[0] + x_relative - 1 if strand == '+' else ov_aln.qry_start[0] - x_relative + 2
    vals_up.qry_end = vals_up.qry_start + 1 if strand == '+' else vals_up.qry_start - 1
    vals_up.ref_coord = vals_up.ref_end
    vals_up.qry_coord = vals_up.qry_end

    vals_down = vals_up.copy().rename(index={'upstream':'downstream'})
    vals_down.loc[:,('ref_start', 'ref_end')] += 1
    qry_shift = 1 if strand == '+' else -1
    vals_down.loc[:,('qry_start', 'qry_end')] += qry_shift
    vals_down.ref_coord = vals_down.ref_start
    vals_down.qry_coord = vals_down.qry_start

    anchors = pd.concat([pd.concat([vals_up, vals_down], axis=0), pd.DataFrame({'qry_strand':strand}, index=['upstream','downstream'])], axis=1)
  else:
    anchor_upstream = closest_anchors.loc[abs(closest_anchors.loc[closest_anchors.ref_coord < x,].ref_coord - x).sort_values().index,].head(1).rename(index=lambda x:'upstream')
    if anchor_upstream.shape[0] > 0:
      anchor_upstream['qry_strand'] = '+' if anchor_upstream.qry_start[0] < anchor_upstream.qry_end[0] else '-'
    anchor_downstream = closest_anchors.loc[abs(closest_anchors.loc[closest_anchors.ref_coord > x,].ref_coord - x).sort_values().index,].head(1).rename(index=lambda x:'downstream')
    if anchor_downstream.shape[0] > 0:
      anchor_downstream['qry_strand'] = '+' if anchor_downstream.qry_start[0] < anchor_downstream.qry_end[0] else '-'
    if not anchor_upstream.shape[0] + anchor_downstream.shape[0] == 2: # if not both up- and downstream anchors were found (e.g. at synteny break points where one side does not have any anchors to the majority chromosome)
      return pd.DataFrame(columns=anchor_cols)
    anchors = pd.concat([anchor_upstream, anchor_downstream]).loc[:,anchor_cols]
  anchors.loc[:,['ref_coord','qry_coord']] = anchors.loc[:,['ref_coord','qry_coord']].astype(int) # convert numeric coordinates to int
  return anchors
  
def get_scaling_factor(genome_size_reference, half_life_distance):
  # function to determine scaling factor that produces a score of 0.5 for a given half_life_distance in the reference species.
  # this scaling factor will be used in all other species in the graph, but scaled to the according respective genome sizes.
  return -int(half_life_distance) / (genome_size_reference * np.log(0.5))

def projection_score(x, anchors, genome_size, scaling_factor=5e-4):
  # anchors must be the locations of the up- and downstream anchors, not the data frame with ref and qry coordinates.
  # the scaling factor determines how fast the function falls when moving away from an anchor.
  # ideally, we define a half-life X_half, i.e. at a distance of X_half, the model is at 0.5.
  # with a scaling factor of 50 kb, X_half is at 20 kb (with 100 kb at 10 kb)
  d = min([abs(x-y) for y in anchors])
  return np.exp(-d / (genome_size * scaling_factor))

### the function takes the input from shortest_path[ref] and returns the values to put into orange
def project_genomic_location(ref, qry, ref_coords, score, pwaln, genome_size, scaling_factor):
  # if not qry in pwaln[ref].keys(): # immediately return in case of a missing pairwise alignment
  #   return 0, np.nan, np.nan, np.nan
  ref_chrom = int(ref_coords.split(':')[0])
  ref_loc = int(ref_coords.split(':')[1])
  anchors = get_anchors(pwaln[ref][qry], ref_chrom, ref_loc)
  if anchors.shape[0] < 2: # if only one anchor is found because of border region, return 0 score and empty coordinate string
    return 0., '', ',', ','
  ref_anchors = ','.join(anchors.apply(lambda x: '{}:{}'.format(x['ref_chrom'], x['ref_coord']), axis=1))
  qry_anchors = ','.join(anchors.apply(lambda x: '{}:{}'.format(x['qry_chrom'],x['qry_coord']), axis=1))
  x_relative_to_upstream = (ref_loc - anchors.ref_coord['upstream']) / max(np.diff(anchors.ref_coord)[0], 1) # the max statement prevents potential zero division when anchors are the same (i.e. when the coord is on an alignment)
  qry_loc = int(anchors.qry_coord['upstream'] + np.diff(anchors.qry_coord)[0] * x_relative_to_upstream)
  qry_chrom = anchors.qry_chrom['upstream']
  # ONLY USE DISTANCE TO CLOSE ANCHOR AT REF SPECIES, because at the qry species it should be roughly the same as it is a projection of the reference.
  score *= projection_score(ref_loc, anchors.ref_coord, genome_size[ref], scaling_factor)
  qry_coords = '{}:{}'.format(qry_chrom, qry_loc)
  return score, qry_coords, ref_anchors, qry_anchors

def get_shortest_path_to_qry(target_species, shortest_path):
  # Backtrace the shortest path from the reference to the given target species.
  current_species = target_species
  l = []
  while current_species: # x='' at ref species, condition gets False, loop stops
    l.append(current_species)
    current_species = shortest_path[current_species][1]
    
  res = pd.DataFrame({k : shortest_path[k] for k in l[::-1]}, index=['score','from','coords','ref_anchors', 'qry_anchors']).T.loc[:,['from','score','coords','ref_anchors', 'qry_anchors']]
  # res.ref_anchors = res.ref_anchors.apply(lambda x: ','.join(x)) # reformat tuple of anchor-pair to a comma-separated string
  # res.qry_anchors = res.qry_anchors.apply(lambda x: ','.join(x)) # reformat tuple of anchor-pair to a comma-separated string  
  return res

def get_shortest_path(ref, qry, ref_coords, species, pwaln, genome_size, scaling_factor, verbose=False):
  shortest_path = {}
  orange = []
  heapq_max.heappush_max(orange, (1, ref, ref_coords))
  shortest_path[ref] = (1.0, '', ref_coords, ',', ',')

  while len(orange) > 0:
    (current_score, current_species, current_coords) = heapq_max.heappop_max(orange)
    if shortest_path.get(current_species,(0,))[0] > current_score:
        continue # the current species was already reached by a faster path, ignore this path and go to the next species
    if verbose:
      print('current species and score: ', current_species, current_score) # remember: this is not necessarily going to the shortest path as it might be a dead end that doesn't lead to the qry. not all printed species are part of the shortest path!
    if current_species == qry:
      break # qry species reached, stop
    for nxt_species in pwaln[current_species].keys():
      nxt_best_score = shortest_path.get(nxt_species,(0,))[0] # current score entry for nxt_species in shortest_path
      if current_score <= nxt_best_score:
        continue # if the score to current_species was lower than any previous path to nxt_species, nxt_species won't be reached faster through current_species. ignore and move on to the next species
      nxt_score, nxt_coords, current_anchors, nxt_anchors = project_genomic_location(current_species, nxt_species, current_coords, current_score, pwaln, genome_size, scaling_factor)
      if nxt_score <= nxt_best_score:
        continue # only save the current path to nxt_species if it was indeed faster than any previous path to it
      shortest_path[nxt_species] = (nxt_score, current_species, nxt_coords, current_anchors, nxt_anchors)
      heapq_max.heappush_max(orange, (nxt_score, nxt_species, nxt_coords))
  shortest_path_to_qry = get_shortest_path_to_qry(qry, shortest_path)
  return shortest_path_to_qry, shortest_path, orange
