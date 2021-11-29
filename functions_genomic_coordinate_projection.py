import numpy as np, pandas as pd, sys, os, glob, heapq, heapq_max, re, itertools, copy
# if 'functions' in sys.modules:
#   del sys.modules['functions'] # deletes previously imported module so that potential changes will be loaded
# from functions import *

# function for computing cosine similarity of 8-mers of reference vs all query projections
def compute_cosine_similarity_ref_qry(seqs):
  
  # compute possible kmers
  ks = [4,6,8]
  possible_kmers = [''.join(p) for k in ks for p in itertools.product(list('ACGT'), repeat=k)]
  kmer_ids = {kmer : get_kmer_id(kmer) for kmer in possible_kmers}
  
  # compute repeat_score for k 4, 6 and 8 and take mean
  # here, taking multiple k is somewhat important to reflect repetitiveness on different scales (e.g. ATATATATATAT vs ACCGTTACCGTTACCGTT)
  rep_scores = pd.DataFrame({i : {j : np.mean([compute_repeat_score(seqs.loc[i,j], possible_kmers, kmer_ids, k=k) for k in ks]) for j in seqs.columns} for i in seqs.index}).T
  rep_scores = rep_scores.loc[:,seqs.columns] # put columns in right order
  
  # count kmers
  k = 8
  kmer_counts_dict_ref = {i : count_kmers(seqs.ref[i], kmer_ids, ks=[8]) for i in seqs.index}
  kmer_counts_dict_qry = {i : {j : count_kmers(seqs.loc[i,j], kmer_ids, ks=[8]) for j in seqs.columns[1:]} for i in seqs.index}
  
  # compute similarity scores
  sim = np.array([[compute_similarity(kmer_counts_dict_ref[i], kmer_counts_dict_qry[i][j], 'cosine') for j in seqs.columns[1:]] for i in seqs.index])
  sim = pd.DataFrame(sim, columns=seqs.columns[1:])
  return sim, rep_scores

# precompute poisson_estimate for computation of repeat_score
def compute_poisson_estimate(seq_len, k):
  nuc_freq = .25
  try:
    return ((seq_len-int(k)+1) * nuc_freq**int(k))**2
  except ValueError:
    return np.nan

def compute_repeat_score(seq, possible_kmers, kmer_ids, k):
  l = len(seq)
  try:
    d = count_kmers(seq, kmer_ids, [int(k)])
    return sum([x * (x-1) / 2. - compute_poisson_estimate(l,k) for x in d.values()]) / (l * (l-1) / 2.)
  except ValueError:
    return np.nan

def read_genome_size(filename):
  return float(pd.read_csv(filename, sep='\t', usecols=[1]).values.sum())

def get_ce_path(cne_dir, s1, s2):
  try:
    return glob.glob(cne_dir + '/ce_*' + s1 + '*' + s2 + '*')[0]
  except IndexError:
    return None

def read_ce_file(filename):
  with open(filename, 'r') as f:
    lines = f.readlines()
  df = pd.DataFrame(np.array([l.strip().split('\t')[:4] for l in lines]), columns=['ref_chrom','ref_start','ref_end','qry'])
  df = df.astype({'ref_start': int, 'ref_end': int})
  df.insert(1, 'ref_center', ((df.ref_start + df.ref_end) / 2).astype(int))
  df['qry_chrom'] = df.qry.apply(lambda x: x.split(':')[0])
  df['qry_center'] = df.qry.apply(lambda x: np.array(x.split(':')[1].split('-')).astype(int).sum() / 2).astype(int)
  return(df.loc[:,['ref_chrom', 'ref_center', 'qry_chrom', 'qry_center']])

def read_cne(ref, qry, grb):
  filename = '/project/wig/tobias/reg_evo/data/CNEs/CNEr/cne_%s_%s_35_50.bed' %(ref,qry)
  df = pd.read_csv(filename, sep='\t', header=None, names=['chrom','start','end'], usecols=range(3))
  df['center'] = pd.Series((df.start + df.end) / 2, dtype=int)
  df = df.loc[(df.chrom == grb.chrom) & (df.center > grb.start) & (df.center < grb.end), ('chrom', 'center')]
  df['qry'] = qry
  return(df)
  
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

def get_anchors(df, chrom, x, return_top_n=False): # return_top_n: do not just return the two closest anchors, but the list of topn anchors (used in get_anchors_between_boundaries())
  # if x lies within an alignment, return its up- and downstream ends as anchors
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
  ref_chrom = ref_coords.split(':')[0]
  ref_loc = int(ref_coords.split(':')[1])
  anchors = get_anchors(pwaln[ref][qry], ref_chrom, ref_loc)
  if anchors.shape[0] < 2: # if only one anchor is found because of border region, return 0 score and empty coordinate string
    return 0., '', ',', ','
  ref_anchors = ','.join(anchors.apply(lambda x: x['ref_chrom'] + ':' + str(x['ref_coord']), axis=1))
  qry_anchors = ','.join(anchors.apply(lambda x: x['qry_chrom'] + ':' + str(x['qry_coord']), axis=1))
  x_relative_to_upstream = (ref_loc - anchors.ref_coord['upstream']) / max(np.diff(anchors.ref_coord)[0], 1) # the max statement prevents potential zero division when anchors are the same (i.e. when the coord is on an alignment)
  qry_loc = int(anchors.qry_coord['upstream'] + np.diff(anchors.qry_coord)[0] * x_relative_to_upstream)
  qry_chrom = anchors.qry_chrom['upstream']
  # ONLY USE DISTANCE TO CLOSE ANCHOR AT REF SPECIES, because at the qry species it should be roughly the same as it is a projection of the reference.
  score *= projection_score(ref_loc, anchors.ref_coord, genome_size[ref], scaling_factor)
  qry_coords = qry_chrom + ':' + str(qry_loc)
  return score, qry_coords, ref_anchors, qry_anchors

def get_shortest_path_to_qry(x, shortest_path):
  l = [x]
  while x: # x='' at ref species, condition gets False, loop stops
    x = shortest_path[x][1]
    if x:
      l.append(x)
  res = pd.DataFrame({k : shortest_path[k] for k in l[::-1]}, index=['score','from','coords','ref_anchors', 'qry_anchors']).T.loc[:,['from','score','coords','ref_anchors', 'qry_anchors']]
  # res.ref_anchors = res.ref_anchors.apply(lambda x: ','.join(x)) # reformat tuple of anchor-pair to a comma-separated string
  # res.qry_anchors = res.qry_anchors.apply(lambda x: ','.join(x)) # reformat tuple of anchor-pair to a comma-separated string  
  return res

def get_shortest_path(ref, qry, ref_coords, species, pwaln, genome_size, scaling_factor, verbose=False):
  if verbose:
    print('current species: (might be a dead end)')
  shortest_path = {}
  orange = []
  heapq_max.heappush_max(orange, (1, ref, ref_coords))
  shortest_path[ref] = (1.0, '', ref_coords, ',', ',')

  while len(orange) > 0:
    (current_score, current_species, current_coords) = heapq_max.heappop_max(orange)
    if shortest_path.get(current_species,(0,))[0] > current_score:
        continue # the current species was already reached by a faster path, ignore this path and go to the next species
    if verbose:
      print(current_species, current_score) # remember: this is not necessarily going to the shortest path as it might be a dead end that doesn't lead to the qry. not all printed species are part of the shortest path!
    if current_species == qry:
      break # qry species reached, stop
    # for nxt_species in species[species!=current_species]:
    for nxt_species in pwaln[current_species].keys(): # this allows incomplete pwalns where some pairwise alignments are missing
      nxt_best_score = shortest_path.get(nxt_species,(0,))[0] # current score entry for nxt_species in shortest_path
      if current_score <= nxt_best_score:
        continue # if the score to current_species was lower than any previous path to nxt_species, nxt_species won't be reached faster through current_species. ignore and move on to the next species
      else:
        nxt_score, nxt_coords, current_anchors, nxt_anchors = project_genomic_location(current_species, nxt_species, current_coords, current_score, pwaln, genome_size, scaling_factor)
      if nxt_score <= nxt_best_score:
        continue # only save the current path to nxt_species if it was indeed faster than any previous path to it
      else:
        shortest_path[nxt_species] = (nxt_score, current_species, nxt_coords, current_anchors, nxt_anchors)
        heapq_max.heappush_max(orange, (nxt_score, nxt_species, nxt_coords))
  shortest_path_to_qry = get_shortest_path_to_qry(qry, shortest_path)
  return shortest_path_to_qry, shortest_path, orange


def get_anchors_between_boundaries(pwaln, current_species, next_species, ob, reference_species, species, rel_ori, verbose=False):
  # this function returns all anchors between two genomic coordinates of a ref to a qry species
  # new: added orientation check: local small inversion rarely give single or few alignments that are oriented the 'wrong' way. discard them. only keep the major_orientation
  # new: instead of only looking at anchors in OB, first get topn closest anchors to the center of the OB from get_anchors() which includes a major_chrom and collinearity check, then filter for anchors in OB-span
  # this better detects local inversions and especially outliers. the disadvantage is that it could ignore true local inversions, but then the species will probably just not be used as a bridge which is better than including wrong anchors.
  if not current_species in species: # in case the current species was eliminated before from the species list because it had contradicting orientation to reference and target (i.e. local inversions)
    return pd.DataFrame(), species, rel_ori
    # get topn anchors from get_anchors() (usually less than top n because of the major_chrom and collinearity check in get_anchors())
  anchor_cols = ['ref_chrom','ref_start','ref_end','qry_chrom','qry_start','qry_end']
  ob_center = int(np.mean(ob[current_species].current_coord))
  df = get_anchors(pwaln[current_species][next_species], ob[current_species].current_chrom[0], ob_center, return_top_n=True)[anchor_cols]
  
  # orientation check: local small inversion rarely give single or few alignments that are oriented the 'wrong' way. discard them. only keep the major_orientation
  # do this check here instead of after the qry-OB-check when the number of anchors might be down to very few
  if df.shape[0] > 0:
    same_orientation = (df.qry_start - df.qry_end) < 0
    majority_has_same_orientation = same_orientation.value_counts().idxmax()
    df = df.loc[same_orientation==majority_has_same_orientation]
    # save relative orientation and remove the current (ref) species from the species list if the relative orientation of current and target contradict the relative orientations of reference to current and to target
    rel_ori = '+' if majority_has_same_orientation else '-'
    strand_sign = {'+':1, '-':-1}
    if not current_species == reference_species and not strand_sign[rel_ori] * strand_sign[rel_ori] == strand_sign[rel_ori]:
      species = species[~np.in1d(species,[current_species])]
      if verbose:
        print('Removed %s in get_anchors_between_boundaries() (relative orientation current-%s contradict reference-%s and reference-%s)' %(current_species,next_species,current_species,next_species))
      return pd.DataFrame(), species, rel_ori
  
  # filter for anchors within ref's OB
  df = df.loc[(df.ref_chrom == ob[current_species]['current_chrom'][0]) & (df.ref_end >= ob[current_species].loc['upstream','current_coord']) & (df.ref_start <= ob[current_species].loc['downstream','current_coord'])]
  
  # filter for anchors within qry's OB
  if df.shape[0] > 0:
    if df.qry_start.values[0] < df.qry_end.values[0]: # if qry is inverted, check for inverted start/end being located within OB span.
      df = df.loc[(df.qry_chrom == ob[next_species]['current_chrom'][0]) & (df.qry_end >= ob[next_species].loc['upstream','current_coord']) & (df.qry_start <= ob[next_species].loc['downstream','current_coord']), :]
    else:
      df = df.loc[(df.qry_chrom == ob[next_species]['current_chrom'][0]) & (df.qry_start >= ob[next_species].loc['upstream','current_coord']) & (df.qry_end <= ob[next_species].loc['downstream','current_coord']), :]
  ### do I need to check if the chromosome is the same as the one the direct zf-mm anchors point to? if so, only for the target?
  ### if I don't, few helpful anchors could be discarded when the number of anchors pointing to outliers is bigger. but probably this is negligible anyways.
  # collinearity check (new updated function is much faster (< 100 ms for 1000 coordinates) so I can check them all for collinearity)
  df = df.iloc[np.intersect1d(df.qry_start.values, longest_sorted_subsequence(df.qry_start.values.astype(int)), return_indices=True)[1]]
  
  df.insert(0, 'ref', current_species)
  df.insert(5, 'qry', next_species)
  return df, species, rel_ori

def set_ob(new_boundary, sp, d, ob, verbose=False, silent=False):
  # this function takes a given coordinate and sets it as the new boundary if it lies within the previously defined ob and thus decreases the span of the boundaries.
  # If not, the currently stored boundary is returned.
  # Run this function separately for each direction
  directions = ['upstream','downstream']
  opposite_direction = {d[0]: d[1] for d in (directions,directions[::-1])}
  current_boundary = ob.get(sp, pd.DataFrame(index=[d]))
  if current_boundary.shape[1] == 0:
    ob[sp] = new_boundary.loc[[d]] # return the new boundary if no value is set yet
    return ob
  else:
    # test if new boundary is on the correct side of a potentially existing opposite-direction boundary (new upstream should be smaller than old downstream boundary, otherwise keep old upstream, and vice versa)
    if opposite_direction[d] in current_boundary.index:
      new_coord = new_boundary.current_coord[0]
      old_opposite_coord = current_boundary.loc[opposite_direction[d],'current_coord']
      correct_side = new_coord <= old_opposite_coord if d == 'upstream' else new_coord >= old_opposite_coord
      if not correct_side:
        if verbose:
          print('New %s boundary (%s) in %s lies outside current boundaries (%s, %s). Kept old values.' %(d,new_coord,sp,*current_boundary['current_coord']))
        return ob
    if not opposite_direction[d] in current_boundary.index:
      raise ValueError('There is both a current and a new value for the %s boundary, but no current value for the %s boundary. Set that first.' %(d, opposite_direction[d]))
    df = pd.concat([new_boundary, current_boundary], sort=False).loc[[d]] # new boundary first because at the beginning there is no old boundary, but there always is a new one. that way, idx=1 is always the old.
    idx = np.argmin(abs(df.current_coord.values - current_boundary.loc[[opposite_direction[d]],'current_coord'].values))
    if idx == 0 and not silent:
      print('New %s boundary in target set at %s' %(d,new_boundary.loc[d,'current_coord']))
    if idx == 1 and verbose:
      print('New %s boundary (%s) is not decreasing the %s anchor span. Kept old value (%s).' %(d,new_boundary.loc[d,'current_coord'],sp,current_boundary.loc[d,'current_coord']))
    ob[sp] = pd.concat([current_boundary.loc[[opposite_direction[d]]], df.iloc[[idx]]], axis=0, sort=False).loc[directions] # identify if the old or the new boundaries are decreasing the span and return that row
    return ob

def get_ob(ob, reference_species, reference_coords, current_species, current_coords, target, species, pwaln, relative_orientation, verbose=False):
  # function for finding outer boundaries for both directions as the direct anchors to the target species. this function is run in the beginning to find the first set of OBs for every species.
  directions = ['upstream','downstream']
  opposite_direction = {d[0]: d[1] for d in (directions,directions[::-1])}
  if not current_species == target:
    try:
      direct_anchors = pd.concat([get_anchors(pwaln[current_species][target], current_coords[d]['chrom'], current_coords[d]['coord']).loc[d] for d in directions], axis=1).T
    except KeyError: # in this case, the current species does not have both up-/ and downstream anchors directly to mouse and will be removed from the species list.
      species = species[~np.in1d(species,[current_species])]
      if verbose:
        print('Removed %s in get_ob() (does not have 2 direct anchors to target)' %current_species)
      return ob, species, relative_orientation
    # check if both anchors have the same orientation. if not, the position is likely at a synteny break position. in that case, remove the current species from the species list
    if len(np.unique(direct_anchors.qry_strand)) == 2:
      species = species[~np.in1d(species,[current_species])]
      if verbose:
        print('Removed %s in get_ob() (anchors have different orientation)' %current_species)
      return ob, species, relative_orientation
    # store information about the relative orientation of the current species to the target
    relative_orientation.loc[current_species,target] = relative_orientation.loc[target,current_species] = direct_anchors.qry_strand[0]
    # remove the current species from the species list if the relative orientation of current and target contradict the relative orientations of reference to current and to target
    strand_sign = {'+':1, '-':-1}
    if not reference_species == 'origin' and not strand_sign[relative_orientation.loc[reference_species,current_species]] * strand_sign[relative_orientation.loc[current_species,target]] == strand_sign[relative_orientation.loc[reference_species,target]]:
      if verbose:
        print('Removed %s in get_ob() (relative orientation current-target contradict reference-current and reference-target)' %current_species) 
      species = species[~np.in1d(species,[current_species])]
      return ob, species, relative_orientation
    # set the current species boundaries
    for d in directions:
      # new_boundary contains the ref-coord of the current->target anchor and the coordinate of the reference species that lead to the current species.
      # note that the ref-coord of the current->target anchor always lies more outwards than the reference coordinate.
      # the latter is only saved for later that I can back-track how I reached the current species, i.e. to identify the final reference anchors involved in the path to the most narrow target anchor span.
      if reference_coords == 'origin':
        new_boundary_prev_coords = pd.DataFrame({'prev_species': reference_species, 'prev_chrom':reference_coords, 'prev_coord':reference_coords}, index=[d])
#         current_d = d
      else:
        new_boundary_prev_coords = pd.DataFrame({'prev_species': reference_species, 'prev_chrom':reference_coords[d]['chrom'], 'prev_coord':reference_coords[d]['coord']}, index=[d])
        # combine upstream with downstream anchor if relative orientation is '-'
#         current_d = d if relative_orientation.loc[reference_species,current_species] =='+' else opposite_direction[d]
      new_boundary = pd.concat([new_boundary_prev_coords, direct_anchors.loc[[d]]], axis=1).rename(columns={'ref_chrom':'current_chrom', 'ref_coord':'current_coord'})
      ob = set_ob(new_boundary, current_species, d, ob, silent=True)
    new_target_boundaries = orient_anchors(pd.concat([pd.DataFrame({'prev_species':current_species}, index=directions), direct_anchors], axis=1).rename(columns={'ref_chrom':'prev_chrom', 'ref_coord':'prev_coord','qry_chrom':'current_chrom', 'qry_coord':'current_coord'}), coord_column='current_coord')
    # switch up-/downstream labels and qry_start/qry_end in case the current-target orientation is '-' (analogous to orient_anchors()). otherwise the target will have inverted coordinates.
    if relative_orientation.loc[current_species,target] == '-':
      new_target_boundaries[['qry_start','qry_end']] = new_target_boundaries[['qry_end','qry_start']].values
#       new_target_boundaries.rename(index={'upstream':'downstream', 'downstream':'upstream'}, inplace=True)
    for d in directions:
      ob = set_ob(new_target_boundaries.loc[[d]], target, d, ob, silent=True)
  return ob, species, relative_orientation

def orient_anchors(anchors, coord_column='coord', reverse_orient=False):
  # This function checks anchors if they are inverted. and switches the indices (upstream to downstream and vice versa).
  directions = ['upstream','downstream']
  if anchors.loc['upstream',coord_column] - anchors.loc['downstream',coord_column] > 0:
    anchors.index = anchors.index[::-1]
  if reverse_orient: # reverse orient
    anchors.index = anchors.index[::-1]
  return(anchors.loc[directions])

def get_next_anchor(available_species, df, x, direction, verbose=False):
  # this function gives back the closest anchor stored in df in the given direction
  # in upstream direction it is the end-coord that needs to be closest, in downstream direction it is the start-coord.
  # new: include the possibility that the current position x is actually ON an alignmentdf = anchors[current_species]
  # note that alignments are not necessarily the same length in ref and qry (not exact matches). thus I can't take the exact points.
  # instead, go to the outer coordinate (upstream: start, downstream: end).
  # if there are multiple alignments, take the one with the outer coordinate the closest.
  # once chosen, check if there is no non-overlapping anchor that is closer to the outer coordinate. if there is, ignore the alignment and take the non-overlapping anchor (set a flag to know)
  if df.shape == (0,0): # return empty data frame if there are no anchors
    return pd.DataFrame()
  df = df.loc[df.qry.apply(lambda x: x in available_species)]
  ov_aln = df.loc[(df.ref_start < x) & (df.ref_end > x)].copy()
  ov_aln.loc[:,'opposite_strand'] = [True if span < 0 else False for span in np.diff(ov_aln[['qry_start','qry_end']])]
  ov_aln_sorted = pd.DataFrame()
  if ov_aln.shape[0] > 0: # check if x overlaps an alignment and take the one where the outer coordinate is closest (upstream: start, downstream: end)
    if direction == 'upstream':
      ov_aln_sorted = ov_aln.sort_values(by='ref_start', ascending=False)
      # resize to 1 bp width on outer coordinate because that is the only thing I can confidently do. everything else (interpolation) might be a move in the inward direction.
      ov_aln_sorted.ref_end = ov_aln_sorted.ref_start
      ov_aln_sorted.qry_end = ov_aln_sorted.qry_start
    else:
      ov_aln_sorted = ov_aln.copy().sort_values(by='ref_end', ascending=True)
      ov_aln_sorted.ref_start = ov_aln_sorted.ref_end
      ov_aln_sorted.qry_start = ov_aln_sorted.qry_end
        
  # check if the outer coordinate of the alignment is still closer than the inner coordinate of the next anchor. if not, set ignore_aln flag to True and take the next non-overlapping anchor.
  if direction == 'upstream':
    next_non_ov_anchors = df.loc[df.ref_end < x,:].sort_values(by='ref_end', ascending=False)
  else:
    next_non_ov_anchors = df.loc[df.ref_start > x,:].sort_values(by='ref_start', ascending=True)
  next_non_ov_anchors.loc[:,'opposite_strand'] = [True if span < 0 else False for span in np.diff(next_non_ov_anchors[['qry_start','qry_end']])]
  if next_non_ov_anchors.shape[0] == 0 and ov_aln_sorted.shape[0] == 0: # in case no anchor was found at all
    res = pd.Series()
  elif next_non_ov_anchors.shape[0] > 0 and ov_aln_sorted.shape[0] == 0: # in case only non-overlapping anchors were found
    res = next_non_ov_anchors.iloc[0]
  elif next_non_ov_anchors.shape[0] == 0 and ov_aln_sorted.shape[0] > 0: # in case only an overlapping alignment was find
    res = ov_aln_sorted.iloc[0]
  else: # both overlapping and non-overlapping anchors found. only take best overlapping alignment if its outer coordinate is closer to x than the inner coordinate of best non-overlapping
    # (we can actually both check for inner coordinate because the potential overlapping alignment's outer coordinate was set to the values of its inner coordinate before)
    if ((direction == 'upstream') and (ov_aln_sorted.iloc[0]['ref_end'] < next_non_ov_anchors.iloc[0]['ref_end'])) or ((direction == 'downsteam') and (ov_aln_sorted.iloc[0]['ref_start'] > next_non_ov_anchors.iloc[0]['ref_start'])): # next non-overlapping anchor is closer than outer coord of overlapping alignment
      if verbose:
        inner_anchor_coord = next_non_ov_anchors.iloc[0]['ref_end'] if direction == 'upstream' else next_non_ov_anchors.iloc[0]['ref_start']
        print('Outer coordinate of overlapping alignment: %s %s-%s to %s %s-%s is further than inner coordinate of next non-overlapping anchor (%s). Returned the latter.' %(tuple(ov_aln.loc[ov_aln_sorted.iloc[0].name,['ref','ref_start','ref_end','qry','qry_start','qry_end']])+(inner_anchor_coord,)))
      res = next_non_ov_anchors.iloc[0]
    else: # the overlapping alignment's outer coords are still closer than the next non-overlapping anchor
      if verbose:
        print('Overlapping alignment: %s %s-%s to %s %s-%s' %tuple(ov_aln.loc[ov_aln_sorted.iloc[0].name,['ref','ref_start','ref_end','qry','qry_start','qry_end']])) # print the original alignment coords, not the resized
      res = ov_aln_sorted.iloc[0]
  
  return res

def is_outside_boundary(x, ob, current_species):
  # this function returns True if x is between the boundaries stored in ob for a particular species
  return x not in np.concatenate([np.arange(*ob[current_species]['current_coord'].values), np.arange(*ob[current_species]['current_coord'].values,-1)])[1:-1] # [1:-1] excludes boundaries, i.e. when x reached the boundary it is considered 'outside'

def move_through_anchors(available_species, reference_species, prev_species, prev_chrom, prev_coord, prev_direction, current_species, current_chrom, current_coord, direction, target, anchors, ob, relative_orientation, verbose=False, silent=False):
  directions = ['upstream','downstream']
  opposite_direction = {d[0]: d[1] for d in (directions,directions[::-1])}
  # set local outer boundary to what is currently stored in global OB, then update the global OB to the current position.
  # that way, when we get back to this species at some point, we'll know to continue until the originally set boundary (local).
  # however, if the current species is reached from anywhere else, it will immediately abort if it arrives outside this newly set boundary.
  try:
    ob_local = ob[current_species].loc[direction,'current_coord']
  except KeyError:
    print(current_chrom, current_coord)
    raise KeyError()
  ob[current_species].loc[direction,['prev_species','prev_chrom','prev_coord','current_chrom','current_coord']] = prev_species, prev_chrom, prev_coord, current_chrom, current_coord
  inside_boundary = ob_local < current_coord if direction == 'upstream' else ob_local > current_coord
  if verbose:
    print('Local boundary %s %s: %s' %(current_species,direction,ob_local))
    print('Global boundary %s %s: %s' %(current_species,direction,ob[current_species].loc[direction,'current_coord']))
    print('Inside boundary: %s' %inside_boundary)
  while inside_boundary:
    if verbose:
      print('Current position: %s %s %s' %(current_species, current_coord, direction))
    next_anchor = get_next_anchor(available_species, anchors[current_species], current_coord, direction, verbose=verbose)
    if next_anchor.shape[0] == 0:
      if verbose:
        print('No next anchor found in %s.' %current_species)
      break
    if verbose:
      print('Next anchor at: %s %s %s' %(next_anchor.ref, next_anchor.ref_start, next_anchor.ref_end))
    current_coord = next_anchor.ref_end if direction == 'upstream' else  next_anchor.ref_start # set current position to the reached anchor
    # set next_species to the one the anchor points to. if it is mouse, update boundaries, otherwise call this function recursively with next_species.
    # important: it is the CURRENT direction that determines if we want the start or the end coord of the qry anchor.
    # reason: if the next species is inverted, start and end coordinates will be switched when running move_through_anchors() with the next species (i.e. we will invert our own view).
    next_species = next_anchor.qry
    next_direction = direction if next_anchor.opposite_strand == False else opposite_direction[direction]
    next_chrom = next_anchor.qry_chrom
    next_coord = next_anchor.qry_end if direction == 'upstream' else next_anchor.qry_start
    if verbose:
      print('Next anchor points to: %s %s %s' %(next_species, next_coord, next_direction))
    # discard the next anchor if it points outside of OB in the next species and move 1 position in the current direction
    if is_outside_boundary(next_coord, ob, next_species):
      if verbose:
        print('Anchor to %s points outside the boundaries. Staying in %s and moving 1 %s' %(next_species,current_species,direction))
      if direction == 'upstream':
        current_coord -= 1
      else:
        current_coord += 1
      continue
    # if it points to the target, save boundary and break the loop because we'll not find a better anchor from this species to mouse via other species.
    elif next_species == target:
      new_target_boundary = pd.DataFrame({'prev_species':current_species,'prev_chrom':current_chrom,'prev_coord':current_coord,'current_chrom':next_chrom,'current_coord':next_coord}, index=[next_direction])
      if not silent:
        print('Reached target from %s at %s' %(current_species,next_coord))
      ob = set_ob(new_target_boundary, target, next_direction, ob, verbose=verbose, silent=silent)
      break
    else:
      next_available_species = np.setdiff1d(available_species, current_species) # keeps track of previously visited species and removes them from the list to prevent loops (happens rarely, I guess for non-syntenic regions)
      ob, relative_orientation = move_through_anchors(next_available_species, reference_species, current_species, current_chrom, current_coord, direction, next_species, next_chrom, next_coord, next_direction, target, anchors, ob, relative_orientation, verbose, silent)
  if prev_species == 'origin' and not silent:
    print('Reached %s boundary of reference species. Done.\n' %direction)
  return ob, relative_orientation

def walk_path(ob, to_species, to_direction, to_strand, relative_orientation, strand_sign):
  # this function is used in get_anchor_path() to walk backwards to the next edge in the path
  df = ob[to_species].loc[[to_direction],['prev_species','prev_chrom','prev_coord','current_chrom','current_coord']]
  df.columns = ['from_species','from_chrom','from_coord','to_chrom','to_coord']
  df.insert(3,'to_species',to_species)
  df.insert(6,'to_strand',to_strand)
  from_strand = strand_sign[relative_orientation.loc[df.from_species[0],to_species]] * to_strand if not df.from_species[0] == 'origin' else 1
  df.insert(3,'from_strand',from_strand)
  return df

def get_anchor_path(direction, ob, reference, target, relative_orientation):
  directions = ['upstream','downstream']
  opposite_direction = {d[0]: d[1] for d in (directions,directions[::-1])}
  strand_sign = {'+':1, '-':-1}
  strand_sign_inv =  {v: k for k,v in strand_sign.items()}
  to_species = target
  to_strand = strand_sign[relative_orientation.loc[reference,target]]
  to_direction = direction if to_strand == 1 else opposite_direction[direction]
  path = walk_path(ob, to_species, to_direction, to_strand, relative_orientation, strand_sign)
  while not path.from_species.values[-1] == 'origin':
    to_species, to_strand = path.iloc[-1][['from_species','from_strand']].values
    to_direction = direction if to_strand == 1 else opposite_direction[direction]
    path = pd.concat([path, walk_path(ob, to_species, to_direction, to_strand, relative_orientation, strand_sign)], axis=0)
  path = path.loc[::-1].replace(strand_sign_inv).reset_index(drop=True)
  path['reference_direction'] = direction
  path.set_index('reference_direction', inplace=True)
  path.columns = pd.MultiIndex.from_tuples([tuple(x.split('_')) for x in path.columns])
  return path

def propagate_anchors(reference, target, coord, coord_id, pwaln, verbose=False, silent=False, test_mode=False):
  # initialize for reference species
  if not silent:
    print('Initializing variables')
  species = np.array(list(pwaln.keys())) # The list of species is determined based on the keys of the supplied pairwise aln (pwaln) dict
  x = dict(zip(['chrom','coord'],coord.split(':')))
  x['coord'] = int(x['coord'])
  directions = ['upstream','downstream']
  opposite_direction = {'upstream':'downstream', 'downstream':'upstream'}
  anchors = {}
  relative_orientation = pd.DataFrame(index=species, columns=species)
  prev_species = 'reference'
  ob = {} # define dict with outer boundaries
  current_species = reference
  current_coords = {d: x for d in directions} # in the reference, starting point x is both the up- and downstream coordinate. every other species will first be reached through an anchor from the reference an those two coordinates will thus differ.
  ob, species, relative_orientation = get_ob(ob, 'origin', 'origin', reference, current_coords, target, species, pwaln, relative_orientation, verbose)
  # if the reference species is not in the species list anymore, that means that not both direct up- and downstream anchors were found from refernce to target.
  # in that case, return an empty path and ob. further handling of this case in ~/dev/propagate_anchors.py (message printed)
  if not reference in species:
    path = pd.DataFrame()
  else:
    # project x to every species and from the coordinates of those anchors get OBs
    if not silent:
      print('Setting outer boundaries')
    for current_species in species:
      if not current_species in [reference,target]:
#         if current_species=='cteIde1:'
        direct_anchors_from_reference = get_anchors(pwaln[reference][current_species], x['chrom'], x['coord'])
        if not direct_anchors_from_reference.shape[0] == 2: # in this case, the reference species does not have anchors in both directions from x to the current species. Remove from species list.
          species = species[~np.in1d(species,[current_species])]
          if verbose:
            print('Removed %s in propagate_anchors() (does not have 2 direct anchors from reference to %s: ' %(current_species,current_species))
          continue
        else:
          relative_orientation.loc[reference,current_species] = relative_orientation.loc[current_species,reference] = direct_anchors_from_reference.qry_strand[0]
          current_coords = orient_anchors(direct_anchors_from_reference.rename(columns={'qry_chrom':'chrom', 'qry_coord':'coord'}).loc[:,['chrom','coord']], coord_column='coord').T.to_dict()
          # reverse_orient reference coords if relative orientation is '-' so that the correct from-coordinate is saved in ob in get_ob()
          reverse_orient = True if relative_orientation.loc[reference, current_species] == '-' else False
          reference_coords = orient_anchors(direct_anchors_from_reference.rename(columns={'ref_chrom':'chrom', 'ref_coord':'coord'}).loc[:,['chrom','coord']], coord_column='coord', reverse_orient=reverse_orient).T.to_dict() 
        ob, species, relative_orientation = get_ob(ob, reference, reference_coords, current_species, current_coords, target, species, pwaln, relative_orientation, verbose)
    ob_original = copy.deepcopy(ob)
#     return ob, relative_orientation, species

    # for every species, call all anchors to any other species between OB
    if not silent:
      print('Calling collinear anchors between outer boundaries')
    for ref in species:
      anchors[ref] = pd.DataFrame()
      for qry in species:
          if not qry in [reference,ref]:
            anchors_ref_qry, species, rel_ori = get_anchors_between_boundaries(pwaln, ref, qry, ob, reference, species, relative_orientation.loc[ref,qry], verbose)
            if not relative_orientation.loc[ref,qry] in ['+','-'] and not rel_ori == relative_orientation.loc[ref,qry]:
              # in this case the new anchor orientation contradicts what was previously saved as the relative orientation. this probably only happens in the reference-target case,
              # because there the initial orientation is determined by the direct anchors of the initial coordinate, and now we look for the anchors in relation to the OB-center
              # if the original X is close to a synteny break, this might give us two different relative orientations. in that case, remove the qry from the species list and if it is the reference, abort
              species = species[~np.in1d(species,[current_species])]
              if verbose:
                print('Removed %s in propagate_anchors() (orientation of anchors to %s between OB contradicts originally determined orientation: ' %(ref,qry))
              if ref == reference:
                print('Removed the reference species! In that case, the original coordinate likely lies at a synteny break point and this method will fail. Aborted.')
                return None, None
            relative_orientation.loc[ref,qry] = rel_ori
            if not ref in species:
              del anchors[ref]
              break
            anchors[ref] = pd.concat([anchors[ref], anchors_ref_qry])
#     return anchors, relative_orientation
    
    # move through anchors
    if not silent:
      print('Propagating anchors through species graph\n')
    for direction in directions:
      ob, relative_orientation = move_through_anchors(species, reference, 'origin', 'origin', 'origin', direction, reference, x['chrom'], x['coord'], direction, target, anchors, ob, relative_orientation, verbose=verbose, silent=silent)
    if test_mode:
      return ob, anchors, relative_orientation, ob_original
      
    # obtain optimal paths
    if not silent:
      print('Fetching optimal paths')
    path = pd.concat([get_anchor_path(direction, ob, reference, target, relative_orientation) for direction in directions], axis=0)
    path['id'] = coord_id
    path = path.reset_index().set_index(['id',*path.index.names]) 
  return ob, path


### OLD ###
  
# function to determine anchors of a given genomic coordinate
# pwaln: if x lies within an alignment, return the coordinate itself as both anchors
def get_anchors_old(df, chrom, x, return_top_n=False): # return_top_n: do not just return the two closest anchors, but the list of topn anchors (used in get_anchors_between_boundaries())
  ### NEW: overlapping alignments are no longer interpolated to the exact location on the alignment. instead, the ends of the alignments are taken as flanking anchors.
  ### this is because the alignments may have different lengths in ref and qry due to gaps and we can not be sure that the nth bp on the ref aligns with the same nth bp on the qry.
  
  # if x lies within an alignment, return its up- and downstream ends as anchors
  anchor_cols = ['ref_chrom','ref_start','ref_end','ref_coord','qry_chrom','qry_start','qry_end','qry_coord','qry_strand']
  ov_aln = df.loc[(df.ref_chrom == chrom) & (df.ref_start < x) & (df.ref_end > x),].reset_index(drop=True) # x lies in an alignment. return x itself as both anchors
  ov_aln['ref_coord'] = (ov_aln.ref_start + ov_aln.ref_end) / 2 # add the center of the alignment as the `ref_coord` column in order to check for collinearity with up- and downstream anchors later
  ov_aln['qry_coord'] = (ov_aln.qry_start + ov_aln.qry_end) / 2 # add the center of the alignment as the `qry_coord` column in order to check for collinearity with up- and downstream anchors later
  
  # first define anchors upstream downstream and ov_aln, then do major_chrom / collinearity test, then either return overlapping anchor or closest anchors.
  # take orientation into account for the anchor definition. if start > end, then the aln is to the '-' strand.
  # OLD AND WRONG COMMENT: in that scenario, if we are looking for the downstream anchor, we are interested in the smaller (WRONG) value, i.e. the end coordinate.
  # actually, we are always interested in the same labeled coordinate. so in the DOWNSTREAM case we are looking for the START coordinate.
  # but if the qry is inverted, the start coordinate will be greater than the end and so we'll save that info in the qry_strand variable later.
  # remember that the REF coords are always on the '+' strand, so for slicing the df we don't need to check for the smaller/bigger value of start/end, it will always be start < end
  # only select the first 100. the rest just takes longer to compute min / max and most likely will (and should) not be an anchor anyways. (and they are sorted by distance to x, so this is fine)
  anchors_upstream = df.loc[abs(df.loc[(df.ref_chrom == chrom) & (df.ref_end < x),].ref_end - x).sort_values().index,['ref_chrom','ref_start','ref_end','qry_chrom','qry_start','qry_end']].iloc[:100,:].reset_index(drop=True) # keeping the index makes creating a new column later very slow
  anchors_upstream.insert(3, 'ref_coord', anchors_upstream.ref_end)
  anchors_downstream = df.loc[abs(df.loc[(df.ref_chrom == chrom) & (df.ref_start > x),].ref_start - x).sort_values().index,['ref_chrom','ref_start','ref_end','qry_chrom','qry_start','qry_end']].iloc[:100,:].reset_index(drop=True) # keeping the index makes creating a new column later very slow
  anchors_downstream.insert(3, 'ref_coord', anchors_downstream.ref_start)
#   anchors_upstream.columns = anchors_downstream.columns = ['ref_chrom','ref_coord','qry_chrom','qry_start','qry_end']
  # abort if less than 5 pwalns to each side (too sparse, not able to ensure collinearity)
  if min(anchors_upstream.shape[0], anchors_downstream.shape[0]) < 5:
    return pd.DataFrame(columns=anchor_cols)
  # OLD AND WRONG COMMENT: set the corresponding start or end coordinate that is closer to the projected coordinate (choosing the max/min assures correct handling of inverted alignments)
  # WRONG: anchors_upstream['qry_coord'] = anchors_upstream.loc[:,('qry_start','qry_end')].apply(max, axis=1)
  # WRONG: anchors_downstream['qry_coord'] = anchors_downstream.loc[:,('qry_start','qry_end')].apply(min, axis=1)
  # new: choosing the SAME LABELED COORDINATE as the reference, i.e. START in case of downsteam and END in case of upstream.
  anchors_upstream['qry_coord'] = anchors_upstream.qry_end
  anchors_downstream['qry_coord'] = anchors_downstream.qry_start
  
  # top 10 produced many locally collinear pwalns that were still non-collinear outliers in the global view of the GRB (is that really true?).
  # also, top 10 still often has local cluster of outliers and is depleted in the true major_chrom pwalns, leading to the 'wrong' major_chrom.
  # since the updated collinearity test is very efficient, take top 20 to each direction.
  topn = 20
  
  # MAJOR CHROMOSOME: retain anchors that point to the majority chromosome in top twenty of both up- and downstream anchors
  try:
    major_chrom = pd.concat([anchors_upstream[:topn], ov_aln, anchors_downstream[:topn]], axis=0, sort=False).qry_chrom.value_counts().idxmax()
  except ValueError:
    print(anchors_upstream.head(2))
    print(anchors_upstream.shape)
    print(anchors_downstream.head(2))
    print(anchors_downstream.shape)
  ov_aln[ov_aln.qry_chrom == major_chrom]
  anchors_upstream = anchors_upstream[anchors_upstream.qry_chrom == major_chrom]
  anchors_downstream = anchors_downstream[anchors_downstream.qry_chrom == major_chrom]
  
  # COLLINEARITY: remove pwalns pointing to outliers by getting the longest sorted subsequence of the top n of both up- and downstream anchors.
  closest_anchors = pd.concat([anchors_upstream[:topn][::-1], ov_aln, anchors_downstream[:topn]], axis=0, sort=False).reset_index(drop=True) # reset_index necessary, otherwise working with duplicate indices messing things up
  idx_collinear = closest_anchors.index[np.intersect1d(closest_anchors.qry_coord.values, longest_sorted_subsequence(closest_anchors.qry_coord.values.astype(int)), return_indices=True)[1]] # this step takes 2 sec
  closest_anchors = closest_anchors.loc[idx_collinear,].dropna(axis=1, how='all') # drop columns if it only contains NaNs (see explanation below)
  
  # return the top anchors if the flag was passed
  if return_top_n:
    return closest_anchors[['ref_chrom','ref_start','ref_end','qry_chrom','qry_start','qry_end']]
  
  # check if the original ov_aln is still present (or ever was) in the filtered closest_anchors (that include a potential ov_aln)
  ov_aln = closest_anchors.loc[(closest_anchors.ref_start < x) & (closest_anchors.ref_end > x),].reset_index(drop=True)
  print(ov_aln)
  if ov_aln.shape[0] == 1:
    x_relative_to_upstream = x - ov_aln.ref_start[0]
    strand = '+' if ov_aln.qry_start[0] < ov_aln.qry_end[0] else '-'
    vals_up = ov_aln.copy().rename(index={0:'upstream'})
    vals_down = ov_aln.copy().rename(index={0:'downstream'})
    vals_up.ref_coord = vals_up.ref_start # here it is the start coord for same-oriented because we are looking at an overlapping alignment
    vals_down.ref_coord = vals_down.ref_end
    # for the qry_coord: always take qry_start for upstream and qry_end for downstream. why? if strand=='-', qry_start > qry_end and thus we don't need to check.
    vals_up.qry_coord = vals_up.qry_start
    vals_down.qry_coord = vals_down.qry_end
    anchors = pd.concat([pd.concat([vals_up, vals_down], axis=0), pd.DataFrame({'qry_strand':strand}, index=['upstream','downstream'])], axis=1)
      
#     # old: interpolated the actual position on the alignment. however, this is inaccurate because alignments can contain gaps and so they are not necessarily the same length in ref and qry
#     if strand == '+'
#       vals_up = [chrom, x, ov_aln.qry_chrom[0], ov_aln.qry_start[0]+x_relative_to_upstream, strand]
#       vals_down = [chrom, x, ov_aln.qry_chrom[0], ov_aln.qry_start[0]+x_relative_to_upstream, strand] # add an offset of 1 between up- and downstream anchor
#     else:
#       vals_up = [chrom, x, ov_aln.qry_chrom[0], ov_aln.qry_start[0]-x_relative_to_upstream, strand]
#       vals_down = [chrom, x, ov_aln.qry_chrom[0], ov_aln.qry_start[0]-x_relative_to_upstream, strand] # add an offset of 1 between up- and downstream anchor
#     anchors = pd.DataFrame.from_dict({'upstream': vals_up, 'downstream': vals_down}, orient='index', columns=anchor_cols)
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
