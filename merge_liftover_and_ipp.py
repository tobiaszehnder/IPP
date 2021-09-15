#! /usr/bin/env python

import numpy as np, pandas as pd, sys

if not len(sys.argv)==8:
    print('Usage: python merge_liftover_and_ipp.py reference target outdir bedfile liftover_projections IPP_projections IPP_unmapped')
    sys.exit(1)
_, reference, target, outdir, bedfile, liftover_projections, ipp_projections, ipp_unmapped = sys.argv

# read input regions
input_regions = pd.read_csv(bedfile + '_centered', sep='\t', header=None, index_col=3, usecols=range(4))
input_regions.columns = ['chrom','start','end']

# read liftover projections
liftover = pd.read_csv(liftover_projections, sep='\t', header=None, index_col=3, usecols=range(4))
liftover.columns = ['chrom','start','end']
liftover.insert(0, 'coords_ref', input_regions.loc[liftover.index].apply(lambda x: '%s:%s' %(x['chrom'], x['start']), axis=1))
liftover.insert(1, 'coords_target', liftover.apply(lambda x: '%s:%s' %(x['chrom'], x['start']), axis=1))
liftover.insert(2, 'score', 1)
liftover = liftover.loc[:,['coords_ref', 'coords_target', 'score']]

# read IPP projections
proj = pd.read_csv(ipp_projections, sep='\t', index_col=0)
proj = proj.loc[:,['coords_ref', 'coords_multi', 'score_multi']]
proj.columns = ['coords_ref', 'coords_target', 'score']
try:
    unmapped = pd.read_csv(ipp_unmapped, header=None)[0].values
except pd.errors.EmptyDataError:
    unmapped = []
order = [x for x in input_regions.index.values if x not in unmapped]

# concatenate liftover and IPP projections to a table with all reference and target coordinates and scores
df = pd.concat([proj, liftover]).reindex(index=order)

# write target coordinates to a bed file
outfile = '%s/%s' %(outdir, bedfile.split('/')[-1].replace('.bed', '.%s.bed' %target))
df_target_bed = df.loc[:,['score']]
df_target_bed.insert(0, 'name', df_target_bed.index)
df_target_bed['strand'] = '.'
df_target_bed.insert(0, 'chrom', [x.split(':')[0] if isinstance(x,str) else np.nan for x in df.coords_target])
df_target_bed.insert(1, 'start', [int(x.split(':')[1])-1 if isinstance(x,str) else np.nan for x in df.coords_target])
df_target_bed.insert(2, 'end', [int(x.split(':')[1]) if isinstance(x,str) else np.nan for x in df.coords_target])
df_target_bed = df_target_bed.astype({'start': int, 'end': int})
df_target_bed.to_csv(outfile, sep='\t', index=False, header=False)
