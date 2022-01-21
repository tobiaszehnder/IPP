import numpy as np, pandas as pd, pyranges as pr

def format_row_table_to_bed(row, which=('ref','qry')):
    # function to convert a row from the results table to bed format
    # the results table contains reference and query coordinates, bridging species, etc.
    # the bed format either contains the coordinates of the ref or the qry species,
    # and the coordinates of the opposite species in the name field.
    colors = {'DC':'127,201,127', 'IC':'253,180,98', 'NC':'141,153,174', 'DC+':'255,0,0', 'IC+':'255,0,0', 'DC-':'30,30,30', 'IC-':'30,30,30', 'NC+':'141,153,174', 'NC-':'141,153,174'}
    if which == 'ref':
      coords = row['coords_ref']
      name_coords = row['coords_multi']
    else:
      coords = row['coords_multi']
      name_coords = row['coords_ref']
    name = '{}_{}_{}'.format(row.name, name_coords, row['conservation'])
    bed_row = pd.Series([coords.chrom,
                         coords.loc,
                         coords.loc+1,
                         name,
                         row['score_multi'],
                         '.',
                         coords.loc,
                         coords.loc+1,
                         colors[row['conservation']]])
    return bed_row

def classify_conservation(df_projections, target_regions=pr.PyRanges(), thresh=.95, maxgap=500):
  # function for determining the conservation of sequence (DC/IC/NC) and function (+/-)
  
  # determine sequence conservation
  sequence_conservation = df_projections.apply(lambda x: 'DC' if x['score_direct'] == 1 else 'IC' if x['score_multi'] >= thresh else 'NC', axis=1).values
  
  # create PyRanges object of projections
  df_multi = pd.DataFrame({'Chromosome': [x.chrom for x in df_projections['coords_multi']],
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
    names_FC = target_regions.overlap(pr_multi).Name.values
    functional_conservation = ['+' if x in names_FC else '-' for x in df_projections.index]
  
  # combine information of sequence and functional conservation
  conservation = [s+f for s,f in zip(sequence_conservation,functional_conservation)]
  return conservation
