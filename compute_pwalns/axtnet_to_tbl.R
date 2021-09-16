#! /usr/bin/env Rscript-4

# This script converts axtNet files to tbl format (data table)
# axtNet from CNEr show coordinates relative to the chromosome END if region is on '-' strand (only for query, reference is always '+')
# This will be changed by subtracting the values from the respective chromosome sizes.
# The result for such regions in the tbl is that the start coordinate is higher than the end coordinate.

args <- commandArgs(trailingOnly=T)
if (length(args) != 3) stop('Usage: ./axtnet_to_tbl.R S1.S2.net.axt.gz S2.chromsizes S1.S2.tbl')

packages <- c('rtracklayer', 'CNEr')
for (pkg in packages) suppressMessages(library(pkg, character.only=T))

# read axtNet
axt <- readAxt(args[1])
df <- cbind(data.frame(axt@first)[,c('seqnames','start','end')], data.frame(axt@second)[,c('seqnames','start','end')])
colnames(df) <- c('ref_chrom', 'ref_start', 'ref_end', 'qry_chrom', 'qry_start', 'qry_end')

# invert coordinates of query on '-' strand to relate to the chromosome start, not the end.
genome_size_qry_df <- read.table(args[2])
genome_size_qry <- genome_size_qry_df$V2
names(genome_size_qry) <- genome_size_qry_df$V1
chrsize <- genome_size_qry[as.vector(seqnames(axt@second[strand(axt@second)=='-']))]
df$qry_start[which(strand(axt@second)=='-')] <- chrsize - start(axt@second[strand(axt@second)=='-'])
df$qry_end[which(strand(axt@second)=='-')] <- chrsize - end(axt@second[strand(axt@second)=='-'])

# write table
write.table(df, args[3], sep='\t', col.names=F, row.names=F, quote=F)
