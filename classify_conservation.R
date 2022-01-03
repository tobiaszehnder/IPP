#! /usr/bin/env Rscript

## This script will classify reference regions from a projection file (.proj) into directly (DC), indirectly (IC) and not conserved (NC) elements based on a given threshold on the projection scores:
## DC: direct score >= thresh
## IC: direct score < thresh && dijkstra score >= thresh
## NC: all the rest
## In addition, the script takes an optional argument specifying the target regions to identify regions of functional conservation (+/-)
## Color coding: DC+ and IC+: red, DC- and IC-: black, NC: grey.
## If target regions are not given, this script simply classifies into DC (green), IC (orange) and NC (grey).

## parse arguments
args <- commandArgs(trailingOnly=T)
if (length(args) < 6) stop('Usage: classify_conservation.R reference outfile reference_bedfile projection_bedfile direct_chain_file threshold <optional: target_bedfile>')
reference <- args[1]
outfile <- args[2]
reference_bedfile <- paste0(args[3], '_centered')
if (!file.exists(reference_bedfile)) stop(sprintf('reference_bedfile not found: %s', reference_bedfile))
projection_bedfile <- args[4]
if (!file.exists(projection_bedfile)) stop(sprintf('projection_bedfile not found: %s', projection_bedfile))
direct_chain_file <- args[5]
if (!file.exists(direct_chain_file)) stop(sprintf('direct_chain_file not found: %s', direct_chain_file))
thresh <- as.numeric(args[6])
target_bedfile <- NA
if (length(args) == 7) {
	target_bedfile <- args[7]
	if (!file.exists(target_bedfile)) stop("target_bedfile doesn't exist.")
}

## packages
packages <- c('rtracklayer')
for (pkg in packages) suppressMessages(library(pkg, character.only=T))

## read reference regions and projections and classify into DC/IC/NC
gr <- import.bed(reference_bedfile)
chain <- import.chain(direct_chain_file)
chain_gr <- suppressWarnings(do.call('c', unname(sapply(names(chain), function(chrom) GRanges(seqnames=chrom, ranges=chain[[chrom]]@ranges)))))
projections <- import.bed(projection_bedfile)
gr$score <- 0
gr$score[gr$name %in% projections$name] <- projections$score
gr$class <- 'NC'
gr$class[gr$score >= thresh] <- 'IC'
gr$class[overlapsAny(gr, chain_gr)] <- 'DC'
gr$target_coord <- NA
gr$target_coord[gr$name %in% projections$name] <- paste(seqnames(projections), start(projections), sep=':')
gr$name <- paste(gr$name, gr$target_coord, sep='_')

## read target bedfile (if given) to classify into +/- and write to file
maxgap <- 500 # maximum gap for overlaps to still be detected
if (!is.na(target_bedfile)) {
    target_bed <- import.bed(target_bedfile)
    gr_projections <- GRanges(rep('chr0:0', length(gr))) # create 'dummy' regions for 'NA' target projectons such that the indices of findOverlaps can be applied back to the GRanges.
    suppressWarnings(gr_projections[!is.na(gr$target_coord)] <- GRanges(na.omit(gr$target_coord)))
    ov <- suppressWarnings(findOverlaps(gr_projections, target_bed, maxgap=maxgap))
    gr$class <- paste0(gr$class, '-')
    gr$class[queryHits(ov)] <- gsub('-', '+', gr$class[queryHits(ov)])
    colors <- c(`DC+`='255,0,0', `IC+`='255,0,0', `DC-`='30,30,30', `IC-`='30,30,30', `NC+`='141,153,174', `NC-`='141,153,174')
} else {
    colors <- c(DC='127,201,127', IC='253,180,98', NC='141,153,174')
}
gr$itemRgb <- sapply(gr$class, function(cl) colors[cl])
gr$name <- paste(gr$name, gr$class, sep='_')
df <- data.frame(gr)[,c('seqnames','start','end','name','score','strand','start','end','itemRgb')]
colnames(df) <- c('seqnames','start','end','name','score','strand','thickStart','thickEnd','itemRgb')
df$start <- df$thickStart <- df$start - 1
df$strand <- '.'
df[,c('start','end','thickStart','thickEnd')] <- format(df[,c('start','end','thickStart','thickEnd')], scientific=F)
write.table(df, outfile, sep='\t', quote=F, col.names=F, row.names=F)
