#! /usr/bin/env python

# This script is used to convert chain files into "chain tables".
# Every row in the chain table will contain an ungapped alignment block between the ref and the qry species.
# The columns will refer to chromosome, start and end coordinates of each alignment block in the ref and the qry species.

def main():
    # import packages and parse arguments
    import os, sys
    if len(sys.argv) != 3:
        print('Usage: python chainToTbl.py chain-file outfile')
        sys.exit(1)
    _, chain, outfile = sys.argv
    
    # read chain file
    with open (chain, 'r') as f:
        lines = f.readlines()

    # write coordinates of alignment blocks from chain file to outfile
    with open(outfile, 'w') as f:
        for line in lines:
            if (line.startswith('#')) or (line == '\n'):
                continue
            if line.startswith('chain'):
                line = line.strip().split()
                chrom_ref, chrom_qry, strand_qry = [line[i] for i in (2,7,9)]
                start_ref, end_ref, chrom_size_qry, start_qry, end_qry = [int(line[i]) for i in (5,6,8,10,11)]
            else:
                line = [int(x) for x in line.strip().split()] # [block_width, gap_after_block_ref, gap_after_block_qry]
                aln_block = [chrom_ref,
                             start_ref,
                             start_ref + line[0] - 1,
                             chrom_qry, 
                             start_qry,
                             start_qry + line[0] - 1]
                # the chain file stores coordinates on the reverse strand from the right end, i.e. chrom_size - x.
                if strand_qry == '-':
                    # the end coordinate is inclusive. an alignment block of
                    # size 3 at the end of a chain [0,5) should have the 
                    # coordinates [2,4] (i.e. the positions 2,3,4).
                    aln_block[4] = chrom_size_qry - start_qry - 1
                    aln_block[5] = chrom_size_qry - start_qry - line[0]
                f.write('\t'.join(map(str, aln_block)) + '\n')
                if len(line) == 3:
                    start_ref += sum(line[:2]) # add block width plus gap size to current position
                    start_qry += sum([line[0],line[2]])

    # sort outfile
    os.system('sort -k1,1 -k2,2n %s > %s.tmp; mv %s.tmp %s' %((outfile,)*4))

if __name__ == '__main__':
    main()
