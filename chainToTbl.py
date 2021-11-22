#! /usr/bin/env python

# This script is used to convert chain files into "chain tables".
# Every row in the chain table will contain an ungapped alignment block between the ref and the qry species.
# The columns will refer to chromosome, start and end coordinates of each alignment block in the ref and the qry species.

def main():
    # import packages and parse arguments
    import os, sys
    if len(sys.argv) != 3:
        print('Usage: python chainToTbl.py chain-file outdir')
        sys.exit(1)
    _, chain, outdir = sys.argv
    outdir = os.path.realpath(outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outfile = os.path.join(outdir, chain.split('/')[-1] + '.tbl')
    
    # read chain file
    with open (chain, 'r') as f:
        lines = f.readlines()

    # get alignment blocks from chain file
    with open(outfile, 'w') as f:
        for line in lines:
            if (line.startswith('#')) or (line == '\n'):
                continue
            if line.startswith('chain'):
                line = line.strip().split()
                chrom_ref, start_ref, end_ref, chrom_qry, start_qry, end_qry = [line[i] for i in (2,5,6,7,10,11)]
                x_ref, x_qry = 0, 0
            else:
                line = list(map(int, line.strip().split()))
                aln_block = [chrom_ref, int(start_ref) + x_ref, int(start_ref) + x_ref + line[0], chrom_qry,  int(start_qry) + x_qry, int(start_qry) + x_qry + line[0]]
                f.write('\t'.join(map(str, aln_block)) + '\n')
                if len(line) == 3:
                    x_ref += sum(line[:2]) # add block width plus gap size to current position
                    x_qry += sum([line[0],line[2]])

    # sort outfile
    os.system('sort -k1,1 -k2,2n %s > tmp; mv tmp %s' %(outfile,outfile))

if __name__ == '__main__':
    main()
