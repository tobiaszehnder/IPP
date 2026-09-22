#!/usr/bin/env python

# This script pickles a set of given pairwise alignment tables (.tbl format)

# parse arguments
import sys


if len(sys.argv) < 5:
    print(
        "Usage: ./pickle_pwalns.py outfile.pkl missing_pwalns_allowed/forbidden sp1.sp2.tbl sp1.sp3.tbl ..."
    )
    sys.exit(0)

outfile, missing_pwalns = sys.argv[1:3]
tbls = sys.argv[3:]
if missing_pwalns == "missing_pwalns_allowed":
    missing_pwalns_allowed = True
elif missing_pwalns == "missing_pwalns_forbidden":
    missing_pwalns_allowed = False
else:
    print(
        "Usage: ./pickle_pwalns.py outfile.pkl missing_pwalns_allowed/forbidden sp1.sp2.tbl sp1.sp3.tbl ..."
    )
    sys.exit(0)

# import packages
import os
import pickle
import sys
from collections import defaultdict

import pandas as pd


# read pairwise alignment tbl files
cols = ["ref_chrom", "ref_start", "ref_end", "qry_chrom", "qry_start", "qry_end"]
# species = set([x_i for x in tbls for x_i in x.split('.')[:2]])
tbl_dict = defaultdict(dict)
for tbl in tbls:
    sp1, sp2, _suffix = tbl.split("/")[-1].split(".")
    tbl_dict[sp1][sp2] = tbl

if missing_pwalns_allowed:
    pwaln = {
        ref: {
            qry: pd.read_csv(tbl, sep="\t", header=None, names=cols)
            for qry, tbl in tbl_dict[ref].items()
            if os.path.isfile(tbl)
        }
        for ref in tbl_dict.keys()
    }
else:
    pwaln = {
        ref: {
            qry: pd.read_csv(tbl, sep="\t", header=None, names=cols)
            for qry, tbl in tbl_dict[ref].items()
        }
        for ref in tbl_dict.keys()
    }

# pickle the dict data frames and write it to file
with open(outfile, "wb") as f:
    pickle.dump(pwaln, f)
