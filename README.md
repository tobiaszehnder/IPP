## IPP - Independent Point Projection

A tool for comparative genomics beyond sequence alignments

> Tobias Zehnder


## Method Summary

For a genomic region with conserved synteny, any non-alignable coordinate can be approximately projected from one genome to another by linear interpolation of its relative position between two alignable anchor points.
The accuracy of such interpolations negatively correlates with the distance to the anchor points. Therefore, projections between species with large evolutionary distances (e.g. > 200 my) tend to be inaccurate due to a low anchor point density. Including so-called bridging species may increase the anchor point density and thus improve projection accuracy.
The optimal choice and combination of bridging species may vary from one genomic location to another. This presents a shortest path problem in a graph where every node is a species and the weighted edges between nodes represent distances of genomic locations to their anchor points (|x - a|). We established a scoring function that exponentially decreases with increasing distances |x - a|. The shortest path problem is solved using Dijkstra’s Shortest Path Algorithm (Dijkstra, 1959).


## Installation

1. Download or clone the repository from github and go to the diretory.
2. Compile the module: `python setup.py build`
3. Tell python where to look for your module.
   Use the directory that is created according to your python version, e.g.:
   `export PYTHONPATH=/path/to/your/IPP_directory/build/lib.linux-x86_64-3.10/`
   Check the last folder, it might be different than what is stated here.
   Write the line to your ~/.bashrc or ~/.bash_profile if you want it to be set in every new shell session.
4. Make sure the following python modules are installed:
   `os, sys, numpy, pandas, argparse, tabulate, tqdm, pyranges`

## Alignment Data

IPP uses a collection of pairwise alignments (pwaln) between the reference, target, and all bridging species stored in binary format.
We provide a set of precomputed pwaln collection files for selected species comparisons.
These are large files stored separately from github and can be downloaded <a href="https://oc-molgen.gnz.mpg.de/owncloud/s/ACWNtKRCiN8BYxi">HERE</a>.

Alternatively, we provide a pipeline to compute your own alignment collections for your choice of species.
For that, run `compute_alignments/compute_pairwise_alignments`
A Snakefile guides through the whole alignment process from fasta to chain file.
With that, you can create the final collection of pairwise alignments between reference, target and a given set of bridging species that is used for running IPP.


## Usage

Project the regions of a bed-file from one genome to another:
```bash
python project.py
```
Type `python project.py -h` for a detailed description.

Questions, suggestions and criticism are welcome at zehnder[at]molgen.mpg.de
