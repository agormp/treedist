# treedist: command line program for quantifying differences between phylogenetic trees

![](https://img.shields.io/badge/version-1.0.0-blue)

The command-line program `treedist` takes as input two or more treefiles and computes various measures of the "distance" between pairs of trees. The program can also output information about how many and what bipartitions differ between trees.


## Availability

The `treedist` source code is available on GitHub: https://github.com/agormp/treedist. The executable can be installed from PyPI: https://pypi.org/project/treedist/

## Installation

```
python3 -m pip install treedist
```

## Dependencies

`treedist` relies on the [phylotreelib library](https://github.com/agormp/phylotreelib), which is automatically included when using pip to install.

## Overview

* Input:
	* Two or more files containing phylogenetic trees in NEXUS or Newick format
	* Trees do not need to have identical leaves, but each pair of trees need to share at least 4 leaves (all leaves that are not shared will be automatically pruned before computing tree distances).
* Output:
	* Different measures of treesimilarity (Robinson-Fould's symmetric distance, normalised RF, normalised similarity)
	* Summaries of what bipartitions that differ or are shared between trees

## Usage

```
Usage: treedist [options] TREEFILE1 TREEFILE2 [TREEFILE3 ...]

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -I FORM, --informat=FORM
                        format of tree files: NEXUS or newick [default:
                        newick]
  -d, --symnormdist     print normalized symmetric distance
  -s, --symdist         print symmetric distance
  -b, --bipnumbers      print no. total, shared, and unique bipartitions for
                        each tree
  -u, --uniquelist      print list of unique bipartitions for each tree
  -c, --commonlist      print list of bipartitions trees have in common
  -l, --labels          also print branch labels for bipartitions (requires -c
                        or -u)
  -o, --oneline         oneline output of normalized symmetric similarities
  -n, --namesummary     Print names of leaves that are shared between trees,
                        and leaves that are discarded
  --collapselen=BRLEN[,BRLEN,...]
                        Collapse branches shorter than this (comma-separated
                        list: individual cutoffs for trees)
  --collapsefrac=FRAC   Collapse branches shorter than FRAC fraction of
                        treeheight (of minimum variance rooted tree)
  --collapsesupport=SUPPORT
                        Collapse branches with less clade support than this
                        (must be number in [0,1])
```

## Usage examples

TBD