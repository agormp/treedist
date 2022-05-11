#! /usr/bin/env python3
# Computes the normalized symmetric similarity between trees listed in two or more files
# (only first tree in each file is considered)
# Can also print other measures, and details concerning how many and which bipartitions that differ between trees

import sys, os.path, itertools
import phylotreelib as treelib
from optparse import OptionParser

# Build commandline parser
parser = OptionParser(usage="usage: %prog [options] TREEFILE1 TREEFILE2 [TREEFILE3 ...]",
                      version="1.1")

parser.add_option("-I", "--informat", type="choice", dest="format",
                  choices=["NEXUS", "nexus", "NEWICK", "newick"], metavar="FORM",
                  help="format of tree files: NEXUS or newick [default: newick]")

parser.add_option("-d", "--symnormdist", action="store_true", dest="symnormdist",
                      help="print normalized symmetric distance")

parser.add_option("-s", "--symdist", action="store_true", dest="symdist",
                      help="print symmetric distance")

parser.add_option("-b", "--bipnumbers", action="store_true", dest="numbers",
                      help="print no. total, shared, and unique bipartitions for each tree")

parser.add_option("-u", "--uniquelist", action="store_true", dest="uniqlist",
                      help="print list of unique bipartitions for each tree")

parser.add_option("-c", "--commonlist", action="store_true", dest="commonlist",
                      help="print list of bipartitions trees have in common")

parser.add_option("-l", "--labels", action="store_true", dest="labels",
                      help="also print branch labels for bipartitions (requires -c or -u)")

parser.add_option("-o", "--oneline", action="store_true", dest="oneline",
                      help="oneline output of normalized symmetric similarities")

parser.add_option("-n", "--namesummary", action="store_true", dest="namesummary",
                      help="Print names of leaves that are shared between trees, and leaves that are discarded")

parser.add_option("--collapselen", type="string", dest="collapselen", metavar="BRLEN[,BRLEN,...]",
                        help="Collapse branches shorter than this (comma-separated list: individual cutoffs for trees)")

parser.add_option("--collapsefrac", type="string", dest="collapsefrac", metavar="FRAC",
                        help="Collapse branches shorter than FRAC fraction of treeheight (of minimum variance rooted tree)")

parser.add_option("--collapsesupport", type="string", dest="collapsesupport", metavar="SUPPORT",
                        help="Collapse branches with less clade support than this (must be number in [0,1])")

parser.set_defaults(format="newick", symnormdist=False, symdist=False,
                    numbers=False, uniqlist=False, commonlist=False, oneline=False, namesummary=False,
                    collapselen=None, collapsefrac=None, collapsesupport=None)

# Parse commandline, open treefiles and get trees (if file names are given and files exist)
(options, args) = parser.parse_args()
if len(args) < 2:
    parser.error("please provide the name of two or more tree files")
else:
    tree_list = []
    for filename in args:
        if not os.path.isfile(filename):
            parser.error("file %s not found." % filename)
        else:
            if options.format.lower() == "nexus":
                treefile = treelib.Nexustreefile(filename)
            else:
                treefile = treelib.Newicktreefile(filename)

            tree = next(treefile)
            tree.name = os.path.basename(filename)        # add .name attribute to tree object for later identification
            tree_list.append(tree)

    # Parse list of branch length cutoffs for collapsing, if present
    if options.collapselen:
        brlencutoffs = []
        collapse_words = options.collapselen.split(",")

        if len(collapse_words) == 1:
            brlencutoffs = len(tree_list) * [float(collapse_words[0])]
        else:
            for word in collapse_words:
                brlencutoffs.append(float(word))

#Prune larger tree(s) so leaf sets become identical

# Find set of leaves that are shared between all trees
shared_leafs = set(tree_list[0].leaves)
for tree in tree_list[1:]:
    shared_leafs.intersection_update(tree.leaves)
sharedlist = list(shared_leafs)
sharedlist.sort()

# If trees share less than four leafs: bail out
if len(shared_leafs) < 4:
    print("Error: trees have less than 4 leafs in common - comparison is not meaningful")
    print("Shared leafs (if any):")
    for leaf in shared_leafs:
        print("    %s" % leaf)
    sys.exit()

# Prune trees so they only contain common set of leafs:
for tree in tree_list:
    if tree.leaves != shared_leafs:
        remove_set = tree.leaves - shared_leafs
        for leaf in remove_set:
            tree.remove_leaf(leaf)
        removelist = list(remove_set)   # This is just from last tree in treelist??? Should perhaps create cumulated list of discarded leaves
        removelist.sort()

# If requested: collapse branches shorter than collapsefrac fraction of treeheight (after minvar rooting)
# First computes branch length cutoff from frac cutoff for each tree, then uses collapselen to do work
if options.collapsefrac:
    options.collapsefrac = float(options.collapsefrac)
    if options.collapsefrac < 0 or options.collapsefrac > 1:
        raise treelib.TreeError("--collapsefrac option requires value between 0.0 and 1.0 (actual value: {})".format(options.collapsefrac))
    options.collapselen = True
    brlencutoffs = []
    for i in range(len(tree_list)):
        tree = tree_list[i]
        tree.rootminvar()
        cutoff = options.collapsefrac * tree.height
        brlencutoffs.append(cutoff)

# If requested: collapse short branches
if options.collapselen:
    for i in range(len(tree_list)):
        tree = tree_list[i]
        cutoff = brlencutoffs[i]
        unchecked_parents = tree.intnodes.copy()
        while(unchecked_parents):
            parent = unchecked_parents.pop()
            unchecked_children = tree.children(parent) & tree.intnodes  # Set intersection: do not check leaves
            while unchecked_children:
                child = unchecked_children.pop()
                blen = tree.tree[parent][child].length
                if blen <= cutoff:
                    grand_children = tree.children(child) & tree.intnodes
                    unchecked_children.update(grand_children)    # Add child's children to parent
                    unchecked_parents.remove(child)
                    tree.remove_branch(parent, child)

# If requested: collapse branches with low branch support
if options.collapsesupport:
    suportcutoff = float(options.collapsesupport)
    for i in range(len(tree_list)):
        tree = tree_list[i]
        unchecked_parents = tree.intnodes.copy()
        while(unchecked_parents):
            parent = unchecked_parents.pop()
            unchecked_children = tree.children(parent) & tree.intnodes  # Set intersection: do not check leaves
            while unchecked_children:
                child = unchecked_children.pop()
                label = tree.tree[parent][child].label
                if not label:
                    raise treelib.TreeError("The tree does not have branch support values: {}".format(tree.name))
                else:
                    try:
                        support = float(label)
                    except ValueError:
                        raise TreeError("This branch label is not a floating point number: {} {}".format(label, tree.name))
                    if support < suportcutoff:
                        grand_children = tree.children(child) & tree.intnodes
                        unchecked_children.update(grand_children)    # Add child's children to parent
                        unchecked_parents.remove(child)
                        tree.remove_branch(parent, child)


# Find set of bipartitions in each tree
# Recall that: Names of leafs on one side of a branch are represented as an immutable set.
# A bipartition is represented as an immutable set of two such (complementary) sets
# The entire tree topology is represented as a set of bipartitions
bipart_list = {}
for tree in tree_list:
    bipart_list[tree] = tree.topology()

# Print header if oneline output requested
# if options.oneline:
#     sys.stdout.write("Normalized symmetric similarity: ")

# Compare all trees against all others
for tree1, tree2 in itertools.combinations(tree_list, 2):

    # Find bipartitions for tree1
    tree1_biparts = bipart_list[tree1]
    if options.labels:
        tree1_bipdict = tree1.bipdict()

    # Find bipartitions for tree2
    tree2_biparts = bipart_list[tree2]
    if options.labels:
        tree2_bipdict = tree2.bipdict()

    # Find bipartitions unique to tree1 and tree2 (using set arithemtic)
    tree1_unique_biparts = tree1_biparts - tree2_biparts
    tree2_unique_biparts = tree2_biparts - tree1_biparts

    # Find shared bipartitions
    shared_biparts = tree1_biparts & tree2_biparts

    # Compute results
    num_shared = len(shared_biparts) - len(tree1.leaves)  # Only internal branches counts!!!
    num_uniq1 = len(tree1_unique_biparts)
    num_uniq2 = len(tree2_unique_biparts)
    num_bip1 = len(tree1_biparts) - len(tree1.leaves)  # Only internal branches counts!!!
    num_bip2 = len(tree2_biparts) - len(tree2.leaves)  # Only internal branches counts!!!
    symdif = num_uniq1 + num_uniq2
    symdifnorm = 1.0*symdif/(num_bip1 + num_bip2)
    symsimnorm = 1.0 - symdifnorm


    # Prints pretty list of bipartitions, one bipartition per line, "#" used as a delimiter between parts
    # NOTE: bipartitions corresponding to external branches (leading to leafs) are not printed
    def print_biparts(biplist, printlab1=False, printlab2=False):
        for bipartition in biplist:
            setlist = list(bipartition)
            names1 = list(setlist[0])
            names1.sort()
            names2 = list(setlist[1])
            names2.sort()

            # Only print bipartition if this is not an external branch (<=> one bipart has one member)
            if (len(names1) !=1) and (len(names2) !=1):
                if len(names1) < len(names2):
                    first = names1
                    second = names2
                else:
                    first = names2
                    second = names1

                for name in first:
                    sys.stdout.write(name + " ")
                sys.stdout.write("  #   ")
                for name in second:
                    sys.stdout.write(name + " ")
                sys.stdout.write("\n")

                # If branch labels are requested, print those too
                if printlab1 and printlab2:
                    lab1 = tree1_bipdict[bipartition].label
                    lab2 = tree2_bipdict[bipartition].label
                    print("\n    Branch labels for bipartition:  %s: %s    %s: %s\n" % (tree1.name, lab1, tree2.name, lab2))
                elif printlab1 and not printlab2:
                    lab1 = tree1_bipdict[bipartition].label
                    print("\n    Branch label for bipartition:  %s: %s\n" % (tree1.name, lab1))
                elif not printlab1 and printlab2:
                    lab2 = tree2_bipdict[bipartition].label
                    print("\n    Branch label for bipartition:  %s: %s\n" % (tree2.name,lab2))


    # Print results
    if options.oneline:
        sys.stdout.write("%s\t%s\t%5.3f\t" % (tree1.name, tree2.name, symsimnorm))
    else:
        print("\n  Comparison of %s and %s:" % (tree1.name, tree2.name))
        print("    Normalized symmetric similarity: %4.2f" % (symsimnorm))
        if options.symnormdist:
            print("    Normalized symmetric distance: %6.2f" % (symdifnorm))
        if options.symdist:
            print("    Symmetric distance: %17d" % (symdif))
        if options.numbers:
            print("\n    Total number of bipartitions:  %s: %3d    %s: %3d" % (tree1.name, num_bip1, tree2.name, num_bip2))
            print("    Number of unique bipartitions: %s: %3d    %s: %3d" % (tree1.name, num_uniq1, tree2.name, num_uniq2))
            print("    Number of shared bipartitions: %d" % (num_shared))

        if options.uniqlist:
            print("\n    List of unique bipartitions in %s:" % (tree1.name))
            if num_uniq1 == 0:
                print("        no unique biparts\n")
            else:
                print_biparts(tree1_unique_biparts, printlab1=options.labels)

            print("\n    List of unique bipartitions in %s" % (tree2.name))
            if num_uniq2 == 0:
                print("        no unique biparts\n")
            else:
                print_biparts(tree2_unique_biparts, printlab2=options.labels)
            print("")
        if options.commonlist:
            print("\n    List of bipartitions present in both trees\n")
            print_biparts(shared_biparts, printlab1=options.labels, printlab2=options.labels)

        if options.namesummary:
            print("\n    Leaves present in both trees (N = {}):\n".format(len(sharedlist)))
            for name in sharedlist:
                print("       {}".format(name))
            print("\n    Leaves discarded from larger tree (N = {}):\n".format(len(removelist)))
            for name in removelist:
                print("       {}".format(name))
        print()

if options.oneline:
    sys.stdout.write("\n")

