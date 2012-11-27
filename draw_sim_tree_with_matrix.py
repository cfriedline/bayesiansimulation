__author__ = 'chris'

import os
import sys
from ete2 import Tree, TreeNode, TreeStyle, AttrFace, TextFace, ClusterTree, ClusterNode, ProfileFace, CircleFace, NodeStyle, PhyloTree
import numpy
import tempfile

def get_color(data, row, col):
    val = float(data[row, col])
    coldata = get_column(data, col)
    colsum = numpy.sum(coldata)
    colmed = numpy.median(coldata)
    colavg = numpy.average(coldata)
    colsd = numpy.std(coldata)
    colmin = numpy.min(coldata)
    colptp = numpy.ptp(coldata)
    n = (val-colmin)/colptp

    if n >= 0.5:
        first = int(round(n*255))
        first = hex(first)[2:]
        third = "00"
    else:
        third = int(round(n*255))
        third = hex(third)[2:]
        first = "00"
    color = "#%s%s%s" % (first, "00", third)
    return color


def get_tree_style(tree_file, abund, rownames):

    with open("matrix.txt", "w") as temp:
        cols = len(abund[0])
        header = "#Names"
        for i in xrange(cols):
            header += "\tOTU%d" % i
        temp.write("%s\n" % header)
        for i, row in enumerate(abund):
            temp.write("%s\t%s\n" % (rownames[i], '\t'.join([str(i) for i in row])))

    t = Tree(tree_file)
    t.convert_to_ultrametric(10)

    assert isinstance(abund, numpy.ndarray)
    assert isinstance(rownames, numpy.ndarray)
    ts = TreeStyle()
    ts.mode = "r"
    ts.show_leaf_name = False
    ts.show_scale = False
    ts.show_branch_length = False
    ts.branch_vertical_margin = 20
    ts.force_topology = True
    ts.optimal_scale_level = "full"
    ts.scale = 50
    ts.draw_guiding_lines = True
    ts.guiding_lines_type = 0
    ts.guiding_lines_color = "black"
    for n in t.traverse():
        if not n.is_leaf():
            nstyle = NodeStyle()
            n.set_style(nstyle)
            nstyle['size'] = 0
            nstyle['hz_line_width'] = 3
            nstyle['vt_line_width'] = 3
        else:
            nstyle = NodeStyle()
            n.set_style(nstyle)
            nstyle['size'] = 0
            nstyle['hz_line_width'] = 3
            nstyle['vt_line_width'] = 3
            nstyle['fgcolor'] = "Black"
            nstyle['shape'] = "square"
            name_face = AttrFace("name", fsize=14, ftype="Arial", fgcolor="black", penwidth=10, text_prefix=" ", text_suffix=" ")
            n.add_face(name_face, column=0, position="aligned")
            row_index = rownames.tolist().index(n.name)
            col = 1
            for i in xrange(10):
                col += 1
                n.add_face(CircleFace(5, color=get_color(abund, row_index, i)), column=col, position="aligned")
    return t, ts

def get_tree(tree_file, abund, rownames):
    t, ts = get_tree_style(tree_file, abund, rownames)
    return t, ts


def read_data_file(file):
    data = []
    rownames = []
    f = open(file)
    for line in f:
        d = line.rstrip().split("\t")
        rownames.append(d[0])
        data.append([int(i) for i in d[1:]])
    return numpy.array(data), numpy.array(rownames)

def get_column(matrix, i):
    return numpy.array([int(row[i]) for row in matrix])


def main():
    dir = "/Users/chris/projects/bsim4/test/"
    tree_file = os.path.join(dir, "tree_8_1000_0.txt")
    abund_file = os.path.join(dir, "abund_1000_0_0.txt")
    gap_file = os.path.join(dir, "gap_1000_0_0.txt")

    abund, abund_names = read_data_file(abund_file)
    gap, gap_names = read_data_file(gap_file)

    tree, tree_style = get_tree(tree_file, abund, abund_names)
#    tree.show(tree_style=tree_style)
    tree.render("sim_tree.pdf", tree_style=tree_style)

if __name__ == '__main__':
    main()
