__author__ = 'chris'

from Bio import Phylo
import dendropy
from cStringIO import StringIO
import Bio
import pylab
from ete2 import Tree, TreeStyle, TreeNode, NodeStyle, TextFace, CircleFace, AttrFace, faces
import ete2

def get_node_style(n, metadata, colors):
    assert isinstance(n, TreeNode)
    line_width=3

    if n.is_root():
        n.dist = 0

    if n.is_leaf():
        name = n.name.split(".")[0]
        n.name = name
        site = metadata[name]
        color = colors[site]
        nstyle = NodeStyle()
        nstyle['fgcolor'] = color
        nstyle['size'] = 10
        nstyle['shape'] = 'circle'
#        nstyle['hz_line_width'] = line_width
#        nstyle['vt_line_width'] = line_width
        n.set_style(nstyle)
        name_face = AttrFace("name", fsize=14, ftype="Arial", fgcolor=color, penwidth=10, text_prefix=" ", text_suffix=" ")
#        n.add_face(name_face, column=1, position="aligned")
    else:
        nstyle = NodeStyle()
        nstyle['size'] = 0
#        nstyle['hz_line_width'] = line_width
#        nstyle['vt_line_width'] = line_width
        nstyle['fgcolor'] = "Black"
        nstyle['shape'] = "square"
        n.set_style(nstyle)

def get_node_styles(tree, metadata, colors):
    assert isinstance(tree, TreeNode)
    for n in tree.traverse():
        n.set_style(get_node_style(n, metadata, colors))


def get_tree_style(tree, metadata, colors):
    assert isinstance(tree, TreeNode)
    ts = TreeStyle()
    ts.mode = "c"
    ts.show_leaf_name = False
    ts.show_scale = False
    ts.show_branch_length = False
    ts.force_topology = True
    ts.optimal_scale_level = "mid"
    get_node_styles(tree, metadata, colors)

    for site, color in colors.items():
        if 'Gastro' in site:
            site = "GI"
        ts.legend.add_face(CircleFace(7, color), column=0)
        ts.legend.add_face(TextFace(" %s" % site, fsize=20), column=1)

    ts.aligned_header.add_face(TextFace("fuck you"), column=0)


    return ts

def get_metadata(map_file):
    bodysitemap = {}
    m = open(map_file)
    m.readline()
    for line in m:
        data = line.rstrip().split("\t")
        bodysitemap[data[4]] = data[13]
    return bodysitemap

def get_colors(metadata):
    sites = set()
    for k, v in metadata.items():
        sites.add(v)
    sites = list(sites)
    sites.sort()
    color_map = {}
    colors = ['Medium blue', "Red", "Aqua", "ForestGreen"]
    for i, s  in enumerate(sites):
        color_map[s] = colors[i]
    return color_map


def main():
    tree_file = '/Users/chris/projects/dacc/110912_190012/out/mb/mb_73_9840_n_0_dacc.nex.con.tre'
    map_file = "ppAll_V69_map.txt"

    metadata = get_metadata(map_file)
    colors = get_colors(metadata)


    tree  = dendropy.Tree.get_from_path(tree_file, "nexus")
    assert isinstance(tree, dendropy.Tree)
    newick = tree.as_newick_string()

    t = Tree("%s;" % newick)
    ts = get_tree_style(t, metadata, colors)
#    t.convert_to_ultrametric(tree_length=10)
    t.show(tree_style=ts)
    t.render("dacc_tree.svg", tree_style=ts)


if __name__ == '__main__':
    main()
