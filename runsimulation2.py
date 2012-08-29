__author__ = 'chris'
import numpy
import app
import rpy2.robjects as robjects
import os
import datetime
from rpy2.robjects import numpy2ri
from app import clockit
import dendropy
import sys

numpy2ri.activate()
sys.setrecursionlimit(1000000)

hostname = os.uname()[1]
project_dir = "/Users/chris/projects/bsim_maria"
n_gen = 100000
mpi = "/opt/local/bin/mpirun"
mb_old = "/Users/chris/src/mrbayes-3.1.2/mb"
mb_new = "/Users/chris/src/mrbayes_3.2.1/src/mb"
mb = mb_old
procs = 4
num_taxa = 4
num_states = 8
runs = 10
cols = [100, 1000, 9000]

if 'godel' in hostname:
    mpi = '/test/riveralab/cfriedline/bin/mpirun'
    mb = '/test/riveralab/cfriedline/src/mrbayes_3.2.1/src/mb'
    procs = 8
    project_dir = '/home/cfriedline/projects/bsim_maria'
elif 'phylogeny' in hostname:
    mpi = '/usr/local/bin/mpirun'
    mb = '/home/cfriedline/src/mrbayes_3.2.1/src/mb'
    procs = 8
    project_dir = '/home/cfriedline/projects/bsim_maria'
    n_gen = 1000000

run_dir_name = datetime.datetime.now().strftime("%m%d%y_%H%M%S")
run_dir_name = "test2"
result_dir = app.create_dir(os.path.join(project_dir, "results"))
run_dir = app.create_dir(os.path.join(result_dir, run_dir_name))
out_dir = app.create_dir(os.path.join(run_dir, "out"))
log_dir = app.create_dir(os.path.join(run_dir, "log"))
log_file = open(os.path.join(log_dir, "log.txt"), "w")
app.log_file = log_file

def create_R():
    r = robjects.r
    r('library(phangorn)')
    r('library(vegan)')

    r("""
        get_valid_matrix = function(cols, numstates) {
            found = 0
            while (found < cols) {
                root = sample(numstates)[1]
                temp = rTraitDisc(tree, model="ER", k=numstates, states=1:numstates, rate=0.1, root.value=root)
                temp = as.matrix(temp)
                if (min(temp) < max(temp)) {
                    if (found == 0) {
                        m = temp
                        roots = c(root)
                    } else {
                        m = cbind(m, temp)
                        roots = c(roots, root)
                    }
                    found = found + 1
                }
            }
        return(list(m, roots))
        }
    """)

    return r


def get_column_ranges(data):
    ranges = []
    for col in data.T:
        ranges.append((min(col), max(col)))
    return ranges


def get_bad_cols(ranges):
    bad = []
    for i, range in enumerate(ranges):
        if range[0] == range[1]:
            bad.append(i)
    return bad


def print_matrix(matrix, prefix, sample_names, num_cols, run, numpy, roots):
    fh = open(os.path.join(log_dir, "%s_%d_%d.txt" % (prefix, num_cols, run)), "w")
    m = matrix
    if not numpy:
        m = numpy.asarray(matrix)
    with fh:
        if roots:
            fh.write("#\t%s\n" % '\t'.join([str(int(elem)) for elem in roots]))
        for i, row in enumerate(m):
            fh.write("%s\t%s\n" % (sample_names[i], '\t'.join([str(int(elem)) for elem in row])))


def print_ranges(ranges, num_cols, run):
    fh = open(os.path.join(log_dir, "ranges_%d_%d.txt" % (num_cols, run)), "w")
    with fh:
        for i, range in enumerate(ranges):
            fh.write("%d\t%s\n" % (i, '\t'.join([str(int(elem)) for elem in range])))


def print_matrices(data, gap, abund, ranges, num_cols, run, sample_names, roots):
    print_matrix(data, "data", sample_names, num_cols, run, True, roots)
    print_matrix(gap, "gap", sample_names, num_cols, run, True, None)
    print_matrix(abund, "abund", sample_names, num_cols, run, True, None)
    print_ranges(ranges, num_cols, run)


@clockit
def print_state_distribution(data, num_cols, run, sample_names, dist_file):
    counts = [None] * len(data)
    for i, row in enumerate(data):
        rowdata = {}
        for elem in row:
            if elem in rowdata:
                rowdata[elem] += 1
            else:
                rowdata[elem] = 1
        counts[i] = rowdata

    keys = set()
    for i, d in enumerate(counts):
        for key, count in d.items():
            keys.add(key)

    keys = list(keys)
    keys.sort()

    dist_file.write("%d, %d\n" % (num_cols, run))
    for i, d in enumerate(counts):
        s = sample_names[i]
        for key in keys:
            count = 0.0
            if key in counts[i]:
                count = counts[i][key]
            s += "\t" + str(int(key)) + ":" + str(int(count))
        print s
        dist_file.write("%s\n" % s)
    dist_file.write("\n")
    dist_file.flush()


def print_sample_trees(r, tree, num_taxa, num_cols, run):
    assert isinstance(tree, dendropy.Tree)
    pdf_file = os.path.join(log_dir, "tree_%d_%d_%d.pdf" % (num_taxa, num_cols, run))
    text_file = os.path.join(log_dir, "tree_%d_%d_%d.txt" % (num_taxa, num_cols, run))
    r("pdf(file='%s')" % pdf_file)
    r('par(mfrow=c(1,2))')
    r("plot(tree, root.edge=T)")
    r("plot(tree, 'un')")
    r('dev.off()')

    fh = open(text_file, "w")
    with fh:
        fh.write("%s;\n" % tree.as_newick_string())


def curate_data_matrix(r, data, ranges):
    badcols = get_bad_cols(ranges)
    drop_string = ','.join([str(elem + 1) for elem in badcols])
    print drop_string
    if len(drop_string) > 0:
        r('data = data[,-c(%s)]' % drop_string)
        data = numpy.asarray(r['data'])
        ranges = get_column_ranges(data)
    return data, ranges


def get_tab_string(tuple):
    return '\t'.join([str(elem) for elem in tuple])


@clockit
def store_valid_matrix_data(r, taxa_tree, num_cols, num_states):
    print "Storing valid matrix data in R session"
    r('matrix_data = get_valid_matrix(%d, %d)' % (num_cols, num_states))
    r('data = matrix_data[[1]]')
    r('roots = matrix_data[[2]]')
    r("data = t(apply(data, 1, as.numeric))")
    robjects.globalenv['colnames'] = sorted(taxa_tree.taxon_set.labels())
    r('colnames(data) = colnames')


def get_unifrac_pcoa(tree, matrix, rownames):
    uni_tree = app.get_unifrac_pcoa_tree(matrix, rownames)
    uni_diff = app.calculate_differences_r(tree, uni_tree)
    return uni_tree, uni_diff


def get_unifrac_cluster(tree, matrix, rownames):
    uni_tree = app.get_py_unifrac_cluster(matrix, rownames)
    uni_diff = app.calculate_differences_r(tree, uni_tree)
    return uni_tree, uni_diff


def get_unifrac_nj(tree, matrix, rownames):
    uni_tree = app.get_unifrac_nj(matrix, rownames)
    uni_diff = app.calculate_differences_r(tree, uni_tree)
    return uni_tree, uni_diff


@clockit
def run_simulation(r, taxa_tree, num_cols, run, out_file, dist_file):
    assert isinstance(taxa_tree, dendropy.Tree)
    r('tree = rtree(%d)' % num_taxa)
    tree = app.ape_to_dendropy(r['tree'])
    store_valid_matrix_data(r, taxa_tree, num_cols, num_states)
    data = numpy.asarray(r['data'])
    roots = list(r['roots'])
    sample_names = numpy.asarray(r('rownames(data)'))
    print_sample_trees(r, tree, num_taxa, num_cols, run)
    print_state_distribution(data, num_cols, run, sample_names, dist_file)
    ranges = get_column_ranges(data)
    data, ranges = curate_data_matrix(r, data, ranges)
    gap = app.restandardize_matrix(data, ranges, num_states = 4)
    abund_ranges = app.get_range_from_gamma(num_cols * 2, 2, 1, 1000, app.compute_smallest_max())
    abund = app.get_abundance_matrix(gap, abund_ranges, "gamma", num_states = 4)
    (u_uni_matrix, u_uni_rownames), (w_uni_matrix, w_uni_rownames) = app.calculate_unifrac(abund, sample_names,
                                                                                           taxa_tree)

    u_unifrac_pcoa_tree, u_unifrac_pcoa_diffs = get_unifrac_pcoa(tree, u_uni_matrix, u_uni_rownames)
    w_unifrac_pcoa_tree, w_unifrac_pcoa_diffs = get_unifrac_pcoa(tree, w_uni_matrix, w_uni_rownames)

    u_unifrac_cluster_tree, u_unifrac_cluster_diffs = get_unifrac_cluster(tree, u_uni_matrix, u_uni_rownames)
    w_unifrac_cluster_tree, w_unifrac_cluster_diffs = get_unifrac_cluster(tree, w_uni_matrix, w_uni_rownames)

    u_unifrac_nj_tree, u_unifrac_nj_diffs = get_unifrac_nj(tree, u_uni_matrix, u_uni_rownames)
    w_unifrac_nj_ree, w_unifrac_nj_diffs = get_unifrac_nj(tree, w_uni_matrix, w_uni_rownames)

    disc = app.get_discrete_matrix_from_standardized(gap, 2, sample_names)
    mb_tree = app.run_mrbayes(run, disc, sample_names, disc.ncol, n_gen, mpi, mb, procs, None, run_dir,
                              len(sample_names),
                              "d")
    mb_diffs = app.calculate_differences_r(tree, mb_tree)

    out_file.write("%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (run, num_cols,
                                                             get_tab_string(mb_diffs),
                                                             get_tab_string(u_unifrac_pcoa_diffs),
                                                             get_tab_string(u_unifrac_cluster_diffs),
                                                             get_tab_string(u_unifrac_nj_diffs),
                                                             get_tab_string(w_unifrac_pcoa_diffs),
                                                             get_tab_string(w_unifrac_cluster_diffs),
                                                             get_tab_string(w_unifrac_cluster_diffs)))
    out_file.flush()
    print_matrices(data, gap, abund, abund_ranges, num_cols, run, sample_names, roots)


def print_taxa_tree(tree, num_cols):
    assert isinstance(tree, dendropy.Tree)
    r = robjects.r
    r("temp = read.tree(text='%s;')" % tree.as_newick_string())
    pdf_file = os.path.join(log_dir, "taxa_tree_%d.pdf" % num_cols)
    text_file = os.path.join(log_dir, "taxa_tree_%d.txt" % num_cols)
    r("pdf(file='%s')" % pdf_file)
    r("plot(temp, root.edge=T)")
    r('dev.off()')

    fh = open(text_file, "w")
    with fh:
        fh.write("%s;\n" % tree.as_newick_string())


def get_header():
    return "run\tcols\t"\
           "mb_topo\tmb_sym\tmb_path\t"\
           "u_pcoa_topo\tu_pcoa_symm\tu_pcoa_path\t"\
           "u_cluster_topo\tu_cluster_symm\tu_cluster_path\t"\
           "u_nj_topo\tu_nj_symm\tu_nj_path\t"\
           "w_pcoa_topo\tw_pcoa_symm\tw_pcoa_path\t"\
           "w_cluster_topo\tw_cluster_symm\tw_cluster_path\t"\
           "w_nj_topo\tw_nj_symm\tw_nj_path"


@clockit
def main():
    r = create_R()
    out_file = open(os.path.join(out_dir, "out.txt"), "w")
    out_file.write("%s\n" % get_header())
    dist_file = open(os.path.join(out_dir, "dists.txt"), "w")
    for col in cols:
        tree = app.create_tree(col, "T")
        print_taxa_tree(tree, col)
        for run in xrange(runs):
            run_simulation(r, tree, col, run, out_file, dist_file)
    out_file.close()
    dist_file.close()


if __name__ == '__main__':
    main()
    print "Done!"


