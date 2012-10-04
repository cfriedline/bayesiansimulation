__author__ = 'chris'

# uses constant branch length on the sample trees

import numpy
import app
import rpy2.robjects as robjects
import os
import datetime
from rpy2.robjects import numpy2ri
from app import clockit
import dendropy
import sys
import multiprocessing as mp
import argparse

numpy2ri.activate()
sys.setrecursionlimit(1000000)

hostname = os.uname()[1]
project_dir = "/Users/chris/projects/bsim3"
n_gen = 100000
mpi = "/opt/local/bin/mpirun"
mb_old = "/Users/chris/src/mrbayes-3.1.2/mb"
mb_new = "/Users/chris/src/mrbayes_3.2.1/src/mb"
mb = mb_old
procs = 4
num_taxa = 8
num_states = 8
bits = 3
runs = 3
distance_runs = 3
rate = 1.0
brlen = 0.5
cols = [100, 1000, 9000]
cols = [1000]
gamma_shape = 1
gamma_scale = 1000

if 'godel' in hostname:
    mpi = '/test/riveralab/cfriedline/bin/mpirun'
    mb = '/test/riveralab/cfriedline/src/mrbayes_3.2.1/src/mb'
    procs = 8
    project_dir = '/home/cfriedline/projects/bsim3'
    mb_new = mb
elif 'phylogeny' in hostname:
    mpi = '/usr/local/bin/mpirun'
    mb = '/home/cfriedline/src/mrbayes_3.2.1/src/mb'
    procs = 8
    project_dir = '/home/cfriedline/projects/bsim3'
    n_gen = 1000000
    mb_new = mb

if num_taxa > 4:
    mb = mb_new
    n_gen = 500000

description = ["%d-%d-%d-%d-%.1f-%.1f" % (num_taxa, bits, num_states, runs, rate, brlen),
               "%d samples, %d bit encoding, %d states, %d runs, rate=%f, brlen=%f" % (num_taxa, bits, num_states, runs, rate, brlen)]
run_dir_name = datetime.datetime.now().strftime("%m%d%y_%H%M%S")
run_dir_name = description[0]
result_dir = app.create_dir(os.path.join(project_dir, "results"))
run_dir = app.create_dir(os.path.join(result_dir, run_dir_name))
out_dir = app.create_dir(os.path.join(run_dir, "out"))
log_dir = app.create_dir(os.path.join(run_dir, "log"))
log_file = open(os.path.join(log_dir, "log.txt"), "w")
app.log_file = log_file

desc = open(os.path.join(run_dir, "readme.txt"), "w")
with desc:
    desc.write("%s\n" % description[1])

range_fh = open(os.path.join(log_dir, "range_time.txt"), "w")


def create_R():
    r = robjects.r
    r('library(phangorn)')
    r('library(vegan)')

    r("""
        gapweight = function(x, min, max, states) {
	        return(round(((x-min)*(states-1))/(max-min)))
        }
    """)

    r("""
        get_valid_matrix = function(cols, numstates, rate) {
            found = 0
            while (found < cols) {
                root = sample(numstates)[1]
                temp = rTraitDisc(tree, model="ER", k=numstates, states=1:numstates, rate=rate, root.value=root)
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


def print_matrix(matrix, prefix, sample_names, num_cols, run, numpy, roots, i):
    fh = open(os.path.join(log_dir, "%s_%d_%d_%d.txt" % (prefix, num_cols, run, i)), "w")
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
            fh.write("%d\t%s\n" % (i, '\t'.join([str(int(elem)) for elem in range[0:2]])))


def print_matrices(data, gap, abund, ranges, num_cols, run, sample_names, roots, i):
    print_matrix(data, "data", sample_names, num_cols, run, True, roots, i)
    print_matrix(gap, "gap", sample_names, num_cols, run, True, None, i)
    print_matrix(abund, "abund", sample_names, num_cols, run, True, None, i)
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
    return '\t'.join([str(int(elem)) for elem in tuple])


@clockit
def store_valid_matrix_data(r, taxa_tree, num_cols, num_states):
    print "Storing valid matrix data in R session"
    r('matrix_data = get_valid_matrix(%d, %d, %f)' % (num_cols, num_states, rate))
    r('data = matrix_data[[1]]')
    r('roots = matrix_data[[2]]')
    r("data = t(apply(data, 1, as.numeric))")
    robjects.globalenv['colnames'] = sorted(taxa_tree.taxon_set.labels())
    r('colnames(data) = colnames')


def get_unifrac_pcoa(tree, matrix, rownames):
    test = app.get_unifrac_pcoa_tree(matrix, rownames)
    diff = app.calculate_differences_r(tree, test)
    return test, diff


def get_unifrac_cluster(tree, matrix, rownames):
    test = app.get_py_unifrac_cluster(matrix, rownames)
    diff = app.calculate_differences_r(tree, test)
    return test, diff


def get_unifrac_nj(tree, matrix, rownames):
    test = app.get_unifrac_nj(matrix, rownames)
    diff = app.calculate_differences_r(tree, test)
    return test, diff


def get_bc_pcoa(tree, abund, sample_names):
    test = app.get_bc_pcoa_tree(abund, sample_names)
    diff = app.calculate_differences_r(tree, test)
    return test, diff


def get_bc_cluster(tree, abund, sample_names):
    test = app.get_bc_cluster(abund, sample_names)
    diff = app.calculate_differences_r(tree, test)
    return test, diff


def get_bc_nj(tree, abund, sample_names):
    test = app.get_bc_nj(abund, sample_names)
    diff = app.calculate_differences_r(tree, test)
    return test, diff


def get_gap_weight(abund, min, max, states):
    try:
        return round(((float(abund) - min) * (states - 1)) / (max - min))
    except:
        print max, min

def fake(i):
    print "fake", i


def get_abundance_ranges(r, gap, num_states):
    num = len(gap[0])
    ranges = []
    for i in xrange(num):
        ranges.append(app.get_random_min_and_max_from_gamma(gamma_shape, gamma_scale, app.compute_smallest_max()))

    for i, range in enumerate(ranges):
        for i in xrange(num_states):
            range[2].append([-1] * 2)
        r("range=matrix(%d:%d)" % (range[0], range[1]))
        r('weights=apply(range, 2, gapweight, min=min(range), max=max(range), states=%d)' % num_states)
        weights = r('table(weights)')
        min = range[0]
        for i, count in enumerate(weights):
            range[2][i][0] = min
            range[2][i][1] = min + count - 1
            min += count
    return ranges

def get_column(matrix, i):
    return [row[i] for row in matrix]

@clockit
def create_abund_pool(r, gap, num_states):
    pool = mp.Pool()
    results = []
    for i in xrange(distance_runs):
        results.append(pool.apply_async(get_abundance_ranges, (r, gap, num_states)))
    pool.close()
    pool.join()
    abund_pool = []
    for r in results:
        abund_pool.append(r.get())
    return abund_pool

@clockit
def create_abund_pool_from_states(r, data):
    data2 = numpy.ndarray(data.shape)
    for i in xrange(len(data)):
        for j in xrange(len(data[i])):
            data2[i][j] = data[i][j] - 1
    pool = mp.Pool()
    results = []
    for i in xrange(distance_runs):
        results.append(pool.apply_async(get_abundance_ranges, (r, data2, num_states)))
    pool.close()
    pool.join()
    abund_pool = []
    for r in results:
        abund_pool.append(r.get())
    return abund_pool, data2

def run_mr_bayes(run, i, disc, sample_names, orig_tree):
    mb_tree = app.run_mrbayes(str(run) + "-" + str(i), disc,
                           sample_names, disc.ncol,
                           n_gen, mpi,
                           mb, procs,
                           None, run_dir,
                           len(sample_names), "d")
    diffs = app.calculate_differences_r(orig_tree, mb_tree)
    return mb_tree, diffs


@clockit
def run_simulation(r, taxa_tree, num_cols, run, out_file, dist_file, abundance_from_states):
    assert isinstance(taxa_tree, dendropy.Tree)
    r('tree = rtree(%d, rooted=F)' % num_taxa)
    r('tree$edge.length = rep(%f, length(tree$edge.length))' % brlen)
    tree = app.ape_to_dendropy(r['tree'])
    store_valid_matrix_data(r, taxa_tree, num_cols, num_states)
    data = numpy.asarray(r['data'])
    roots = list(r['roots'])
    sample_names = numpy.asarray(r('rownames(data)'))
    print_sample_trees(r, tree, num_taxa, num_cols, run)
    print_state_distribution(data, num_cols, run, sample_names, dist_file)
    gap = None
    if not abundance_from_states:
        ranges = get_column_ranges(data)
        data, ranges = curate_data_matrix(r, data, ranges)
        gap = app.restandardize_matrix(data, ranges, num_states)


    u_pcoa_diffs_total = []
    u_cluster_diffs_total = []
    u_nj_diffs_total = []

    w_pcoa_diffs_total = []
    w_cluster_diffs_total = []
    w_nj_diffs_total = []

    bc_pcoa_diffs_total = []
    bc_cluster_diffs_total = []
    bc_nj_diffs_total = []
    mb_diffs_total = []

    abund_pool = None
    if not abundance_from_states:
        abund_pool = create_abund_pool(r, gap, num_states)
    else:
        abund_pool, gap = create_abund_pool_from_states(r, data)


    for i in xrange(distance_runs):
        abund_ranges = abund_pool[i]
        abund = app.get_abundance_matrix(gap, abund_ranges, "gamma", num_states)
        print_matrices(data, gap, abund, abund_ranges, num_cols, run, sample_names, roots, i)
        (u_matrix, u_names), (w_matrix, w_names) = app.calculate_unifrac(abund, sample_names, taxa_tree)

        # unifrac tests
        u_pcoa_tree, u_pcoa_diffs = get_unifrac_pcoa(tree, u_matrix, u_names)
        u_pcoa_diffs_total.append(u_pcoa_diffs)
        u_cluster_tree, u_cluster_diffs = get_unifrac_cluster(tree, u_matrix, u_names)
        u_cluster_diffs_total.append(u_cluster_diffs)
        u_nj_tree, u_nj_diffs = get_unifrac_nj(tree, u_matrix, u_names)
        u_nj_diffs_total.append(u_nj_diffs)

        w_pcoa_tree, w_pcoa_diffs = get_unifrac_pcoa(tree, w_matrix, w_names)
        w_pcoa_diffs_total.append(w_pcoa_diffs)
        w_cluster_tree, w_cluster_diffs = get_unifrac_cluster(tree, w_matrix, w_names)
        w_cluster_diffs_total.append(w_cluster_diffs)
        w_nj_tree, w_nj_diffs = get_unifrac_nj(tree, w_matrix, w_names)
        w_nj_diffs_total.append(w_nj_diffs)

        # bray-curtis tests
        bc_pcoa_tree, bc_pcoa_diffs = get_bc_pcoa(tree, abund, sample_names)
        bc_pcoa_diffs_total.append(bc_pcoa_diffs)
        bc_cluster_tree, bc_cluster_diffs = get_bc_cluster(tree, abund, sample_names)
        bc_cluster_diffs_total.append(bc_cluster_diffs)
        bc_nj_tree, bc_nj_diffs = get_bc_nj(tree, abund, sample_names)
        bc_nj_diffs_total.append(bc_nj_diffs)

        if abundance_from_states:
            gap_from_abund = app.restandardize_matrix(abund, abund_ranges, num_states)
            disc = app.get_discrete_matrix_from_standardized(gap_from_abund, bits, sample_names)
            mb_tree, mb_diffs = run_mr_bayes(run, i, disc, sample_names, tree)
            mb_diffs_total.append(mb_diffs)


    # mrbayes
    disc = None
    if not abundance_from_states:
        disc = app.get_discrete_matrix_from_standardized(gap, bits, sample_names)
        mb_tree, mb_diffs = run_mr_bayes(run, 0, disc, sample_names, tree)
        mb_diffs_total.append(mb_diffs)

    # output
    for i in xrange(distance_runs):
        if not abundance_from_states:
            out_file.write("%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (run, num_cols, i,
                                                                     get_tab_string(mb_diffs_total[0]),
                                                                     get_tab_string(u_pcoa_diffs_total[i]),
                                                                     get_tab_string(u_cluster_diffs_total[i]),
                                                                     get_tab_string(u_nj_diffs_total[i]),
                                                                     get_tab_string(w_pcoa_diffs_total[i]),
                                                                     get_tab_string(w_cluster_diffs_total[i]),
                                                                     get_tab_string(w_nj_diffs_total[i]),
                                                                     get_tab_string(bc_pcoa_diffs_total[i]),
                                                                     get_tab_string(bc_cluster_diffs_total[i]),
                                                                     get_tab_string(bc_nj_diffs_total[i])))
        else:
            out_file.write("%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (run, num_cols, i,
                                                                     get_tab_string(mb_diffs_total[i]),
                                                                     get_tab_string(u_pcoa_diffs_total[i]),
                                                                     get_tab_string(u_cluster_diffs_total[i]),
                                                                     get_tab_string(u_nj_diffs_total[i]),
                                                                     get_tab_string(w_pcoa_diffs_total[i]),
                                                                     get_tab_string(w_cluster_diffs_total[i]),
                                                                     get_tab_string(w_nj_diffs_total[i]),
                                                                     get_tab_string(bc_pcoa_diffs_total[i]),
                                                                     get_tab_string(bc_cluster_diffs_total[i]),
                                                                     get_tab_string(bc_nj_diffs_total[i])))
        out_file.flush()


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
    return "run\tcols\titer\t"\
           "mb_topo\tmb_sym\tmb_path\t"\
           "u_pcoa_topo\tu_pcoa_symm\tu_pcoa_path\t"\
           "u_cluster_topo\tu_cluster_symm\tu_cluster_path\t"\
           "u_nj_topo\tu_nj_symm\tu_nj_path\t"\
           "w_pcoa_topo\tw_pcoa_symm\tw_pcoa_path\t"\
           "w_cluster_topo\tw_cluster_symm\tw_cluster_path\t"\
           "w_nj_topo\tw_nj_symm\tw_nj_path\t"\
           "bc_pcoa_topo\tbc_pcoa_symm\tbc_pcoa_path\t"\
           "bc_cluster_topo\tbc_cluster_symm\tbc_cluster_path\t"\
           "bc_nj_topo\tbc_nj_symm\tbc_nj_path"

def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--abundance_from_states", help="create abundance matrix directly from states", action="store_true", default=False)

    args = p.parse_args()
    args.abundance_from_states = True
    return args

@clockit
def main():
    args = get_args()
    r = create_R()
    out_file = open(os.path.join(out_dir, "out.txt"), "w")
    out_file.write("%s\n" % get_header())
    dist_file = open(os.path.join(out_dir, "dists.txt"), "w")
    for col in cols:
        tree = app.create_tree(col, "T")
        print_taxa_tree(tree, col)
        for run in xrange(runs):
            run_simulation(r, tree, col, run, out_file, dist_file, args.abundance_from_states)
    out_file.close()
    dist_file.close()

if __name__ == '__main__':
    main()
    range_fh.close()
    print "Done!"


