import traceback

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
import argparse
import celery
from celery import Celery
from runsimulation4 import *  #workaround for celery
import tempfile

numpy2ri.activate()
sys.setrecursionlimit(1000000)

cluster = Celery('runsimulation4')
cluster.config_from_object('localconfig')

hostname = os.uname()[1]
project_dir = "/Users/chris/projects/asmw4"
n_gen = 100000
mpi = "/usr/local/openmpi/bin/mpirun"
mb_old = "/Users/chris/src/mrbayes-3.1.2/mb"
mb_new = "/Users/chris/src/mrbayes_3.2.1/src/mb"
mb = mb_old
procs = 8
num_taxa = 8
num_states = 8
bits = 3
rate = 1.0
gamma_shape = 1
gamma_scale = 1000
sigma = 0.5  #for BM

if 'godel' in hostname:
    mpi = '/usr/global/openmpi-1.5.3-w-psm/bin/mpirun'
    mb = '/home/cfriedline/src/mrbayes_3.2.1/src/mb'
    procs = 8
    project_dir = '/home/cfriedline/projects/asmw4'
    mb_new = mb
    cluster.config_from_object('clusterconfig')
elif 'phylogeny' in hostname:
    mpi = '/usr/local/bin/mpirun'
    mb = '/home/cfriedline/src/mrbayes_3.2.1/src/mb'
    procs = 8
    project_dir = '/home/cfriedline/projects/asmw2'
    n_gen = 1000000
    mb_new = mb

if num_taxa > 4:
    mb = mb_new
    n_gen = 500000


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

    r("""
        get_continuous_matrix = function(cols, shape, scale, sigma) {
            r = rgamma(1, shape=shape, scale=scale)
            m = replicate(cols, rTraitCont(tree, ancestor=F, sigma=r*0.5, root=r))
            m = round(as.matrix(m))
            return(list(m, r))
        }

    """)

    return r


def get_column_ranges(data):
    ranges = []
    for col in data.T:
        minval = min(col)
        if minval < 0:
            minval = 0
        ranges.append((minval, max(col)))
    return ranges


def get_bad_cols(ranges):
    bad = []
    for i, range in enumerate(ranges):
        if range[0] == range[1]:
            bad.append(i)
    return bad


def print_matrix(matrix, prefix, sample_names, num_cols, tree_num, numpy, roots, i, filedata):
    fh = open(os.path.join(filedata['log_dir'], "%s_%d_%d_%d.txt" % (prefix, num_cols, tree_num, i)), "w")
    m = matrix
    if not numpy:
        m = numpy.asarray(matrix)
    with fh:
        if roots:
            fh.write("#\t%s\n" % '\t'.join([str(int(elem)) for elem in roots]))
        for i, row in enumerate(m):
            fh.write("%s\t%s\n" % (sample_names[i], '\t'.join([str(int(elem)) for elem in row])))


def print_ranges(ranges, prefix, num_cols, run, filedata):
    fh = open(os.path.join(filedata['log_dir'], "%s_%d_%d.txt" % (prefix, num_cols, run)), "w")
    with fh:
        for i, range in enumerate(ranges):
            fh.write("%d\t%s\n" % (i, '\t'.join([str(int(elem)) for elem in range[0:2]])))


def print_matrices(data, gap, abund, ranges, gap_abund, new_ranges, cont_abund, cont_gap, cont_ranges, sub_abund,
                   sub_abund_ranges, num_cols, tree_num, sample_names, roots, i,
                   filedata):
    print_matrix(data, "data", sample_names, num_cols, tree_num, True, roots, i, filedata)
    print_matrix(gap, "gap", sample_names, num_cols, tree_num, True, None, i, filedata)
    print_matrix(abund, "abund", sample_names, num_cols, tree_num, True, None, i, filedata)
    print_matrix(gap_abund, "gap_abund", sample_names, num_cols, tree_num, True, None, i, filedata)
    print_matrix(cont_abund, "cont_abund", sample_names, num_cols, tree_num, True, None, i, filedata)
    print_matrix(cont_gap, "cont_gap", sample_names, num_cols, tree_num, True, None, i, filedata)
    print_matrix(sub_abund, "sub_abund", sample_names, num_cols, tree_num, True, None, i, filedata)

    print_ranges(ranges, "ranges", num_cols, tree_num, filedata)
    print_ranges(new_ranges, "new_ranges", num_cols, tree_num, filedata)
    print_ranges(cont_ranges, "cont_ranges", num_cols, tree_num, filedata)
    print_ranges(sub_abund_ranges, "sub_ranges", num_cols, tree_num, filedata)


@clockit
def print_state_distribution(name_key, data, num_cols, tree_num, sample_names, dist_file):
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

    dist_file.write("%d, %d, %s\n" % (num_cols, tree_num, name_key))
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


def print_sample_trees(r, tree, num_taxa, num_cols, tree_num, filedata):
    assert isinstance(tree, dendropy.Tree)
    pdf_file = os.path.join(filedata['log_dir'], "tree_%d_%d_%d.pdf" % (num_taxa, num_cols, tree_num))
    text_file = os.path.join(filedata['log_dir'], "tree_%d_%d_%d.txt" % (num_taxa, num_cols, tree_num))
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
def store_valid_matrix_data(r, taxa_tree, num_cols, num_states, rate):
    print "Storing valid matrix data in R session"
    r('matrix_data = get_valid_matrix(%d, %d, %f)' % (num_cols, num_states, rate))
    r('cont_matrix_data = get_continuous_matrix(%d, %f, %f, %f)' % (num_cols, gamma_shape, gamma_scale, sigma))
    r('data = matrix_data[[1]]')
    r('roots = matrix_data[[2]]')
    r('data_cont = cont_matrix_data[[1]]')
    r('data_cont = ifelse(data_cont < 0, 0, data_cont)')
    r('roots_cont = cont_matrix_data[[2]]')
    r("data = t(apply(data, 1, as.numeric))")
    robjects.globalenv['colnames'] = sorted(taxa_tree.taxon_set.labels())
    r('colnames(data) = colnames')
    r('colnames(data_cont) = colnames')


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
        for j in xrange(num_states):
            range[2].append([-1] * 2)
        r("range=matrix(%d:%d)" % (range[0], range[1]))
        r('weights=apply(range, 2, gapweight, min=min(range), max=max(range), states=%d)' % num_states)
        weights = r('table(weights)')
        min_val = range[0]
        for j, count in enumerate(weights):
            range[2][j][0] = min_val
            range[2][j][1] = min_val + count - 1
            min_val += count
    return ranges


def get_column(matrix, i):
    return [row[i] for row in matrix]


def create_abund_pool_from_states(r, data):
    data2 = numpy.ndarray(data.shape)
    for i in xrange(len(data)):
        for j in xrange(len(data[i])):
            data2[i][j] = data[i][j] - 1
    abund_pool = list()
    abund_pool.append(get_abundance_ranges(r, data2, num_states))
    return abund_pool, data2


def run_mr_bayes(key, tree_num, i, disc, sample_names, orig_tree, filedata, mrbayes_timeout):
    mb_tree = app.run_mrbayes(key, str(tree_num) + "-" + str(i), disc,
                              sample_names, disc.ncol,
                              n_gen, mpi,
                              mb, procs,
                              None, filedata['run_dir'],
                              len(sample_names), "d", filedata['hostfile'], mrbayes_timeout)
    diffs = app.calculate_differences_r(orig_tree, mb_tree)
    return mb_tree, diffs


def is_float(num):
    try:
        float(num)
        return True
    except ValueError:
        return False


@cluster.task()
def run_simulation(taxa_tree, taxa_tree_fixedbr, sample_tree, tree_num, num_cols, out_file, dist_file,
                   abundance_from_states, filedata, brlen, mrbayes_timeout):
    assert isinstance(taxa_tree, dendropy.Tree)
    out_file_name = out_file
    dist_file_name = dist_file
    is_celery = False
    r = create_R()
    if not out_file:
        is_celery = True
        out_file = tempfile.NamedTemporaryFile(dir=".", delete=False)
        out_file_name = out_file.name
    if not dist_file:
        dist_file = tempfile.NamedTemporaryFile(dir=".", delete=False)
        dist_file_name = dist_file.name

    r('temp=rtree(%d, rooted=F)' % num_taxa)
    r('edges = runif(length(temp$edge.length), min=0.1)')
    r('tree=read.tree(text="%s")' % sample_tree)
    r('tree$edge.length=edges')
    if is_float(brlen):
        r('tree$edge.length = rep(%f, length(tree$edge.length))' % float(brlen))

    tree = app.ape_to_dendropy(r['tree'])
    print tree
    store_valid_matrix_data(r, taxa_tree, num_cols, num_states, filedata['rate'])
    data = numpy.asarray(r['data'])
    roots = list(r['roots'])
    sample_names = numpy.asarray(r('rownames(data)'))
    print_sample_trees(r, tree, num_taxa, num_cols, tree_num, filedata)
    print_state_distribution("state", data, num_cols, tree_num, sample_names, dist_file)
    gap = None
    if not abundance_from_states:
        print "not abundance from states!"
        ranges = get_column_ranges(data)
        data, ranges = curate_data_matrix(r, data, ranges)
        gap = app.restandardize_matrix(data, ranges, num_states)

    abund_pool, gap = create_abund_pool_from_states(r, data)

    abund_ranges = abund_pool[0]

    abund = app.get_abundance_matrix(gap, abund_ranges, "gamma", num_states)
    sub_abund = app.subsample_abundance_matrix(abund, 10)
    sub_abund_ranges = get_column_ranges(numpy.array(sub_abund))
    gap_from_sub = app.restandardize_matrix(sub_abund, sub_abund_ranges, num_states)
    print_state_distribution("sub", gap_from_sub, num_cols, tree_num, sample_names, dist_file)

    #continus matrix to gap to discrete matrix
    cont_abund = app.get_continuous_abundance_matrix(r)
    cont_ranges = get_column_ranges(numpy.array(cont_abund))
    cont_gap = app.restandardize_matrix(cont_abund, cont_ranges, num_states)
    print_state_distribution("cont", cont_gap, num_cols, tree_num, sample_names, dist_file)

    (u_matrix, u_names), (w_matrix, w_names) = app.calculate_unifrac(abund, sample_names, taxa_tree)
    (u_matrix_norm, u_names_norm), (w_matrix_norm, w_names_norm) = app.calculate_unifrac(abund, sample_names,
                                                                                         taxa_tree_fixedbr)

    # unifrac tests
    u_pcoa_tree, u_pcoa_diffs = get_unifrac_pcoa(tree, u_matrix, u_names)
    u_cluster_tree, u_cluster_diffs = get_unifrac_cluster(tree, u_matrix, u_names)
    u_nj_tree, u_nj_diffs = get_unifrac_nj(tree, u_matrix, u_names)
    w_pcoa_tree, w_pcoa_diffs = get_unifrac_pcoa(tree, w_matrix, w_names)
    w_cluster_tree, w_cluster_diffs = get_unifrac_cluster(tree, w_matrix, w_names)
    w_nj_tree, w_nj_diffs = get_unifrac_nj(tree, w_matrix, w_names)

    # unifrac tests (normalized branch lengths on taxa tree)
    u_pcoa_tree_norm, u_pcoa_diffs_norm = get_unifrac_pcoa(tree, u_matrix_norm, u_names_norm)
    u_cluster_tree_norm, u_cluster_diffs_norm = get_unifrac_cluster(tree, u_matrix_norm, u_names_norm)
    u_nj_tree_norm, u_nj_diffs_norm = get_unifrac_nj(tree, u_matrix_norm, u_names_norm)
    w_pcoa_tree_norm, w_pcoa_diffs_norm = get_unifrac_pcoa(tree, w_matrix_norm, w_names_norm)
    w_cluster_tree_norm, w_cluster_diffs_norm = get_unifrac_cluster(tree, w_matrix_norm, w_names_norm)
    w_nj_tree_norm, w_nj_diffs_norm = get_unifrac_nj(tree, w_matrix_norm, w_names_norm)

    # bray-curtis tests
    bc_pcoa_tree, bc_pcoa_diffs = get_bc_pcoa(tree, abund, sample_names)
    bc_cluster_tree, bc_cluster_diffs = get_bc_cluster(tree, abund, sample_names)
    bc_nj_tree, bc_nj_diffs = get_bc_nj(tree, abund, sample_names)

    new_ranges = get_column_ranges(numpy.array(abund))
    gap_from_abund = app.restandardize_matrix(abund, new_ranges, num_states)
    disc = app.get_discrete_matrix_from_standardized(gap_from_abund, bits, sample_names)
    disc2 = app.get_discrete_matrix_from_standardized2(gap_from_sub, num_states, sample_names)
    cont_disc = app.get_discrete_matrix_from_standardized2(cont_gap, num_states, sample_names)

    print_matrices(data, gap, abund, abund_ranges, gap_from_abund, new_ranges, cont_abund, cont_gap, cont_ranges,
                   sub_abund, sub_abund_ranges, num_cols, tree_num, sample_names, roots, 0, filedata)

    mb_tree, mb_diffs = run_mr_bayes("state", tree_num, 0, disc, sample_names, tree, filedata, mrbayes_timeout)
    mb_tree2, mb_diffs2 = run_mr_bayes("sub", tree_num, 0, disc2, sample_names, tree, filedata, mrbayes_timeout)
    mb_tree_cont, mb_diffs_cont = run_mr_bayes("cont", tree_num, 0, cont_disc, sample_names, tree, filedata,
                                               mrbayes_timeout)

    # output
    try:
        out_file.write("%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                       (tree_num, num_cols,
                        get_tab_string(mb_diffs),
                        get_tab_string(mb_diffs2),
                        get_tab_string(mb_diffs_cont),
                        get_tab_string(u_pcoa_diffs),
                        get_tab_string(u_cluster_diffs),
                        get_tab_string(u_nj_diffs),
                        get_tab_string(w_pcoa_diffs),
                        get_tab_string(w_cluster_diffs),
                        get_tab_string(w_nj_diffs),
                        get_tab_string(u_pcoa_diffs_norm),
                        get_tab_string(u_cluster_diffs_norm),
                        get_tab_string(u_nj_diffs_norm),
                        get_tab_string(w_pcoa_diffs_norm),
                        get_tab_string(w_cluster_diffs_norm),
                        get_tab_string(w_nj_diffs_norm),
                        get_tab_string(bc_pcoa_diffs),
                        get_tab_string(bc_cluster_diffs),
                        get_tab_string(bc_nj_diffs)))
    except:
        print "Unexpected error in writing output", sys.exc_info()
        sys.exit()

    if is_celery:
        out_file.close()
        dist_file.close()

    return (out_file_name, dist_file_name, tree, sample_names, data, gap, abund, abund_ranges,
            gap_from_abund, new_ranges,
            mb_tree, mb_diffs,
            (u_matrix, u_names), (w_matrix, w_names),
            (u_matrix_norm, u_names_norm), (w_matrix_norm, w_names_norm))


def print_taxa_tree(tree, num_cols, filedata):
    assert isinstance(tree, dendropy.Tree)
    r = robjects.r
    r("temp = read.tree(text='%s;')" % tree.as_newick_string())
    pdf_file = os.path.join(filedata['log_dir'], "taxa_tree_%d.pdf" % num_cols)
    text_file = os.path.join(filedata['log_dir'], "taxa_tree_%d.txt" % num_cols)
    r("pdf(file='%s')" % pdf_file)
    r("plot(temp, root.edge=T)")
    r('dev.off()')

    fh = open(text_file, "w")
    with fh:
        fh.write("%s;\n" % tree.as_newick_string())


def get_header():
    return "tree_num\tcols\t" \
           "mb_topo\tmb_symm\tmb_path\t" \
           "mb_topo2\tmb_symm2\tmb_path2\t" \
           "mb_topo_cont\tmb_symm_cont\tmb_path_cont\t" \
           "u_pcoa_topo\tu_pcoa_symm\tu_pcoa_path\t" \
           "u_cluster_topo\tu_cluster_symm\tu_cluster_path\t" \
           "u_nj_topo\tu_nj_symm\tu_nj_path\t" \
           "w_pcoa_topo\tw_pcoa_symm\tw_pcoa_path\t" \
           "w_cluster_topo\tw_cluster_symm\tw_cluster_path\t" \
           "w_nj_topo\tw_nj_symm\tw_nj_path\t" \
           "u_pcoa_topo_norm\tu_pcoa_symm_norm\tu_pcoa_path_norm\t" \
           "u_cluster_topo_norm\tu_cluster_symm_norm\tu_cluster_path_norm\t" \
           "u_nj_topo_norm\tu_nj_symm_norm\tu_nj_path_norm\t" \
           "w_pcoa_topo_norm\tw_pcoa_symm_norm\tw_pcoa_path_norm\t" \
           "w_cluster_topo_norm\tw_cluster_symm_norm\tw_cluster_path_norm\t" \
           "w_nj_topo_norm\tw_nj_symm_norm\tw_nj_path_norm\t" \
           "bc_pcoa_topo\tbc_pcoa_symm\tbc_pcoa_path\t" \
           "bc_cluster_topo\tbc_cluster_symm\tbc_cluster_path\t" \
           "bc_nj_topo\tbc_nj_symm\tbc_nj_path"


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--tree_file", help="tree file, one newick per line")
    p.add_argument("--cols", help="cols in the matrix", type=int)
    p.add_argument("--run", help="run id (0-9 for 10 runs)")
    p.add_argument("--brlen", help="branch lengths for sample tree", default=0.5)
    p.add_argument("--project_dir", help="root dir for the project", default="../asmw8")
    p.add_argument('--mrbayes_timeout', help="timeout for mrbayes instance", type=float)
    p.add_argument('--rate', help="rate to simulate matrix", type=float, default=1.0)
    p.add_argument('--celery', help="use celery parallelism", action='store_true')
    p.add_argument('--task', help="the task to run", choices=['full', 'bm_sigma', 'rate', 'subsample'], required=True,
                   action='append')

    args = p.parse_args()

    if not args.tree_file or not len(sys.argv) > 1:
        p.print_help()
        exit()

    args.abundance_from_states = True
    return args


def get_trees(tree_file, out_dir):
    trees = []
    for line in open(tree_file):
        trees.append(line.rstrip())
    return trees


def create_file_data(args):
    data = dict()
    data['project_dir'] = args.project_dir
    data['rate'] = args.rate
    data['description'] = ["%d-%d-%d-%.1f-%s" % (num_taxa, bits, num_states, data['rate'], args.brlen),
                           "%d samples, %d bit encoding, %d states, rate=%.2f, brlen=%s" % (
                               num_taxa, bits, num_states, rate, args.brlen)]
    data['run_dir_name'] = datetime.datetime.now().strftime("%m%d%y_%H%M%S")
    data['run_dir_name'] = str(args.cols) + "-" + str(args.run)
    data['result_dir'] = app.create_dir(os.path.join(args.project_dir, "results"))
    data['run_dir'] = app.create_dir(os.path.join(data['result_dir'], data['run_dir_name']))
    data['out_dir'] = app.create_dir(os.path.join(data['run_dir'], "out"))
    data['log_dir'] = app.create_dir(os.path.join(data['run_dir'], "log"))
    data['log_file'] = open(os.path.join(data['log_dir'], "log.txt"), "w")
    app.log_file = data['log_file']

    data['desc'] = open(os.path.join(data['run_dir'], "readme.txt"), "w")
    with data['desc']:
        data['desc'].write("%s\n" % data['description'][1])

    data['range_fh'] = open(os.path.join(data['log_dir'], "range_time.txt"), "w")
    data['hostfile'] = None

    if os.path.isfile("hostfile"):
        data['hostfile'] = os.path.abspath("hostfile")

    data['rate'] = args.rate
    data['celery'] = args.celery
    return data


def create_uniform_brlen_tree(taxa_tree, brlen):
    assert isinstance(taxa_tree, dendropy.Tree)

    r = robjects.r
    temp_name = app._get_random_string(10)
    r("%s = read.tree(text='%s;')" % (temp_name, taxa_tree.as_newick_string()))
    r('%s$edge.length = rep(%f, length(%s$edge.length))' % (temp_name, brlen, temp_name))
    return app.ape_to_dendropy(r[temp_name])


def write_file_data(filedata):
    f = os.path.join(filedata['project_dir'], 'params.log')
    with open(f, "w") as out:
        for k, v in filedata.items():
            out.write("%s: %s\n" % (k, v))


def run_full_simulation(sample_trees, filedata, args, taxa_tree, taxa_tree_fixedbr, col, out_file, dist_file):
    submit_count = 0
    celery_results = []
    for tree_num, sample_tree in enumerate(sample_trees):

        if filedata['celery']:
            print tree_num
            res = run_simulation.delay(taxa_tree, taxa_tree_fixedbr, sample_tree, tree_num, col, None, None,
                                       args.abundance_from_states,
                                       filedata, args.brlen, args.mrbayes_timeout)
            celery_results.append(res)
            submit_count += 1
            if submit_count == 1:
                break
        else:
            run_simulation(taxa_tree, taxa_tree_fixedbr, sample_tree, tree_num, col, out_file, dist_file,
                           args.abundance_from_states,
                           filedata, args.brlen, args.mrbayes_timeout)
            if tree_num == 10:
                break

    if filedata['celery']:
        for i, res in enumerate(celery_results):
            res = res.get()
            print "getting results for task %d" % i
            o = open(res[0])
            d = open(res[1])
            for line in o:
                out_file.write("%s" % line)
            for line in d:
                dist_file.write("%s" % line)
            o.close()
            d.close()
            os.unlink(res[0])
            os.unlink(res[1])
    out_file.close()
    dist_file.close()
    filedata['range_fh'].close()


def run_bm_sigma_simulation():
    pass


def run_rate_simulation():
    pass


def run_subsample_simulation():
    pass


@clockit
def main():
    args = get_args()
    filedata = create_file_data(args)
    write_file_data(filedata)
    #    sys.stdout = open(os.path.basename("%s_stdout.txt" % args.project_dir), "w")
    sample_trees = get_trees(args.tree_file, filedata['out_dir'])
    out_file = open(os.path.join(filedata['out_dir'], "out.txt"), "w", 0)
    out_file.write("%s\n" % get_header())
    dist_file = open(os.path.join(filedata['out_dir'], "dists.txt"), "w", 0)
    col = args.cols
    taxa_tree = app.create_tree(col, "T")
    taxa_tree_fixedbr = create_uniform_brlen_tree(taxa_tree, 0.5)
    print_taxa_tree(taxa_tree, col, filedata)
    if 'full' in args.task:
        run_full_simulation(sample_trees, filedata, args, taxa_tree, taxa_tree_fixedbr, col, out_file, dist_file)

    if 'bm_sigma' in args.task:
        run_bm_sigma_simulation()

    if 'rate' in args.task:
        run_rate_simulation()

    if 'subsample' in args.task:
        run_subsample_simulation()


if __name__ == '__main__':
    main()
    print "Done!"


