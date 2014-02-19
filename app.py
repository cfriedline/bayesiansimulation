import os
import string
import traceback

import rpy2.rinterface as rinterface


rinterface.set_initoptions(('rpy2', '--vanilla', '--max-ppsize=500000', '--quiet'))

import numpy
from subprocess import Popen, PIPE, STDOUT
import dendropy
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import tempfile
from dendropy.interop import ape
import rpy2
from rpy2.robjects import numpy2ri
import math
import stopwatch
import multiprocessing as mp
from multiprocessing import Pool, current_process, Manager
from threading import Timer
import platform
import random
from itertools import izip
#mp_logger = mp.log_to_stderr()
#mp_logger.setLevel(mp.SUBDEBUG)

numpy2ri.activate()
from cogent import LoadTree
from cogent.maths.unifrac.fast_unifrac import fast_unifrac, UNIFRAC_DIST_MATRIX
from scipy import stats

_author_ = 'chris'
log_file = None


def clockit(func):
    """
    Function decorator that times the evaluation of *func* and prints the
    execution time. Modified from stopwatch package
    """

    def new(*args, **kw):
        t = stopwatch.Timer()
        retval = func(*args, **kw)
        t.stop()
        output = '\t%s in %s' % (func.__name__, t)
        print output
        if log_file is not None:
            log_file.write("%s\n" % output)
            log_file.flush()
        del t
        return retval

    return new


def compute_smallest_max():
    return 20


@clockit
def get_py_unifrac_cluster(uni_matrix, uni_names):
    """
    clusters unifrac distance matrix created in pycogent
    @param uni_matrix: pycogent dist matrix
    @param uni_names: vector of row names
    @return a dendropy tree of the cluster
    @rtype: dendropy.Tree
    """
    assert isinstance(uni_matrix, numpy.ndarray)
    r = robjects.r
    robjects.globalenv['unimatrix'] = robjects.conversion.py2ri(numpy.array(uni_matrix))
    robjects.globalenv['uninames'] = robjects.conversion.py2ri(numpy.array(uni_names))
    r('rownames(unimatrix)=uninames')
    r('colnames(unimatrix)=uninames')
    r('unimatrix = as.dist(unimatrix)')
    r("unimatrixclust = hclust(unimatrix, 'ave')")
    tree = r('multi2di(as.phylo(unimatrixclust))')
    return ape_to_dendropy(tree)


@clockit
def calculate_unifrac(abund, sample_names, taxa_tree):
    """
    calculates the unifrac distance between samples both
    weighted and unweighted
    @param abund: the abundance matrix
    @param sample_names: the sample names
    @param taxa_tree: the tree of data
    @return: (unweighted matrix, row names), (weighted matrix, row names)
    @rtype: tuple
    """
    unifrac_dict = _create_unifrac_dict(abund, sample_names, taxa_tree)
    tree = dendropy_to_cogent(taxa_tree)
    unweighted = fast_unifrac(tree, unifrac_dict, modes={UNIFRAC_DIST_MATRIX}, is_symmetric=True, weighted=False)
    un_matrix = unweighted[UNIFRAC_DIST_MATRIX][0]
    un_rows = unweighted[UNIFRAC_DIST_MATRIX][1]

    weighted = fast_unifrac(tree, unifrac_dict, modes={UNIFRAC_DIST_MATRIX}, is_symmetric=True, weighted=True)
    w_matrix = weighted[UNIFRAC_DIST_MATRIX][0]
    w_rows = weighted[UNIFRAC_DIST_MATRIX][1]
    return (un_matrix, un_rows), (w_matrix, w_rows)


def _create_unifrac_dict(abund, sample_names, taxa_tree):
    """
    creates a unifrac dictionary
    @param abund:
    @param sample_names:
    @param taxa_tree:
    @return a dictionary of {taxa:{sample:count}}
    """
    assert isinstance(taxa_tree, dendropy.Tree)
    taxa_names = sorted(taxa_tree.taxon_set.labels())
    data = {}
    for i, row in enumerate(abund):
        sample = sample_names[i]
        for j, elem in enumerate(row):
            taxa = taxa_names[j]
            if not taxa in data:
                data[taxa] = {}
            data[taxa][sample] = row[j]
    return data


def make_tree_binary(tree_string):
    """
    takes an input string and uses R to make it binary using multi2di
    @param tree_string: newick repr of the tree
    @return: a dendropy tree
    @rtype: dendropy.Tree
    """
    r = robjects.r
    robjects.globalenv['temptree'] = tree_string + ";"
    tree = r('multi2di(read.tree(text=temptree))')
    f = tempfile.NamedTemporaryFile()
    r['write.nexus'](tree, file=f.name)
    tree = dendropy.Tree.get_from_path(f.name, "nexus")
    f.close()
    return tree


def get_paralinear_cluster():
    """
    Returns the paralinear cluster from the R session
    @return: a dendropy tree
    @rtype: dendropy.Tree
    """
    r = robjects.r
    tree = r('multi2di(as.phylo(paralinear_cluster))')
    f = tempfile.NamedTemporaryFile()
    r['write.nexus'](tree, file=f.name)
    tree = dendropy.Tree.get_from_path(f.name, "nexus")
    f.close()
    return tree


def create_dir(dir):
    """
    Creates a directory if it doesn't exist
    @param dir: directory name to create (us os.path.join...)
    @return: the path
    @rtype: string
    """
    try:
        os.makedirs(dir)
    except OSError:
        pass
    return dir


@clockit
def create_R(dir):
    """
    creates the r environment
    @param dir: the directory for the output files
    """
    r = robjects.r
    importr("phangorn")
    importr("picante")
    importr("MASS")
    importr("vegan")
    r("options(expressions=500000)")
    robjects.globalenv['outfile'] = os.path.abspath(os.path.join(dir, "trees.pdf"))
    r('pdf(file=outfile, onefile=T)')
    r("par(mfrow=c(2,3))")

    r("""
        get_discrete_matrix = function(tree, needed) {
            matrix = replicate(needed, rTraitDisc(tree, model="ER", k=2,states=0:1))
            matrix = t(apply(matrix, 1, as.numeric))
            return(matrix)                    
        }
    """)

    r("""
        get_valid_triplets = function(numsamples, needed, bits) {
            #needed is needed cols, not needed triplets, so + bits b/c generating bit-lets
            tryCatch({
                found = 0
                while (found < needed) {
                    triplet = replicate(bits, rTraitDisc(tree, model="ER", k=2,states=0:1))
                    triplet = t(apply(triplet, 1, as.numeric))
                    sums = rowSums(triplet)
                    if (length(which(sums==0)) > 0 && length(which(sums==3)) > 0) {
                        if (found == 0) {
                            m = triplet
                        } else {
                            m = cbind(m, triplet)
                        }
                        found = found + bits
                    }

                }
            return(m)
            }, error = function(e){print(message(e))}, warning = function(e){print(message(e))})
        }
    """)


def is_binary_tree(tree):
    """
    Determines if a tree is binary
    @param tree: dendropy tree
    @return: true/false
    @rtype: bool
    """
    r_tree = ape.as_ape_object(tree)
    return robjects.r['is.binary.tree'](r_tree)[0]


@clockit
def _create_paralin_matrix(a):
    """
    creates a paralinear distance matrix
    @param a: a rojbects.Matrix object of the discrete character matrix
    @return: a numpy.ndarray of the distances, whether it's valid or not
    given that it may contain na/nan/inf if not all patterns in the seqs are
    observed (less of a problem with bigger matrices)
    @rtype tuple
    """
    assert isinstance(a, robjects.Matrix)
    dist = numpy.zeros(shape=(a.nrow, a.nrow))
    valid = True
    for i in xrange(a.nrow):
        for j in range(i + 1):
            d = get_paralin_distance(a.rx(i + 1, True), a.rx(j + 1, True))
            if math.isnan(d) or math.isinf(d) or d == -1:
                valid = False
            dist[i][j] = d
            dist[j][i] = d
    return dist, valid


def get_paralin_distance(seq1, seq2):
    """
    Computes the paralinear distance between two sequences from
    a discrete character matrix
    @param seq1: 0/1 string
    @param seq2: 0/1 string
    @return: the paralinear distance
    @rtype: float
    """
    j = numpy.zeros(shape=(2, 2))
    d1 = numpy.zeros(shape=(2, 2))
    d2 = numpy.zeros(shape=(2, 2))
    for pos, c in enumerate(seq1):
        c1 = int(c)
        c2 = int(seq2[pos])
        j[c1][c2] += 1
        d1[c1][c1] += 1
        d2[c2][c2] += 1
    x = abs(numpy.linalg.det(j))
    y = math.sqrt(numpy.linalg.det(d1))
    if x == 0 or y == 0:
        return -1
    z = math.sqrt(numpy.linalg.det(d2))
    d = -math.log(x / (y * z))
    return d


def _ret_valid_triplets(results):
    return results


def _get_valid_triplets(num_samples, num_triplets, bits, q):
    try:
        r = robjects.r
        name = current_process().name.replace("-", "_")
        timer = stopwatch.Timer()
        log("\trunning %s (%d cols/%d triplets), pid %d, ppid %d" % (
            name, num_triplets, num_triplets / bits, current_process().pid, os.getppid()),
            log_file)
        r('%s = get_valid_triplets(%d, %d, %d)' % (name, num_samples, num_triplets, bits))
        q.put((name, r[name]))
        timer.stop()
        log("\t%s complete (%s)" % (name, str(timer)), log_file)
    except Exception, e:
        q.put("DEATH")
        traceback.print_exc()


def _generate_candiate_free_discrete_matrix(sample_tree, usable_cols):
    assert isinstance(sample_tree, dendropy.Tree)
    print "Creating discrete triplet character matrix"
    r = robjects.r
    treename = _get_random_string(10)
    matrixname = _get_random_string(10)
    newick = sample_tree.as_newick_string()
    num_samples = len(sample_tree.leaf_nodes())
    robjects.globalenv[treename] = newick + ";"
    r("%s = read.tree(text=%s)" % (treename, treename))
    r("%s = get_discrete_matrix(%s, %d)" % (matrixname, treename, usable_cols))
    a = r[matrixname]
    n = r('rownames(%s)' % matrixname)
    return a, n


def _generate_candidate_triplet_discrete_matrix(num_cols, num_samples, sample_tree, bits, usable_cols):
    assert isinstance(sample_tree, dendropy.Tree)
    print "Creating discrete triplet character matrix"
    r = robjects.r
    newick = sample_tree.as_newick_string()
    num_samples = len(sample_tree.leaf_nodes())
    robjects.globalenv['numcols'] = usable_cols
    robjects.globalenv['newick'] = newick + ";"
    r("tree = read.tree(text=newick)")
    r('m = matrix(nrow=length(tree$tip.label))')  #create empty matrix
    r('m = m[,-1]')  #drop the first NA column
    num_procs = mp.cpu_count()
    args = []
    div, mod = divmod(usable_cols, num_procs)
    [args.append(div) for i in range(num_procs)]
    args[-1] += mod
    for i, elem in enumerate(args):
        div, mod = divmod(elem, bits)
        args[-1] += mod
        args[i] -= mod
    manager = Manager()
    pool = Pool(processes=num_procs, maxtasksperchild=1)
    q = manager.Queue(maxsize=num_procs)
    for arg in args:
        pool.apply_async(_get_valid_triplets, (num_samples, arg, bits, q))
    pool.close()
    pool.join()

    while not q.empty():
        name = data = None
        q_data = q.get()
        if q_data == 'DEATH':
            return None, None
        else:
            name, data = q_data
        robjects.globalenv[name] = data
        r('m = cbind(m, %s)' % name)

    r('m = m[,1:%d]' % usable_cols)
    r('m = m[order(rownames(m)),]')  # consistently order the rows (for unifrac compatibility)
    r('m = t(apply(m, 1, as.numeric))')  # convert all factors given by rTraitDisc to numeric
    a = r['m']
    n = r('rownames(m)')
    return a, n


@clockit
def create_discrete_matrix(num_cols, num_samples, sample_tree, bits, accept_cols=False):
    """
    Creates a discrete char matrix from a tree
    @param num_cols: number of columns to create
    @param sample_tree: the tree
    @return: a r object of the matrix, and a list of the row names
    @rtype: tuple(robjects.Matrix, list)
    """
    r = robjects.r
    usable_cols = num_cols
    if not accept_cols:
        usable_cols = find_usable_length(num_cols, bits)
    a, n = _generate_candidate_triplet_discrete_matrix(num_cols, num_samples, sample_tree, bits, usable_cols)
    b, m = _generate_candiate_free_discrete_matrix(sample_tree, usable_cols)

    if a is None:
        log("triplet discrete matrix is none, tryin again", log_file)
        #probably because rpy2 flooded the R pointer stack b/c R sucks.
        return create_discrete_matrix(num_cols, num_samples, sample_tree, bits)

    if b is None:
        log("free discrete matrix is none, tryin again", log_file)
        #probably because rpy2 flooded the R pointer stack b/c R sucks.
        return create_discrete_matrix(num_cols, num_samples, sample_tree, bits)

    assert isinstance(a, robjects.Matrix)
    assert a.ncol == usable_cols

    assert isinstance(b, robjects.Matrix)
    assert b.ncol == usable_cols

    triplet_paralin_matrix, triplet_valid = _create_paralin_matrix(a)
    free_paralin_matrix, free_valid = _create_paralin_matrix(b)

    if triplet_valid is False or free_valid is False:
        sample_tree = create_tree(num_samples, type="S")
        return create_discrete_matrix(num_cols, num_samples, sample_tree, bits)
    else:
        robjects.globalenv['paralin_matrix'] = triplet_paralin_matrix
        r('rownames(paralin_matrix) = rownames(m)')
        r('paralin_dist = as.dist(paralin_matrix, diag=T, upper=T)')
        r("paralinear_cluster = hclust(paralin_dist, method='average')")
    return sample_tree, a, n, b, m


@clockit
def create_tree(num_tips, type):
    """
    creates the taxa tree in R
    @param num_tips: number of taxa to create
    @param type: type for naming (e.g., 'taxa')
    @return: a dendropy Tree
    @rtype: dendropy.Tree
    """
    r = robjects.r
    print "Creating %s tree in %s" % (type, __name__)
    robjects.globalenv['numtips'] = num_tips
    robjects.globalenv['treetype'] = type
    if type == "T":
        r("t = rtree(numtips, rooted=T, tip.label=paste(treetype, seq(1:(numtips)), sep=''))")
    else:
        r("t = rtree(numtips, rooted=F, tip.label=paste(treetype, seq(1:(numtips)), sep=''))")
    tree = r['t']
    return ape_to_dendropy(tree)


def find_usable_length(num_cols, bits):
    """
    Determines usable length given a the number of columns
    in the matrix and the number of bits used in the gap engineering
    @param num_cols: number of columns
    @param bits: number of bits
    @return: max usable length
    @rtype: int
    """
    return max([x for x in range(num_cols + 1) if x % bits ** 2 == 0])


@clockit
def _append_gap_ranges(data, states):
    # pre-populate to avoid finding later, with [-1, -1]
    for range in data:
        for i in xrange(states):
            range[2].append([-1] * 2)

    for range in data:
        d = range[2]
        for i in xrange(int(range[0]), int(range[1]) + 1):
            weight = int(compute_weight(states, i, range[1], range[0]))
            if d[weight][0] == -1:
                d[weight][0] = i
            else:
                d[weight][1] = i


@clockit
def get_range_from_gamma(num_cols, bits, gamma_shape, gamma_scale, smallest_max, accept_cols=False):
    """
    gets a range of random column totals from a gamma distribution
    of given shape and scale
    @param num_cols: number of cols in the discrete matrix
    @param bits: num bits
    @param gamma_shape: gamma shape
    @param gamma_scale: gamma scale
    @return: numpy.ndarray of ranges[[min,max]...]
    @rtype: numpy.ndarray
    """
    data = []
    cols = num_cols
    if not accept_cols:
        cols = find_usable_length(num_cols, bits) / bits
    for x in xrange(0, cols):
        nums = get_random_min_and_max_from_gamma(gamma_shape, gamma_scale, smallest_max)
        data.append(nums)
    states = int("1" * bits, 2) + 1
    _append_gap_ranges(data, states)
    return data


#    return numpy.array(data)

@clockit
def get_range_from_normal(num_cols, bits, mean, sd, smallest_max, accept_cols=False):
    """
    gets a range of random column totals from a normal distribution
    of given shape and scale
    @param num_cols: number of cols in the discrete matrix
    @param bits: num bits
    @param mean: mean
    @param sd: sd
    @param: smallest_max: the smallest value for a max value in an OTU
    @return: numpy.ndarray of ranges[[min,max]...]
    @rtype: numpy.ndarray
    """
    data = []
    if not accept_cols:
        cols = find_usable_length(num_cols, bits) / bits
    for x in xrange(0, cols):
        nums = get_random_min_and_max_from_normal(mean, sd, smallest_max)
        data.append(nums)
    states = int("1" * bits, 2) + 1
    _append_gap_ranges(data, states)
    return numpy.array(data)


def get_random_min_and_max_from_normal(mean, sd, min):
    """
    Gets random min and max from a normal distrubiton
    @param mean: desired mean
    @param sd: desired standard deviation
    @return: tuple of min and max
    @rtype: tuple
    """
    max = get_random_from_normal(mean, sd)
    while max < min:
        max = get_random_from_normal(mean, sd)
    min = round(numpy.random.uniform(low=0.0, high=max / 8.0))
    return [min, max, []]


def get_random_min_and_max_from_gamma(gamma_shape, gamma_scale, min):
    """
    Gets random min and max from a gamma distrubiton
    @param gamma_shape: desired shape
    @param gamma_scale: desired scale
    @return: tuple of min and max
    @rtype: tuple
    """
    max = get_random_from_gamma(gamma_shape, gamma_scale)
    while max < min:
        max = get_random_from_gamma(gamma_shape, gamma_scale)
    min = round(numpy.random.uniform(low=0.0, high=max / 8.0))
    return [min, max, []]


def get_random_from_normal(mean, sd):
    """
    gets a random value from a normal distribution
    @param mean: mean of the normal dist
    @param sd: std dev of the normal dist
    @return: int
    """
    return abs(round(numpy.random.normal(mean, sd)))


def get_random_from_gamma(gamma_shape, gamma_scale):
    """
    gets a random value from a gamma distribution
    @param gamma_shape: shape param of the gamma
    @param gamma_scale: scale param of the gamma
    @return: int
    """
    return abs(round(numpy.random.gamma(gamma_shape, gamma_scale)))


@clockit
def get_range_standardized_matrix_from_discrete(matrix, bits, num_cols):
    """
    Transforms the discrete char matrix into a gap-weighted
    matrix of n bits
    @param matrix: the discrete binary matrix
    @param bits: number of bits to collapse into a single column (OTU)
    @param num_cols: num cols in the matrix
    @return: a new matrix of proper usable length given by bits used
    @rtype: numpy.ndarray
    """
    print "Getting gap-weighted matrix"
    gap = []
    assert isinstance(matrix, rpy2.robjects.vectors.Matrix)
    for rownum in xrange(matrix.nrow):
        row = matrix.rx[rownum + 1, True]
        data = []
        for i in xrange(0, matrix.ncol, bits):
            data.append(int(''.join([str((int(elem))) for elem in row[i:i + bits]]), 2))
        gap.append(data)
    return gap


def compute_weight(num_states, abund, max, min):
    """
    Computes the gap weight
    @param num_states: max char possible
    @param abund: abundance to weight
    @param max: max in column
    @param min: min in column
    @return: the gap weight
    @rtype: int
    """
    w = round(((float(abund) - min) * (num_states - 1)) / (max - min))
    return int(w)


def find_range_limit(weight, abund, max, min, limit, num_states):
    """
    finds the range limit of a gap weight given a starting abundance
    by interating high and low from the starting abundance and finding
    the boundaries that still give the same weight
    @param weight: the weight from the gap matrix
    @param abund: the starting abundance
    @param max: the max of the OTU across samples
    @param min: the min of the OTU across samples
    @param limit: whether we're looking at high or low case
    @param num_states: max char used in gap formula (3 bits = 7)
    @return: the highest minimum abundance or lowest maximum abundance at a given weight
    @rtype: int
    """
    test = abund
    if limit == 'low':
        test -= 1
    elif limit == 'high':
        test += 1

    if test < min:
        return min

    if test > max:
        return max

    test_weight = compute_weight(num_states, test, max, min)

    if test_weight == weight:
        return find_range_limit(weight, test, max, min, limit, num_states)
    else:
        return abund


def find_range(weight, abund, max, min, num_states):
    """
    finds the upper and lower boundary of a weight given and abundance
    @param weight: the gap weight
    @param abund: the starting abundance
    @param max: the max value of an OTU
    @param min: the min value of an OTU
    @param num_states: the max char in the formula (3 bits = 7)
    @return: the upper and lower bound
    @rtype: tuple
    """
    lower = find_range_limit(weight, abund, max, min, "low", num_states)
    upper = find_range_limit(weight, abund, max, min, "high", num_states)
    return lower, upper


def get_random_abundance(weight, col_range):
    """
    gets a random abundance given for an OTU by drawing from
    a uniform distribution on the range interval
    @param weight: the gap weight
    @param abund: the staring aboundance
    @param max: the max value of the OTU
    @param min: the min value of the OTU
    @param num_states: the max char in the formula (3 bits = 7)
    @return: the abundance
    @rypte: int
    """
    r = col_range[2][int(weight)]
    rand = numpy.random.uniform(low=r[0], high=r[1])
    return round(rand)


@clockit
def get_continuous_abundance_matrix(r):
    return numpy.array(r('data_cont')).tolist()


@clockit
def get_abundance_matrix(gap, ranges, dist, num_states):
    """
    computes and returns an abundance matrix given a set
    of column ranges, the type of distribution from which the ranges
    were sampled, and the ranges themselves
    @param gap: the gap matrix (0-7) computed from the discrete character matrix
    @param ranges: the ranges, as a list of tuples [(min, max)...]
    @param dist: the name of the sampling distribution
    @return: a two dimensional list
    @rtype: list
    """
    print "Getting %s abundance matrix" % dist
    data = []
    col_min = [False] * len(gap[0])  # stores whether min has been found for a col
    col_max = [False] * len(gap[0])  # stores whether max has been found for a col

    for i, row in enumerate(gap):
        data.append([None] * len(col_min))
        for j, weight in enumerate(row):
            if weight == 0.0:
                if col_min[j] is False:
                    data[i][j] = ranges[j][0]
                    col_min[j] = True
                else:
                    data[i][j] = get_random_abundance(weight, ranges[j])
            elif weight == (num_states - 1):
                if col_max[j] is False:
                    data[i][j] = ranges[j][1]
                    col_max[j] = True
                else:
                    data[i][j] = get_random_abundance(weight, ranges[j])
            else:
                data[i][j] = get_random_abundance(weight, ranges[j])
    return data


@clockit
def restandardize_matrix(abund, ranges, num_states):
    """
    computes a range-standardized matrix from the computed
    abunance matrix
    @param abund: the abundance matrix
    @return: a two dimensional list
    @raise:
    @rtype: list
    """
    print "Restandardizing abundance matrix"
    data = [None] * len(abund)
    #    ranges = _get_range_for_columns(abund)
    for i, row in enumerate(abund):
        data[i] = [None] * len(row)
        for j, val in enumerate(row):
            range = ranges[j]
            if range[0] == range[1] and range[0] > 0:
                raise Exception("dupe range found at row/col %d/%d %s" % (i, j, range), [row[j] for row in abund])
            if range[1] > 0:
                data[i][j] = compute_weight(num_states, abund=val, max=range[1], min=range[0])
            else:
                data[i][j] = 0
    return data


def _get_range_for_columns(matrix):
    cols = len(matrix[0])
    data = []
    for col in range(cols):
        coldata = [row[col] for row in matrix]
        data.append((min(coldata), max(coldata)))
    return data


def _get_range_for_column(matrix, col):
    return [row[col] for row in matrix]


def get_unifrac_cluster(matrix, rownames):
    name = _get_random_string(20)
    rows = name + "rownames"
    clust = name + "clust"
    phylo = name + "phylo"
    print name, rows, clust, phylo
    r = robjects.r
    assert isinstance(matrix, numpy.ndarray)
    nr, nc = matrix.shape
    matrix_vec = robjects.FloatVector(matrix.transpose().reshape(matrix.size))
    matrix_r = r.matrix(matrix_vec, nrow=nr, ncol=nc)
    robjects.globalenv[name] = matrix_r
    robjects.globalenv[rows] = rownames
    print r[name]
    r("rownames(%s) = %s" % (name, rows))
    print matrix
    print r[name]
    r("%s = hclust(%s, method='average')" % (clust, name))
    r("%s = as.phylo(%s)" % (phylo, clust))
    tree = r('multi2di(%s)' % phylo)
    return ape_to_dendropy(tree)


#def get_unifrac_cluster():
#    """
#    Gets unifrac cluster from the r environment
#    @return: unifrac cluster
#    @rtype: dendropy.Tree
#    """
#    r = robjects.r
#    r("uniclust = hclust(unidist, method='average')")
#    r("uniphylo = as.phylo(uniclust)")
#    tree = r('multi2di(uniphylo)')
#    return ape_to_dendropy(tree)


def print_matrices(abund, ranges, gap, gap2, matrix, matrix2, log_dir, i, dist, num_samples, num_cols):
    """
    Prints the matrices
    @param abund: abundance matrix
    @param ranges: matrix of ranges
    @param gap: gap weighted matrix
    @param matrix: discrete matrix
    @param log_dir: directory for files
    @param i: iteration
    @param dist: name of distribution used to generate ranges
    """
    output_matrix(abund, log_dir, "abund_%d_%d_%s_%d.txt" % (num_samples, num_cols, dist, i), False)
    output_matrix(ranges, log_dir, "ranges_%d_%d_%s_%d.txt" % (num_samples, num_cols, dist, i), False)
    output_matrix(gap, log_dir, "gap_orig_%d_%d_%s_%d.txt" % (num_samples, num_cols, dist, i), False)
    output_matrix(gap2, log_dir, "gap_recon_%d_%d_%s_%d.txt" % (num_samples, num_cols, dist, i), False)
    output_matrix(matrix, log_dir, "matrix_orig_%d_%d_%s_%d.txt" % (num_samples, num_cols, dist, i), True)
    output_matrix(matrix2, log_dir, "matrix_recon_%d_%d_%s_%d.txt" % (num_samples, num_cols, dist, i), True)


def output_matrix(data, folder, file_name, is_r):
    """
    writes a matrix out to a file
    @param data: the matrix to write
    @param folder: the folder
    @param file_name: the filenamer
    @param is_r: whether the matrix is from rpy2
    """
    with open(os.path.join(folder, file_name), 'w') as f:
        if is_r is False:
            for row in data:
                line = '\t'.join([str(int(num)) for num in row])
                f.write(line + "\n")
        else:
            assert isinstance(data, robjects.Matrix)
            for i in xrange(data.nrow):
                line = '\t'.join([str(int(num)) for num in data.rx(i + 1, True)])
                f.write(line + '\n')


def create_mrbayes_file(file, log_file, matrix, sample_names, num_cols, n_gen):
    print "Creating MrBayes file"
    file.write("#NEXUS\n\n")
    file.write("BEGIN DATA;\n")
    file.write("DIMENSIONS NTAX=%d NCHAR=%d;\n" % (matrix.nrow, num_cols))
    file.write("FORMAT DATATYPE=STANDARD;\n")
    file.write("MATRIX\n")
    row_num = 0
    assert isinstance(matrix, robjects.Matrix)
    for i in xrange(matrix.nrow):
        file.write("%s %s\n" % (sample_names[row_num], ''.join([str(int(elem)) for elem in matrix.rx(i + 1, True)])))
        row_num += 1
    file.write(";\n")
    file.write("END;\n\n")

    file.write("begin mrbayes;\n")
    file.write("log start filename=%s replace;\n" % os.path.join(os.path.dirname(file.name), log_file))
    file.write("set autoclose=yes nowarn=yes;\n")
    file.write("lset rates=equal coding=all;\n")
    #file.write("mcmcp checkpoint=yes;\n")
    file.write("mcmcp stoprule=YES stopval=0.01 minpartfreq=0.05;\n")
    file.write("mcmc ngen=%d;\n" % n_gen)
    file.write("sump;\n")
    file.write("sumt;\n")
    file.write("end;\n")
    file.close()


@clockit
def _run_mrbayes_cmd(cmd_string, timeout):
    p = Popen(cmd_string, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    timer = None
    if timeout:
        kill_proc = lambda p: p.kill()
        timer = Timer(timeout, kill_proc, [p])
        timer.start()

    stdout, stderr = p.communicate()
    if p.returncode != 0:
        print "MrBayes timeout, killing on %s!" % platform.uname()[1]
        if timeout:
            timer.cancel()
        return _run_mrbayes_cmd(cmd_string, timeout)

    if timeout:
        timer.cancel()


#    for line in iter(p.stdout.readline, ''):
#        print line.rstrip()
#


@clockit
def run_mrbayes(key, i, matrix, sample_names, num_cols, n_gen, mpi, mb, procs, dist, out_dir, num_samples, name_flag,
                hostfile, timeout):
    """
    Function to run mrbayes and return a tree
    @param i: run iteration
    @param matrix: original binary matrix
    @param sample_names: list of sample names
    @param num_cols: number of columns in the matrix
    @param n_gen: number of mr bayes generations
    @param mpi: path to mpi executable
    @param mb: path to mrbayes executable
    @param procs: number of processors for mrbayesf
    @param dist: the name of the distribution being used to generate totals
    @return: a dendropy Tree object of the mrbayes results
    @rtype: dendropy.Tree
    """
    assert isinstance(matrix, robjects.Matrix)
    print "Running MrBayes"
    dist_name = "n"
    if dist == "gamma":
        dist_name = "g"
    mb_dir = os.path.join(out_dir, "mb")
    if not os.path.exists(mb_dir):
        os.mkdir(mb_dir)
    mb_file = os.path.join(mb_dir,
                           "mb_%d_%d_%s_%s_%s_%s.nex" % (num_samples, matrix.ncol, dist_name, i, name_flag, key))
    log_file = mb_file.replace(".nex", ".log")
    create_mrbayes_file(open(mb_file, "w"), log_file, matrix, sample_names, matrix.ncol, n_gen)
    #cmd = [mpi, "-mca", "pml", "ob1", "-mca", "btl", "self,tcp", "-np", procs, mb, os.path.abspath(mb_file)]
    cmd = [mpi, "-np", procs, mb, os.path.abspath(mb_file)]

    temp_file = None
    if hostfile:
        hosts = []
        for line in open(hostfile):
            hosts.append(line.rstrip())
        random.shuffle(hosts)
        temp_file = tempfile.NamedTemporaryFile(delete=False)
        with temp_file:
            for host in hosts:
                temp_file.write("%s\n" % host)
        cmd = [mpi, "-mca", "pml", "ob1", "-mca", "btl", "self,tcp", "-np", procs, "--hostfile", temp_file.name, mb,
               os.path.abspath(mb_file)]
        #cmd = [mpi, "-np", procs, "--hostfile", temp_file.name, mb, os.path.abspath(mb_file)]

    cmd_string = " ".join([str(elem) for elem in cmd])
    print cmd_string

    _run_mrbayes_cmd(cmd_string, timeout)

    mbresult = os.path.abspath(mb_file) + ".con.tre"

    if not os.path.exists(mbresult):
        mbresult = os.path.abspath(mb_file) + ".con"

    if not os.path.exists(mbresult):
        "can't get %s, restarting run" % mbresult

    tree = dendropy.Tree.get_from_path(mbresult, "nexus")
    assert isinstance(tree, dendropy.Tree)
    tree = make_tree_binary(tree.as_newick_string())
    assert isinstance(tree, dendropy.Tree)
    tree.is_rooted = False
    if temp_file:
        os.unlink(temp_file.name)
    return tree


def calculate_differences_r(orig_tree, test_tree):
    """
    Calcs differences between two trees in R
    @param orig_tree: original dendropy tree of samples
    @param test_tree: dendropy tree of samples to test
    @return: a tuple of distances
    @rtype: tuple
    """

    # makes sure that test_tree actually exists, as can happen
    # when unifrac dist mats cannot do pcoa b/c of all 0
    # eigenvalues
    if test_tree is None:
        return -1, -1, -1

    assert isinstance(orig_tree, dendropy.Tree)
    assert isinstance(test_tree, dendropy.Tree)
    if orig_tree.is_rooted is True:
        orig_tree.is_rooted = False

    if test_tree.is_rooted is True:
        test_tree.is_rooted = False

    r = robjects.r
    robjects.globalenv['origtree'] = orig_tree.as_newick_string() + ";"
    robjects.globalenv['testtree'] = test_tree.as_newick_string() + ";"
    r('o_tree = read.tree(text=origtree)')
    r('t_tree = read.tree(text=testtree)')
    r("topo = dist.topo(o_tree, t_tree, method=\"PH85\")")
    r("symm = phangorn::treedist(o_tree, t_tree)[[1]]")
    r("path =  round(phangorn::treedist(o_tree, t_tree)[[3]], digits=2)")
    return r['topo'][0], r['symm'][0], r['path'][0]


def print_trees_to_pdf(taxa_tree, sample_tree,
                       free_mb_tree, mb_tree, mb_tree2,
                       u_unifrac_tree, w_unifrac_tree,
                       free_mb_diff, mb_diff, mb_diff2,
                       u_uni_diff, w_uni_diff,
                       dist, iter):
    assert isinstance(taxa_tree, dendropy.Tree)
    assert isinstance(sample_tree, dendropy.Tree)
    assert isinstance(mb_tree, dendropy.Tree)
    assert isinstance(mb_tree2, dendropy.Tree)
    assert isinstance(free_mb_tree, dendropy.Tree)
    assert isinstance(u_unifrac_tree, dendropy.Tree)
    assert isinstance(w_unifrac_tree, dendropy.Tree)
    r = robjects.r

    robjects.globalenv['free_mb_diff'] = list(free_mb_diff)
    robjects.globalenv['mb_diff'] = list(mb_diff)
    robjects.globalenv['mb_diff2'] = list(mb_diff2)

    robjects.globalenv['u_uni_diff'] = list(u_uni_diff)
    robjects.globalenv['w_uni_diff'] = list(w_uni_diff)

    robjects.globalenv['taxatree'] = taxa_tree.as_newick_string() + ";"
    robjects.globalenv['sampletree'] = sample_tree.as_newick_string() + ";"
    robjects.globalenv['freeMbResult'] = free_mb_tree.as_newick_string() + ";"
    robjects.globalenv['mbResult'] = mb_tree.as_newick_string() + ";"
    robjects.globalenv['mbResult2'] = mb_tree2.as_newick_string() + ";"

    robjects.globalenv['u_uniclust'] = u_unifrac_tree.as_newick_string() + ";"
    robjects.globalenv['w_uniclust'] = w_unifrac_tree.as_newick_string() + ";"

    robjects.globalenv['distname'] = dist
    robjects.globalenv['iter'] = iter

    r('taxatree = read.tree(text=taxatree)')
    r('sampletree = read.tree(text=sampletree)')
    r('freeMbResult = read.tree(text=freeMbResult)')
    r('mbResult = read.tree(text=mbResult)')
    r('mbResult2 = read.tree(text=mbResult2)')
    r('u_uniclust = read.tree(text=u_uniclust)')
    r('w_uniclust = read.tree(text=w_uniclust)')

    r('plot(sampletree)')
    r("title(main=paste(iter, ': sample tree (', length(sampletree$tip.label), ')', sep=''), "
      "sub=paste(length(taxatree$tip.label), ' taxa', sep=''))")

    r('plot(freeMbResult)')
    r("title(main=paste(distname, ' free mrbayes', sep=''), "
      "sub=paste('diffs=', free_mb_diff[1], '/', free_mb_diff[2], '/', free_mb_diff[3], sep=''))")

    r('plot(mbResult)')
    r("title(main=paste(distname, ' orig mrbayes', sep=''), "
      "sub=paste('diffs=', mb_diff[1], '/', mb_diff[2], '/', mb_diff[3], sep=''))")

    r('plot(mbResult2)')
    r("title(main=paste(distname, ' recon mrbayes', sep=''), "
      "sub=paste('diffs=', mb_diff2[1], '/', mb_diff2[2], '/', mb_diff2[3], sep=''))")

    r("plot(u_uniclust)")
    r("title(main=paste(distname, ' u unifrac', sep=''), "
      "sub=paste('diffs=', u_uni_diff[1], '/', u_uni_diff[2], '/', u_uni_diff[3], sep=''))")

    r("plot(w_uniclust)")
    r("title(main=paste(distname, ' w unifrac', sep=''), "
      "sub=paste('diffs=', w_uni_diff[1], '/', w_uni_diff[2], '/', w_uni_diff[3], sep=''))")

    r('plot.new()')


def get_bc_cluster(abund, sample_names):
    """
    Get average linkage clutering of bray curtis distances derived from
    simulated abunance matrix
    @param abund: abundance matrix
    @param sample_names: array of sample names (rownames)
    @return: a tree representing the cluster
    @rtype: dendropy.Tree
    """
    r = robjects.r
    robjects.globalenv['abund'] = robjects.conversion.py2ri(numpy.array(abund))
    robjects.globalenv['sample_names'] = sample_names
    r('rownames(abund) = sample_names')
    r('bc = vegdist(abund)')
    r("bc_clust = hclust(bc, 'ave')")
    tree = r('multi2di(as.phylo(bc_clust))')
    return ape_to_dendropy(tree)


def get_paralinear_nj():
    """
    Gets an nj tree from paralinear distance matrix
    @return: a dendropy Tree
    @rtype: dendropy.Tree
    """
    r = robjects.r
    r('paralin_nj = nj(paralin_dist)')
    tree = r('multi2di(paralin_nj)')
    return ape_to_dendropy(tree)


@clockit
def get_unifrac_nj(matrix, rownames):
    """
    returns a neighbor joining tree from a unifrac distance matrix
    @param matrix:
    @param rownames:
    @return:
    """
    assert isinstance(matrix, numpy.ndarray)
    r = robjects.r
    name = _get_random_string(20)
    rows = name + "_rows"
    nj = name + "_nj"
    robjects.globalenv[name] = robjects.conversion.py2ri(numpy.array(matrix))
    robjects.globalenv[rows] = robjects.conversion.py2ri(numpy.array(rownames))
    r('rownames(%s) = %s' % (name, rows))
    r('%s = nj(as.dist(%s))' % (nj, name))
    tree = r('multi2di(%s)' % nj)
    return ape_to_dendropy(tree)


def get_bc_nj(abund, sample_names):
    """
    Gets an nj tree from bray curtis distance matrix in r
    @return: a denropy tree
    @rtype: dendropy.Tree
    """
    r = robjects.r
    robjects.globalenv['abund'] = robjects.conversion.py2ri(numpy.array(abund))
    robjects.globalenv['sample_names'] = sample_names
    r('rownames(abund) = sample_names')
    r('bc = vegdist(abund)')
    r('bc_nj = nj(bc)')
    tree = r('multi2di(bc_nj)')
    return ape_to_dendropy(tree)


def get_bc_pcoa_tree(abund, sample_names):
    """
    Gets a tree from the first two pcoa ordination axes of the bray curtis matrix
    The axes are used to create a euclidian distance matrix, this matrix is
    then sent to hclust
    @return: denropy Tree of the results
    @rtype: dendropy.Tree
    """
    r = robjects.r
    robjects.globalenv['abund'] = robjects.conversion.py2ri(numpy.array(abund))
    robjects.globalenv['sample_names'] = sample_names
    r('rownames(abund) = sample_names')
    r('bc = vegdist(abund)')
    r('bc_pcoa = pcoa(bc)')
    r("bc_pcoa_euclid = vegdist(bc_pcoa$vectors[,1:2], 'euc')")
    r("bc_pcoa_clust = hclust(bc_pcoa_euclid, 'ave')")
    tree = r('as.phylo(bc_pcoa_clust)')
    return ape_to_dendropy(tree)


def get_unifrac_pcoa_tree(matrix, rownames):
    r = robjects.r
    assert isinstance(matrix, numpy.ndarray)
    name = _get_random_string(10)
    pcoa = name + "_pcoa"
    euclid = name + "_euclid"
    clust = name + "_clust"
    rows = name + "_rownames"
    robjects.globalenv[name] = robjects.conversion.py2ri(numpy.array(matrix))
    robjects.globalenv[rows] = robjects.conversion.py2ri(numpy.array(rownames))
    r('rownames(%s) = %s' % (name, rows))
    r('%s = as.dist(%s, diag=T, upper=T)' % (name, name))
    try:
        r('%s = pcoa(%s)' % (pcoa, name))
        r("%s = vegdist(%s$vectors[,1:2], 'euc')" % (euclid, pcoa))
        r("%s = hclust(%s, 'ave')" % (clust, euclid))
        tree = r('as.phylo(%s)' % clust)
        return ape_to_dendropy(tree)
    except:
        return None


def get_paralin_pcoa_tree():
    """
    Gets a tree from the first two pcoa ordination axes of the paralinear dist matrix
    The axes are used to create a euclidian distance matrix, this matrix is
    then sent to hclust
    @return: denropy Tree of the results
    @rtype: dendropy.Tree
    """
    r = robjects.r
    r('paralin_pcoa = pcoa(paralin_dist)')
    r("paralin_pcoa_euclid = vegdist(paralin_pcoa$vectors[,1:2], 'euc')")
    r("paralin_pcoa_clust = hclust(paralin_pcoa_euclid, 'ave')")
    tree = r('as.phylo(paralin_pcoa_clust)')
    return ape_to_dendropy(tree)


def ape_to_dendropy(phylo):
    """
    converts an ape tree to dendropy tree
    @param phylo: ape instance from rpy2
    @return: a dendropy tree
    @rtype: dendropy.Tree
    """
    f = tempfile.NamedTemporaryFile()
    robjects.r['write.nexus'](phylo, file=f.name)
    tree = dendropy.Tree.get_from_path(f.name, "nexus")
    f.close()
    return tree


def dendropy_to_cogent(tree):
    """
    converts a dendropy tree to a cogent tree
    @param tree: dendropy tree
    @return a cogent tree
    @rtype: cogent.Tree
    """
    assert isinstance(tree, dendropy.Tree)
    ctree = LoadTree(treestring=tree.as_newick_string())
    return ctree


def close_R():
    """
    Closes down the r session
    """
    robjects.r("dev.off()")


def log(text, writer):
    """
    log out to command line and write easily
    @param text: the text to print
    @param writer: a file object
    """
    assert isinstance(writer, file)
    print text
    writer.write("%s\n" % text)
    writer.flush()


def get_discrete_matrix_from_standardized2(gap, num_states, sample_names):
    """
    creates a binary matrix from a standardized matrix (7 = 11111110) and returns
    that matrix as an rpy2.Matrix
    @param gap: the gap standard matrix
    @param num_states: the nubmer of states to encode the integers
    @param sample_names: the names of the samples (rows in the matrix)
    @return: the rpy2 Matrix
    @rtype: robjects.Matrix
    """
    print "Creating new discrete matrix"
    r = robjects.r
    disc = []
    for row in gap:
        row_string = ''.join([''.join(["1"] * int(elem)).ljust(num_states, "0") for elem in row])
        disc.append([int(elem) for elem in row_string])
    disc = numpy.array(disc)
    r_name = _get_random_string(20)
    robjects.globalenv[r_name] = disc
    robjects.globalenv[r_name + "_rownames"] = sample_names
    r('rownames(%s) = %s' % (r_name, r_name + "_rownames"))
    return r[r_name]


def get_discrete_matrix_from_standardized(gap, bits, sample_names):
    """
    creates a discrete matrix from a standardized matrix (7 = 111) and returns
    that matrix as an rpy2.Matrix
    @param gap: the gap standard matrix
    @param bits: the nubmer of bits to encode the integers
    @param sample_names: the names of the samples (rows in the matrix)
    @return: the rpy2 Matrix
    @rtype: robjects.Matrix
    """
    print "Creating new discrete matrix"
    r = robjects.r
    disc = []
    for row in gap:
        row_string = ''.join([bin(int(elem)).replace("0b", "").rjust(bits, "0") for elem in row])
        disc.append([int(elem) for elem in row_string])
    disc = numpy.array(disc)
    r_name = _get_random_string(20)
    robjects.globalenv[r_name] = disc
    robjects.globalenv[r_name + "_rownames"] = sample_names
    r('rownames(%s) = %s' % (r_name, r_name + "_rownames"))
    return r[r_name]


def _get_random_string(length):
    """
    gets a random string of letters/numbers, ensuring that it does not start with a
    number
    @param length: length of the string
    @return: the random string
    @rtype: string
    """
    s = ''.join(random.choice(string.letters + string.digits) for i in xrange(length))
    if not s[0] in string.letters:
        return _get_random_string(length)
    return s


def correlate_matrices(matrix1, matrix2):
    """
    correlates two matrices
    @param matrix1:
    @param matrix2:
    @return:
    """
    r = robjects.r
    assert isinstance(matrix1, robjects.Matrix)
    assert isinstance(matrix2, robjects.Matrix)
    m1 = []
    m2 = []
    for i in xrange(matrix1.nrow):
        [m1.append(elem) for elem in matrix1.rx(i + 1, True)]
        [m2.append(elem) for elem in matrix2.rx(i + 1, True)]
    return stats.pearsonr(numpy.asarray(m1), numpy.asarray(m2))


@clockit
def subsample_abundance_matrix(matrix, perc):
    """
    :param matrix: abunance matrix as list of lists
    :param perc: percent of items to keep
    """
    sub = []
    prob = perc / 100.0
    for i, row in enumerate(matrix):
        pool = []
        for j, col in enumerate(row):
            pool.extend([j] * int(col))
        pool_sub = numpy.random.choice(pool, size=len(pool) * prob, replace=False)
        hist = _count_unique(pool_sub)
        row_sub = [0] * len(row)
        for index, val in izip(hist[0], hist[1]):
            row_sub[index] = val
        sub.append(row_sub)
    return sub


def _count_unique(items):
    """
    count the unique items in a list
    :param items: a list of items
    :return: a tuple of numpy arrays: (keys, values)
    """
    uniq_keys = numpy.unique(items)
    bins = uniq_keys.searchsorted(items)
    return uniq_keys, numpy.bincount(bins)
