import traceback
import datetime

__author__ = 'chris'

# read results from runsimulation4 to compile results
# that are being computed on many separate nodes

import sys
import argparse
import os
import numpy
from rpy2 import robjects
from rpy2.robjects import numpy2ri
import dendropy
from cStringIO import StringIO

numpy2ri.activate()

def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--dir", help = "root dir of simulation")

    if not len(sys.argv) > 1:
        p.print_help()
        exit()

    args = p.parse_args()
    return args


def get_out_files(dir):
    out_files = []
    out = os.path.join(dir, 'out_files.txt')
    if not os.path.isfile(out):
        with open(out, "w") as f:
            for path, dirs, files in os.walk(dir):
                for file in files:
                    if file == 'out.txt':
                        out_files.append(os.path.join(path, file))
                        #                        print out_files[:-1]
                        f.write(os.path.join(path, file) + "\n")
    else:
        for line in open(out):
            out_files.append(line.rstrip())

    return out_files

def add_to_dict_list(d, key, val):
    if not key in d:
        d[key] = []
    d[key].append(val)

def group_out_files(files):
    group = {}
    for file in files:
        cols = os.path.split(os.path.dirname(os.path.dirname(file)))[1].split("-")[0]
        fh = open(file)
        for line in fh:
            add_to_dict_list(group, cols, line.rstrip())
    return group


def get_column(data, i):
    return [row[i] for row in data]


def process_group(cols, data, group_files):
    f = open(cols + ".txt", "w")
    group_files.append(f.name)
    with f:
        for row in data:
            if not "-1" in row:
                f.write(row + "\n")


def get_file_data(file):
    data = []
    with open(file) as f:
        for line in f:
            if not "tree" in line:
                data.append([int(i) for i in line.rstrip().split("\t")])
    return numpy.array(data)


def get_percent_correct(col_data):
    sum = 0
    for i in col_data:
        if i == 0:
            sum += 1
    return round(float(sum) * 100 / len(col_data), 1)


def get_percent_distance_at(col_data, dist):
    sum = 0
    for i in col_data:
        if i == dist:
            sum += 1
    return round(float(sum) * 100 / len(col_data), 1)


def get_percent_over_(col_data, dist):
    sum = 0
    for i in col_data:
        if i > dist:
            sum += 1
    return round(float(sum) * 100 / len(col_data), 1)

def collapse_missing(missing):
    collapsed = {}
    for k, v in missing.items():
        cols = k.split("/")[6][0:4]
        if not cols in collapsed:
            collapsed[cols] = 0
        collapsed[cols] += v
    return collapsed


def summarize_groups(dir, group_files, missing):
    for k, v in missing.items():
        if v > 0:
            print k, v
    collapsed_missing = collapse_missing(missing)
    dirname = os.path.basename(dir)
    env_key = os.path.basename(dir)
    env = {"asmw4": "brlen=0.5",
           "asmw5":"brlen = U(0.0, 1.0)",
           "asmw6":"brlen=U(0.0, 1.0)",
           "asmw7":"brlen=U(0.1, 1.0)",
           "asmw8":"brlen=0.5",
           "bsim4":"brlen=0.5",
           "bsim5": "brlen=U(0.1, 1.0)",
           "bsim6": "brlen=U(0.1, 1.0)"}
    total_sims = 0
    for file in group_files:
        data = get_file_data(file)
        total_sims += len(data)
        trees = get_column(data, 0)
        print file[:-4], len(set(trees)), "missing =", collapsed_missing[file[:-4]]
        f = open(file)
        header = f.readline().rstrip().split("\t")
        rownames = []
        matrix = []
        for i, elem in enumerate(header):
            if '_symm' in elem:
                row = []
                matrix.append(row)
                rownames.append(elem)
                col_data = get_column(data, i)
#                print col_data
                perc_correct = get_percent_correct(col_data)
                perc_at_2 = get_percent_distance_at(col_data, 2)
                perc_gt_2 = get_percent_over_(col_data, 2)
                print i, elem, len(col_data), perc_correct, perc_at_2, perc_gt_2
                row.append(perc_correct)
                row.append(perc_at_2)
                row.append(perc_gt_2)

        r = robjects.r
        robjects.globalenv['data'] = robjects.conversion.py2ri(numpy.array(matrix))
        robjects.globalenv['rownames'] = robjects.conversion.py2ri(numpy.array(rownames))
        try:
            r('rownames(data) = rownames')
            r("colnames(data) = c('0', '2', '>2')")
            r('pdf("%s-%d_bars.pdf")' % (dirname, int(file[:-4])))
            r('par(mar=c(5.1, 5.1, 5.1, 8.1), xpd=T)')
            r("barcolors=c('#FF0000', '#4169E1', '#FA8072', '#2E8B57', '#6A5ACD', '#708090', '#40E0D0', '#EE82EE', '#F5DEB3', '#FFA500')")
            #r("barplot(as.matrix(data), beside=T, col=barcolors, main='Path difference at %s cols', sub='%d simulations on %d trees', xlab='Difference', ylab='Percent', ylim=c(0, 100))" % (file[:-4], len(col_data), len(set(trees))))
            r("barplot(as.matrix(data), beside=T, col=barcolors, main='Symmetric difference at %s OTUs', xlab='Difference', ylab='Percent', ylim=c(0, 100))" % file[:-4])
            r("mtext('%s')" % env[env_key])
            r("barnames=c('Bayesian', 'PCoA (UU)', 'UPGMA (UU)', 'NJ (UU)', 'PCoA (WU)', 'UPGMA (WU)', 'NJ (WU)', 'PCoA (BC)', 'UPGMA (BC)', 'NJ (BC)')")
            r("legend('topright', inset=c(-0.2, 0), barnames, fill=barcolors, bty='n')")
            r('dev.off()')
        except:
            traceback.print_exc()
            print "Error with %s" % file

    print "total sims = %d (%.2f%%)" % (total_sims, float(total_sims*100)/(10395*4*10))


def compute_diff_for_trees(out_files):
    r = robjects.r
    for file in out_files:
        diffs = []
        mins = []
        maxs = []
        log_dir = os.path.join(os.path.dirname(os.path.dirname(file)), "log")
        f = open(file)
        f.readline()
        for line in f:
            line = line.rstrip().split("\t")
            tree = open(os.path.join(log_dir, "tree_8_%s_%s.txt" % (line[1], line[0])))
            newick = tree.readline().rstrip()
            r("tree = read.tree(text='%s')" % newick)
            edges = r('tree$edge.length')
            diffs.append(line[3])
            mins.append(min(edges))
            maxs.append(max(edges))
        out = os.path.join(os.path.dirname(file), "min_diff.txt")
        print out
        with open(out, "w") as o:
            o.write("min\tmax\tmb_path_diff\n")
            for i in xrange(len(diffs)):
                o.write("%s\t%s\t%s\n" % (mins[i], maxs[i], diffs[i]))
        r("data=read.table('%s', header=T)" % out)
        r("pdf('%s_min.pdf')" % out)
        r("plot(data$min, data$mb_path_diff, xlab='Minimum branch length', ylab='Bayesian symmetric difference', xlim=c(0,1))")
        r("dev.off()")
        r("pdf('%s_max.pdf')" % out)
        r("plot(data$max, data$mb_path_diff, xlab='Maximum branch length', ylab='Bayesian symmetric difference', xlim=c(0,1))")
        r("dev.off()")

def modification_date(filename):
    m = os.path.getmtime(filename)
    return datetime.datetime.fromtimestamp(m)

def compute_mod_time(out_files):
    secs = []
    now = datetime.datetime.now()
    for file in out_files:
        s = (now - modification_date(file)).total_seconds()
        secs.append(s)
    secs.sort()
    return secs


def get_missing(out_files):
    missing = {}
    for file in out_files:
        key = os.path.abspath(file)
        missing[key] = 0
        expected = 0
        total = 10395
        f = open(file)
        f.readline() #skip header
        for line in f:
            l = line.rstrip().split("\t")
            num = int(l[0])
            assert num == expected
            expected += 1
        if expected != total:
            missing[key] = total-expected
    return missing

def print_missing(missing):
    for k, v in missing.items():
        print k, v


def get_bad_mrbayes(out_files, root_dir):
    bad_file = open(os.path.join(root_dir, "bad_mrbayes.txt"), "w")
    with bad_file as bad:
        bad.write("%s\t%s" % ("file", open(out_files[0]).readline()))
        for outfile in out_files:
            f = open(outfile)
            f.readline() # skip header
            for line in f:
                line_data = line.rstrip().split("\t")
                if sum([int(i) for i in line_data[2:5]]) > 0:
                    bad.write("%s\t%s\n" % (outfile, line))
    return bad_file.name


def main():
    robjects.r("library(ape)")
    args = get_args()
    out_files = get_out_files(args.dir)
    bad_mrbayes_file = get_bad_mrbayes(out_files, args.dir)
    missing = get_missing(out_files)
    file_mods = compute_mod_time(out_files)

    out_groups = group_out_files(out_files)
    group_files = []
    for k in out_groups:
        process_group(k, out_groups[k], group_files)
    group_files.sort()
    summarize_groups(args.dir, group_files, missing)

    oldest_mod_time = divmod(file_mods[-1], 60)
    print "oldest file is %d mins and %.2f secs old" % (int(oldest_mod_time[0]), oldest_mod_time[1])


if __name__ == '__main__':
    main()
    print "Done!"


