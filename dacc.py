import argparse
import sys

__author__ = 'chris'
import numpy
from stopwatch import clockit
import app
import rpy2.robjects as robjects
import os
import datetime

hostname = os.uname()[1]
project_dir = "/Users/chris/projects/bayessim_dacc"
n_gen = 100000
mpi = "/opt/local/bin/openmpirun"
mb = "/Users/chris/src/mrbayes_3.2.1/src/mb"
procs = 4

if 'godel' in hostname:
    mpi = '/test/riveralab/cfriedline/bin/mpirun'
    mb = '/test/riveralab/cfriedline/src/mrbayes_3.2.1/src/mb'
    procs = 8
    project_dir = '/home/cfriedline/projects/bsim_dacc'
elif 'phylogeny' in hostname:
    mpi = '/usr/local/bin/mpirun'
    mb = '/home/cfriedline/src/mrbayes_3.2.1/src/mb'
    procs = 8
    project_dir = '/home/cfriedline/projects/dacc'
    n_gen = 100000000

run_dir_name = datetime.datetime.now().strftime("%m%d%y_%H%M%S")
result_dir = app.create_dir(os.path.join(project_dir, "results"))
run_dir = app.create_dir(os.path.join(result_dir, run_dir_name))

out_dir = app.create_dir(os.path.join(run_dir, "out"))
log_dir = app.create_dir(os.path.join(run_dir, "log"))
log_file = open(os.path.join(log_dir, "log.txt"), "w")
app.log_file = log_file


class Sample:
    def __init__(self, data):
        self.name = data[0]
        self.otus = [int(i) for i in data[1:]]

    def __str__(self):
        return "%s (%s)" % (self.name, len(self.otus))

    def has_body_site(self, b):
        self.has_body_site = b


@clockit
def get_samples(args):
    map_file = args.map_file
    bodysitemap = {}

    f = open("hmp1.v69.hq.otu.counts")
    num = 0
    header = ""
    allsamples = []
    samples_with_bodysite = []
    for line in f:
        if num > 0:
            line = line.strip()
            allsamples.append(Sample(line.split()))
            print "added sample", allsamples[-1], num
        else:
            header = line.split()
        num += 1

    m = open(map_file)
    m.readline()
    for line in m:
        data = line.rstrip().split("\t")
        bodysitemap[data[4]] = data[12]

    for sample in allsamples:
        name = sample.name.split(".")[0]
        if name in bodysitemap:
            samples_with_bodysite.append(sample)

    print "%d samples, %d samples with metadata" % (len(allsamples), len(samples_with_bodysite))

    return samples_with_bodysite, header[1:]


def get_sample_names(samples):
    data = []
    for s in samples:
        data.append(s.name)
    return data


@clockit
def get_otu_data(samples):
    data = []
    for sample in samples:
        data.append(sample.otus)
    return numpy.asarray(data)


@clockit
def get_column_ranges(data):
    ranges = []
    for col in data.T:
        ranges.append((min(col), max(col)))
    return ranges


def print_gap(gap, header, sample_names):
    out = open("dacc_gap.txt", "w")
    with out:
        out.write("%s\t%s\n" % ("collection", '\t'.join(header)))
        for i, row in enumerate(gap):
            out.write("%s\t%s\n" % (sample_names[i], "\t".join([str(int(elem)) for elem in row])))


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--hostfile')
    parser.add_argument('--mrbayes_timeout', type=float)
    parser.add_argument('--map_file', help="file with mappings for samples (used for filtering by body site)")

    args = parser.parse_args()

    if not len(sys.argv) > 1:
        parser.print_help()
        sys.exit()

    return parser.parse_args()


def main():
    args = get_args()
    samples, header = get_samples(args)
    print "num samples = %d, otus = %d" % (len(samples), len(header))
    sample_names = get_sample_names(samples)
    data = get_otu_data(samples)
    ranges = get_column_ranges(data)
    gap = app.restandardize_matrix(data, ranges, num_states=8)
    print_gap(gap, header, sample_names)
    disc = app.get_discrete_matrix_from_standardized(gap, 3, sample_names)
    assert isinstance(disc, robjects.Matrix)
    #run_mrbayes(i, matrix, sample_names, num_cols, n_gen, mpi, mb, procs, dist, out_dir, num_samples, name_flag):
    app.run_mrbayes(0, disc, sample_names, disc.ncol, n_gen, mpi, mb, procs, None, out_dir, len(sample_names), "dacc",
        args.hostfile, args.mrbayes_timeout)

if __name__ == '__main__':
    main()
    print "Done!"


