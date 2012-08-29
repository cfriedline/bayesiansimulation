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
mpi = "/opt/local/bin/mpirun"
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
    project_dir = '/home/cfriedline/projects/bsim_dacc'
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

@clockit
def get_samples():
    f = open("hmp1.v69.hq.otu.counts")
    num = 0
    header = ""
    samples = []
    for line in f:
        if num > 0:
            line = line.strip()
            samples.append(Sample(line.split()))
            print "added sample", samples[-1], num
        else:
            header = line.split()
        num += 1
    return samples, header[1:]

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

def main():
    samples, header = get_samples()
    print "num samples = %d, otus = %d" % (len(samples), len(header))
    sample_names = get_sample_names(samples)
    data = get_otu_data(samples)
    ranges = get_column_ranges(data)
    gap = app.restandardize_matrix(data, ranges, num_states=8)
    print_gap(gap, header, sample_names)
    disc = app.get_discrete_matrix_from_standardized(gap, 3, sample_names)
    assert isinstance(disc, robjects.Matrix)
    #run_mrbayes(i, matrix, sample_names, num_cols, n_gen, mpi, mb, procs, dist, out_dir, num_samples, name_flag):
    app.run_mrbayes(0, disc, sample_names, disc.ncol, n_gen, mpi, mb, procs, None, out_dir, len(sample_names), "dacc")

if __name__ == '__main__':
    main()
    print "Done!"


