import os
import sys
import socket
from stopwatch import clockit
import gzip

__author__ = 'chris'

import argparse
import numpy
import multiprocessing as mp
import multiprocessing.pool as mpp
import logging




# From http://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic

class NoDaemonProcess(mp.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False

    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)


class MyPool(mpp.Pool):
    Process = NoDaemonProcess


def percent(abund, total):
    return float(abund) / total


def sub_sample(data, num):
    # print  "\t", mp.current_process().name, mp.current_process().pid, mp.current_process()._parent_pid
    otus = [int(float(i)) for i in data[1:]]
    pool = []
    total = sum(otus)
    for otu, abund in enumerate(otus):
        pool += [otu] * abund

    assert len(pool) == total

    numpy.random.shuffle(pool)
    random_indices = set()
    while len(random_indices) < num:
        random_indices.add(numpy.random.randint(0, len(pool)))
    counts = numpy.zeros(len(otus))
    for i in random_indices:
        otu = pool[i]
        counts[otu] += 1
    return data[0], counts


def sub_sample_wrapper(data):
    return sub_sample(data, 10000)


def sample_abund_file(file):
    print file
    # print file,  mp.current_process().name, mp.current_process().pid, mp.current_process()._parent_pid
    f = gzip.open(file)
    o = gzip.open(os.path.join(os.path.dirname(file), os.path.basename(file) + ".sub.gz"), "wb")

    #if os.path.exists(o.name) and os.path.getsize(o.name) > 0:
    #    return

    pool = mp.Pool(8)
    jobs = []
    for line in f:
        d = line.strip().split("\t")
        jobs.append(d)

    outputs = pool.map(sub_sample_wrapper, jobs)
    for d in outputs:
        o.write("%s\t%s\n" % (d[0], '\t'.join([str(i) for i in d[1]])))
    o.close()
    pool.close()
    pool.join()


def sample_abund_files(abund_files):
    pool = MyPool(10)
    jobs = []
    for file in sorted(abund_files):
        if not 'gap' in file:
            if file.endswith(".gz"):
                jobs.append(file)
    result = pool.map(sample_abund_file, jobs)
    pool.close()
    pool.join()
    print result


def main():
    # logger = mp.log_to_stderr()
    # logger.setLevel(logging.DEBUG)

    print "running on %s\n" % socket.gethostname()
    args = get_options()
    abund_files = get_files(args.dir)
    sample_abund_files(abund_files)


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir')

    if not len(sys.argv) > 1:
        parser.print_help()
        sys.exit(0)

    return parser.parse_args()


@clockit
def get_files(dir):
    out_file = os.path.join(dir, "abund_files.txt")

    if not os.path.exists(out_file):
        with open(out_file, "w") as o:
            for root, dirs, files in os.walk(dir):
                print "root =", root, ", dirs = ", dirs
                if 'mb' in dirs:
                    dirs.remove('mb')

                if 'exclude' in dirs:
                    dirs.remove('exclude')

                for file in files:
                    if 'abund' in file and 'gap' not in file and 'sub' not in file:
                        o.write("%s\n" % os.path.abspath(os.path.join(root, file)))

    files = []
    for line in open(out_file):
        files.append(line.strip())
    print "found %d files" % len(files)
    os.unlink(out_file)
    return files


if __name__ == '__main__':
    main()
    print "Done!"
