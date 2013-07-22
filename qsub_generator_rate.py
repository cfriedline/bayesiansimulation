__author__ = 'chris'
import subprocess
import re
import os
import argparse
import sys
import numpy

def find_good_servers():
    d = []
    proc = subprocess.Popen(['ssh', 'godel', 'qhost'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in proc.stdout:
        line = line.rstrip()
        if 'godel' in line:
            data = re.split('\s+', line)
            if int(data[2]) > 4:
                d.append(data[0])
    return d

def write_file(prefix, rate, args):
    name = None
    run = 0
    with open("%s_%.2f.sh" % (prefix, rate), "w") as f:
        name = f.name
        f.write("#!/bin/bash\n")
#        f.write("export PSM_SHAREDCONTEXTS_MAX=8\n")
        f.write("p runsimulation4.py --tree_file 8_taxa.newick_shuffled.txt --cols %d --run %d --brlen %s --project_dir %s_%s --mrbayes_timeout %f --rate %.2f\n" % (args.col, run, args.brlen, args.project_dir, str(rate), args.mrbayes_timeout, rate))
    return name

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--brlen', default="0.5")
    parser.add_argument('--project_dir')
    parser.add_argument('--mrbayes_timeout', default=600, type=float)
    parser.add_argument('--col', default=3000, type=int)

    if not len(sys.argv) > 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()

    if not args.project_dir:
        print "Must have project_dir"
        parser.print_help()
        sys.exit()

    return args


def main():
    args = get_args()
    servers = find_good_servers()
    with open("hostfile", "w") as f:
        for server in servers:
            if 'godel97' not in server:
                f.write("%s\n" % server)

    queues = [("all.q@" + server) for server in servers]
    queue_string = ','.join(queues)

    with open("qsubs.sh", "w") as f:
        rates = numpy.arange(0.1, 1.0, 0.1)
        for rate in rates:
            script = write_file("ratesim", rate, args)
            cmd = "qsub -V -N rate%.2f -cwd -j y -q %s %s" % (rate, queue_string, script)
            f.write("%s\n" % cmd)
if __name__ == '__main__':
    main()
    print "Done!"


