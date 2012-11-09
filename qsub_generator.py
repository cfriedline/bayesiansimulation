__author__ = 'chris'
import subprocess
import re
import os
import argparse
import sys


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

def write_file(prefix, col, run, args):
    name = None
    with open("%s_%d_%d.sh" % (prefix, col, run), "w") as f:
        name = f.name
        f.write("#!/bin/bash\n")
#        f.write("export PSM_SHAREDCONTEXTS_MAX=8\n")
        f.write("p runsimulation4.py --tree_file 8_taxa.newick_shuffled.txt --cols %d --run %d --brlen %s --project_dir %s --mrbayes_timeout %f\n" % (col, run, args.brlen, args.project_dir, args.mrbayes_timeout))
    return name

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--brlen', default="0.5")
    parser.add_argument('--project_dir')
    parser.add_argument('--mrbayes_timeout', default=600, type=float)

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
        cols = [1000, 3000, 6000, 9000]
        for col in cols:
            for run in xrange(10):
                script = write_file("sim", col, run, args)
                cmd = "qsub -V -N %s-%d-%d -cwd -j y -q %s %s" % (os.path.basename(args.project_dir),col, run, queue_string,  script)
                f.write("%s\n" % cmd)
if __name__ == '__main__':
    main()
    print "Done!"


