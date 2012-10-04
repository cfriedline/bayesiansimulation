__author__ = 'chris'
import subprocess
import re
import os


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

def write_file(prefix, col, run):
    name = None
    with open("%s_%d_%d.sh" % (prefix, col, run), "w") as f:
        name = f.name
        f.write("#!/bin/bash\n")
#        f.write("export PSM_SHAREDCONTEXTS_MAX=8\n")
        f.write("p runsimulation4.py --tree_file 8_taxa.newick_shuffled.txt --cols %d --run %d\n" % (col, run))
    return name

def main():
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
                script = write_file("sim", col, run)
                cmd = "qsub -V -N bayessim -cwd -j y -q %s %s" % (queue_string,  script)
                f.write("%s\n" % cmd)
if __name__ == '__main__':
    main()
    print "Done!"


