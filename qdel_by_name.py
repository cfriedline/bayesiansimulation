import re
import subprocess

__author__ = 'chris'

import sys
import argparse

def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--name", help = "name of job")

    if not len(sys.argv) > 1:
        p.print_help()
        exit()

    args = p.parse_args()
    return args


def run_qstat(name):
    d = []
    proc = subprocess.Popen(['qstat'], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    for line in proc.stdout:
        line = line.rstrip()
        data = re.split('\s+', line)
        if len(data) > 1 and data[2] == name:
            d.append(data[0])
        else:
            print "skipping %s" % line
    return d


def run_qdel_for_job(job):
    print "killing %s" % job
    proc = subprocess.Popen(['qdel', job], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = proc.communicate()


def run_qdel(jobs):
    for job in jobs:
        run_qdel_for_job(job)


def main():
    args = get_args()
    jobs = run_qstat(args.name)
    run_qdel(jobs)

if __name__ == '__main__':
    main()
    print "Done!"


