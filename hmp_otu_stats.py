__author__ = 'chris'

import sys
import argparse

def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--file")

    if not len(sys.argv) > 1:
        p.print_help()
        exit()

    args = p.parse_args()
    return args

def get_stats(file):
    f = open(file)
    f.readline()
    rows = []
    for line in f:
        d = line.rstrip().split("\t")
        rows.append([int(elem) for elem in d[1:]])

    print "samples =", len(rows)
    print "otus =", len(rows[0])



def main():
    args = get_args()
    get_stats(args.file)

if __name__ == '__main__':
    main()
    print "Done!"


