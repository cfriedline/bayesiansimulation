__author__ = 'chris'

import sys
import argparse
import re
import subprocess

def get_args():
    p = argparse.ArgumentParser()

#    if not len(sys.argv) > 1:
#        p.print_help()
#        exit()

    args = p.parse_args()
    return args

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

def print_current_time(servers):
    for server in servers:
        proc = subprocess.Popen(['ssh', server, 'date'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for line in proc.stdout:
            line = line.rstrip()
        print server, line

def main():
    args = get_args()
    servers = find_good_servers()
    print_current_time(servers)


if __name__ == '__main__':
    main()
    print "Done!"


