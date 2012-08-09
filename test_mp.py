from time import sleep
import os
import traceback
__author__ = 'chris'

import multiprocessing as mp
from multiprocessing import Pool, Manager

def find_usable_length(num_cols, bits):
    """
    Determines usable length given a the number of columns
    in the matrix and the number of bits used in the gap engineering
    @param num_cols: number of columns
    @param bits: number of bits
    @return: max usable length
    @rtype: int
    """
    return max([x for x in range(num_cols + 1) if x % bits ** 2 == 0])

def __get_valid_triplets(num_samples, arg, bits, q):
    try:
        data = []
        print "\tpid=%d, ppid=%d" % (mp.current_process().pid, os.getppid())
        for i in range(arg):
            data.append(i)

        sleep(5)
        q.put(data)
    except:
        traceback.print_exc()

def do_work(num_samples, num_cols, num_procs, bits):
    args = []
    usable_cols = find_usable_length(num_cols, bits)
    div, mod = divmod(usable_cols, num_procs)
    [args.append(div) for i in range(num_procs)]
    args[-1] += mod
    for i, elem in enumerate(args):
        div, mod = divmod(elem, bits)
        args[-1] += mod
        args[i] -= mod
    manager = Manager()
    pool = Pool(processes=num_procs)
    q = manager.Queue(maxsize=num_procs)
    for arg in args:
        pool.apply_async(__get_valid_triplets, (num_samples, arg, bits, q))
    pool.close()
    pool.join()

    data = []
    while not q.empty():
        data.append(q.get())
    return data

def main():
    num_tests = 1000
    num_samples = 10
    num_cols = 600
    num_procs = mp.cpu_count()
    bits = 3

    for i in xrange(num_tests):
        print "At", i
        do_work(num_samples, num_cols, num_procs, bits)

if __name__ == '__main__':
    main()


