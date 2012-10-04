__author__ = 'chris'
import os
import pexpect
import multiprocessing as mp

tp = "/Users/chris/bnfo/tree-puzzle/bin/puzzle"

def get_files(dir):
    data = []
    for top, dirs, files in os.walk(dir):
        for file in files:
            if file.endswith(".phy"):
                data.append(os.path.join(top, file))

    return data

def run_tree_puzzle(file):
    cmd = "%s %s" % (tp, file)
    print cmd
    child = pexpect.spawn(cmd)
    child.expect("settings: ")
    child.sendline("b")
    child.expect("settings: ")
    child.sendline("y")
    child.close()

    outfile = file + ".puzzle"
    if not os.path.exists(outfile):
        return False
    return True


def main():
    files = get_files("/Users/chris/simulations")

    pool = mp.Pool()
    results = []
    for file in files:
        pool.apply_async(run_tree_puzzle, (file,), callback=results.append)
    pool.close()
    pool.join()

    assert len(results) == len(files)
    assert False not in results

if __name__ == '__main__':
    main()
    print "Done!"


