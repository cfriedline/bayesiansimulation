import math

__author__ = 'chris'

import unittest
import app

class Tester(unittest.TestCase):
    def setUp(self):
        self.t1 = [1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0]
        self.t2 = [0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0]
        self.mb = {"321": "/Users/chris/src/mrbayes_3.2.1/src/mb",
                   "322" :"/Users/chris/src/mrbayes_3.2.2/src/mb"}

    def test_paralin(self):
        # Test per Lake (1994)
        assert round(app.get_paralin_distance(self.t1, self.t2), 4) == 1.0617

    def test_weight(self):
        num_states = 8.0
        max = 100
        min = 0
#        assert self._compute_weight(num_states, max, max, min) == 7
#        assert self._compute_weight(num_states, min, max, min) == 0
        d = {}
        range = (max-min) + 1
        print "range", range
        for i in xrange(min, max + 1):
            weight = self._compute_weight(num_states, i, max, min)
            if weight in d:
                d[weight][0] += 1
                d[weight][2] = i
            else:
                d[weight] = [None] * 3
                d[weight][0] = 1
                d[weight][1] = i
        print d

    def _compute_weight(self, num_states, abund, max, min):
        w = round(((abund-min)*(num_states-1))/(max-min))
        return w

    def test_mb_ver(self):
        import os
        import app
        import shutil
        import dendropy
        test_dir = "mb_test"
        if not os.path.exists(test_dir):
            os.mkdir(test_dir)
        procs = 8
        mb_file = "mb_test.nex"
        mpi = "/usr/local/openmpi/bin/mpirun"
        mb_con = {}
        for ver, mb in self.mb.items():
            mb_ver_file = os.path.join(test_dir, "%s_%s" % (ver, mb_file))
            shutil.copy(mb_file, mb_ver_file)
            cmd = [mpi, "-mca", "pml", "ob1", "-mca", "btl", "self,tcp",
               "-np", str(procs),  mb, os.path.abspath(mb_ver_file)]
            app._run_mrbayes_cmd(" ".join(cmd), 600)
            shutil.copy(os.path.join(test_dir, "mb.log"), os.path.join(test_dir, "%s_%s" % (ver, "mb.log")))
            mb_con[ver] = ("%s.con.tre" % mb_ver_file)
        trees = [dendropy.Tree.get_from_path(x, "nexus") for x in mb_con.values()]
        diffs = app.calculate_differences_r(trees[0], trees[1])
        print diffs
        assert diffs[0] == 0.0



