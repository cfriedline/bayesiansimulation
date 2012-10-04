import math

__author__ = 'chris'

import unittest
import app

class Tester(unittest.TestCase):
    def setUp(self):
        self.t1 = [1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0]
        self.t2 = [0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0]

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


