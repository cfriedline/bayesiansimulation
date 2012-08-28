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
        max = 10900
        min = 0
#        assert self._compute_weight(num_states, max, max, min) == 7
#        assert self._compute_weight(num_states, min, max, min) == 0
        d = {}
        for i in xrange(max + 1):
            weight = self._compute_weight(num_states, i, max, min)
            if weight == num_states:
#                weight = num_states - 1
                pass
            if weight in d:
                d[weight] += 1
            else:
                d[weight] = 1
        print d

    def _compute_weight(self, num_states, abund, max, min):
        w = (num_states*(abund - min))/(max - min)
#        if w == num_states:
#            w -= 1
        return math.floor(w)


