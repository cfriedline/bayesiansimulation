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
        max_char = 7
        max = 10000
        min = 0
        assert app.compute_weight(max_char, max, max, min) == 7
        assert app.compute_weight(max_char, min, max, min) == 0
        assert app.compute_weight(max_char, max - 1, max, min) == 6
        assert app.compute_weight(max_char, min + 1, max, min) == 0



