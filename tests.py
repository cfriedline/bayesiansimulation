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
        max = 1000
        min = 50
        assert app.compute_weight(7, max, max, min) == 7
        assert app.compute_weight(7, min, max, min) == 0
        app.compute_weight(7, max-1, max, min) == 6



