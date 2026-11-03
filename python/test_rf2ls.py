"""Tests for rf2ls.py: run with  python -m unittest  (from the python/ folder)."""

import unittest
from fractions import Fraction as F

from rf2ls import rf2ls, polyvec


def reconstruct(M, lam):
    return [sum(M[i][j] * lam[j] for j in range(len(lam))) for i in range(len(M))]


class TestRF2LS(unittest.TestCase):
    def check(self, N, factors, expected):
        M, n, labels, lam = rf2ls(N, factors)
        self.assertEqual(len(M), len(M[0]))                 # square
        self.assertEqual(reconstruct(M, lam), n)            # M * lambda = n
        self.assertEqual(lam, [F(x) for x in expected])
        return M

    def test_article_examples(self):
        M = self.check([5, 7], [([1, 1], 1), ([1, 2], 1)], [2, 3])
        self.assertEqual(M, [[1, 1], [2, 1]])
        M = self.check([4, 7], [([1, 0], 1), ([1, 1], 2)], [7, -7, -3])
        self.assertEqual(M, [[1, 1, 0], [2, 1, 1], [1, 0, 0]])
        M = self.check([3, 1], [([1, 1], 1), ([1, 0, 4], 1)], ["-2/5", "2/5", "13/5"])
        self.assertEqual(M, [[1, 1, 0], [0, 1, 1], [4, 0, 1]])

    def test_personalized_item(self):
        a, b, c, d = 22, 33, 98, 76
        self.check([a, b], [([1, c], 2), ([1, d], 1)],
                   [F(a * d - b, (d - c) ** 2), F(b - a * c, d - c), F(b - a * d, (d - c) ** 2)])

    def test_worksheet(self):
        self.check([1, 3], [([1, 1], 1), ([1, 2], 1)], [2, -1])
        self.check([1, 1, 1], [([1, -4], 1), ([1, 2], 1), ([1, 3], 1)], ["1/2", "-1/2", 1])
        self.check([1, "1/2"], [([1, 0], 1), ([1, 2], 1), ([1, 5], 1)], ["1/20", "1/4", "-3/10"])
        self.check([2, 0, 3], [([1, 0], 1), ([1, -1], 2)], [3, -1, 5])
        self.check([1, 0], [([1, 0, 1], 1), ([1, 0, 4], 1)], ["1/3", 0, "-1/3", 0])
        self.check([2, 1, 1], [([1, 1], 2), ([1, 0, 1], 1)], ["-1/2", 1, "1/2", "1/2"])
        self.check([3], [([1, 0], 2), ([1, 0, 1], 1)], [0, 3, 0, -3])

    def test_polyvec(self):
        self.assertEqual(polyvec([1, 1], 3), [0, 1, 1])
        self.assertEqual(polyvec([F(1, 2), 3], 2), [F(1, 2), 3])

    def test_improper_rejected(self):
        with self.assertRaises(ValueError):
            rf2ls([1, 0, 0], [([1, 1], 1), ([1, 2], 1)])


if __name__ == "__main__":
    unittest.main()
