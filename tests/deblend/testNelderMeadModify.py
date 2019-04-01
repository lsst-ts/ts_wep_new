import numpy as np
import unittest

from lsst.ts.wep.deblend.nelderMeadModify import feval, nelderMeadModify


class TestNelderMeadModify(unittest.TestCase):
    """Test the nelderMeadModify function."""

    def setUp(self):
        self.func = lambda x, y, c: abs(x**2 + y**2 - c)

    def testFunc(self):
        vars = (1, 2, 1)
        self.assertEqual(feval(self.func, vars), 4)

        xopt = nelderMeadModify(self.func, np.array([2.1]), args=(2, 8,))
        self.assertEqual(xopt[0], 2)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
