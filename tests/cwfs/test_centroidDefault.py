import unittest
import numpy as np
from scipy.ndimage import generate_binary_structure, iterate_structure

from lsst.ts.wep.cwfs.CentroidDefault import CentroidDefault


class TestCentroidDefault(unittest.TestCase):
    """Test the CentroidDefault class."""

    def setUp(self):

        self.centroid = CentroidDefault()

    def testGetImgBinary(self):

        self.assertRaises(NotImplementedError, self.centroid.getImgBinary,
                          np.zeros((2, 2)))

    def testGetCenterAndRfromImgBinary(self):

        structOri = generate_binary_structure(2, 1).astype(int)

        iterations = 7
        donut = iterate_structure(structOri, iterations)
        dY, dX = donut.shape

        cornerX = 10
        cornerY = 20
        imgBinary = np.zeros((120, 120), dtype=int)
        imgBinary[cornerY:cornerY+dY, cornerX:cornerX+dX] = donut

        x, y, r = self.centroid.getCenterAndRfromImgBinary(imgBinary)
        self.assertEqual(x, cornerX+iterations)
        self.assertEqual(y, cornerY+iterations)
        self.assertAlmostEqual(r, 5.9974, places=3)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
