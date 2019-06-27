import os
import numpy as np
import unittest

from lsst.ts.wep.cwfs.CentroidRandomWalk import CentroidRandomWalk
from lsst.ts.wep.Utility import getModulePath


class TestCentroidRandomWalk(unittest.TestCase):
    """Test the CentroidRandomWalk class."""

    def setUp(self):

        self.centroid = CentroidRandomWalk()

    def testGetCenterAndR(self):

        imgDonut = self._prepareDonutImg(1000)

        realcx, realcy, realR = self.centroid.getCenterAndR(imgDonut)
        self.assertAlmostEqual(realcx, 59.7495, places=3)
        self.assertAlmostEqual(realcy, 59.3421, places=3)
        self.assertAlmostEqual(realR, 47.3616, places=3)

    def _prepareDonutImg(self, seed):

        # Read the image file
        imgFile = os.path.join(getModulePath(), "tests", "testData",
                               "testImages", "LSST_NE_SN25",
                               "z11_0.25_intra.txt")
        imgDonut = np.loadtxt(imgFile)
        # This assumes this "txt" file is in the format
        # I[0,0]   I[0,1]
        # I[1,0]   I[1,1]
        imgDonut = imgDonut[::-1, :]

        # Add the noise to simulate the amplifier image
        np.random.seed(seed=seed)
        d0, d1 = imgDonut.shape
        noise = np.random.rand(d0, d1) * 10

        return imgDonut+noise


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
