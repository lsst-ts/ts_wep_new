import os
import unittest
import numpy as np

from lsst.ts.wep.Utility import getModulePath
from lsst.ts.wep.deblend.DeblendDefault import DeblendDefault


class TestDeblendDefault(unittest.TestCase):
    """Test the DeblendDefault class."""

    def setUp(self):

        self.deblend = DeblendDefault()

    def testGenerateMultiDonutWithWrongInputs(self):

        img = np.random.rand(4, 4)
        self.assertRaises(ValueError, self.deblend.generateMultiDonut,
                          img, -1.3, 0.3, 45.0)
        self.assertRaises(ValueError, self.deblend.generateMultiDonut,
                          img, 0, 0.3, 45.0)

        self.assertRaises(ValueError, self.deblend.generateMultiDonut,
                          img, 1.3, -0.3, 45.0)
        self.assertRaises(ValueError, self.deblend.generateMultiDonut,
                          img, 1.3, 1.3, 45.0)

    def testGenerateMultiDonut(self):

        template = self._getTemplate()
        image, imageMain, imageNeighbor, neighborX, neighborY = \
            self.deblend.generateMultiDonut(template, 1.3, 0.3, 45.0)

        self.assertEqual(np.rint(np.sum(image)), 7534)
        self.assertAlmostEqual(neighborX, 130.58742823, places=7)
        self.assertAlmostEqual(neighborY, 130.58742823, places=7)

    def _getTemplate(self):

        imageFilePath = os.path.join(getModulePath(), "tests", "testData",
                                     "testImages", "LSST_NE_SN25",
                                     "z11_0.25_intra.txt")
        return np.loadtxt(imageFilePath)

    def testDeblendDonut(self):
        self.assertRaises(NotImplementedError, self.deblend.deblendDonut, [],
                          [])


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
