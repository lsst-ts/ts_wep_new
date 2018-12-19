import numpy as np
import unittest

from lsst.ts.wep.DefocalImage import DefocalImage


class TestDefocalImage(unittest.TestCase):
    """Test the defocal image class."""

    def setUp(self):
        
        self.defocalImg = DefocalImage(np.arange(2), np.arange(3))

    def testGetIntraImg(self):

        intraImg = self.defocalImg.getIntraImg()
        self.assertEqual(np.sum(intraImg), 1)

    def testGetExtraImg(self):

        extraImg = self.defocalImg.getExtraImg()
        self.assertEqual(np.sum(extraImg), 3)

    def testSetImg(self):
        
        intraImg = np.arange(1)
        extraImg = np.arange(2)
        self.defocalImg.setImg(intraImg=intraImg, extraImg=extraImg)
        
        newIntraImg = self.defocalImg.getIntraImg()
        newExtraImg = self.defocalImg.getExtraImg()
        self.assertEqual(np.sum(newIntraImg), 0)
        self.assertEqual(np.sum(newExtraImg), 1)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
