import numpy as np
import unittest

from lsst.ts.wep.DefocalImage import DefocalImage


class TestDefocalImage(unittest.TestCase):
    """Test the defocal image class."""

    def setUp(self):

        self.intraImg = np.arange(2)
        self.extraImg = np.arange(3)

        self.defocalImg = DefocalImage(self.intraImg, self.extraImg)

    def testGetIntraImg(self):

        intraImg = self.defocalImg.getIntraImg()
        self.assertEqual(np.sum(intraImg), np.sum(self.intraImg))

    def testGetExtraImg(self):

        extraImg = self.defocalImg.getExtraImg()
        self.assertEqual(np.sum(extraImg), np.sum(self.extraImg))

    def testSetImg(self):

        intraImg = np.arange(1)
        extraImg = np.arange(2)
        self.defocalImg.setImg(intraImg=intraImg, extraImg=extraImg)

        newIntraImg = self.defocalImg.getIntraImg()
        newExtraImg = self.defocalImg.getExtraImg()
        self.assertEqual(np.sum(newIntraImg), np.sum(intraImg))
        self.assertEqual(np.sum(newExtraImg), np.sum(extraImg))


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
