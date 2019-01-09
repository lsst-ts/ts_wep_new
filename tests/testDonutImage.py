import numpy as np
import unittest

from lsst.ts.wep.DonutImage import DonutImage


class TestDonutImage(unittest.TestCase):
    """Test the donut image class."""

    def setUp(self):

        self.starId = 0
        self.pixelX = 1
        self.pixelY = 2
        self.fieldX = 3
        self.fieldY = 4
        self.donutImg = DonutImage(self.starId, self.pixelX, self.pixelY,
                                   self.fieldX, self.fieldY)

    def testGetStarId(self):

        self.assertEqual(self.donutImg.getStarId(), self.starId)

    def testGetPixelPos(self):

        pixelX, pixelY = self.donutImg.getPixelPos()

        self.assertEqual(pixelX, self.pixelX)
        self.assertEqual(pixelY, self.pixelY)

    def testGetFieldPos(self):

        fieldX, fieldY = self.donutImg.getFieldPos()

        self.assertEqual(fieldX, self.fieldX)
        self.assertEqual(fieldY, self.fieldY)

    def testGetWfErr(self):

        self.assertEqual(len(self.donutImg.getWfErr()), 0)

    def testSetWfErr(self):

        wfErr = np.arange(19)
        self.donutImg.setWfErr(wfErr)

        recordedWfErr = self.donutImg.getWfErr()
        self.assertEqual(np.sum(np.abs(recordedWfErr-wfErr)), 0)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
