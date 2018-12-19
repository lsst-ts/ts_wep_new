import unittest

from lsst.ts.wep.DonutImage import DonutImage


class TestDonutImage(unittest.TestCase):
    """Test the donut image class."""

    def setUp(self):
        
        starId = 0 
        pixelX = 1 
        pixelY = 1 
        fieldX = 2 
        fieldY = 2
        self.donutImg = DonutImage(starId, pixelX, pixelY, fieldX, fieldY)

    def testFunction(self):
        
        self.assertEqual(self.donutImg.pixelX, self.donutImg.pixelY)

        zer4UpNm = 10
        self.donutImg.setWfErr(zer4UpNm)
        self.assertEqual(self.donutImg.zer4UpNm, zer4UpNm)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
