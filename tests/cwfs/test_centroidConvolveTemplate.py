import unittest
import numpy as np

from lsst.ts.wep.cwfs.CentroidConvolveTemplate import CentroidConvolveTemplate


class TestCentroidConvolveTemplate(unittest.TestCase):
    """Test the CentroidConvolveTemplate class."""

    def setUp(self):

        self.centroidConv = CentroidConvolveTemplate()

    def testGetCenterFromTemplateConv(self):

        radiusInner = 20
        radiusOuter = 40

        singleDonut = np.zeros((160, 160))
        doubleDonut = np.zeros((160, 160))

        for x in range(160):
            for y in range(160):
                if np.sqrt((80 - x) ** 2 + (80 - y) ** 2) <= radiusOuter:
                    singleDonut[x, y] += 1
                if np.sqrt((80 - x) ** 2 + (80 - y) ** 2) <= radiusInner:
                    singleDonut[x, y] -= 1
                if np.sqrt((50 - x) ** 2 + (80 - y) ** 2) <= radiusOuter:
                    doubleDonut[x, y] += 1
                if np.sqrt((50 - x) ** 2 + (80 - y) ** 2) <= radiusInner:
                    doubleDonut[x, y] -= 1
                if np.sqrt((100 - x) ** 2 + (80 - y) ** 2) <= radiusOuter:
                    doubleDonut[x, y] += 1
                if np.sqrt((100 - x) ** 2 + (80 - y) ** 2) <= radiusInner:
                    doubleDonut[x, y] -= 1
        # Make binary image
        doubleDonut[doubleDonut > 0.5] = 1

        # Test recovery of single donut
        singleCX, singleCY = self.centroidConv.getCenterFromTemplateConv(
            singleDonut, singleDonut, 1
        )
        self.assertEqual(singleCX, [80])
        self.assertEqual(singleCY, [80])

        # Test recovery of two donuts at once
        doubleCX, doubleCY = self.centroidConv.getCenterFromTemplateConv(
            doubleDonut, singleDonut, 2
        )
        self.assertCountEqual(doubleCX, [50, 100])
        self.assertEqual(doubleCY, [80, 80])


if __name__ == "__main__":

    unittest.main()
