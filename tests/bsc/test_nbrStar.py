import numpy as np
import unittest

from lsst.ts.wep.bsc.StarData import StarData
from lsst.ts.wep.bsc.NbrStar import NbrStar
from lsst.ts.wep.Utility import FilterType


class TestNbrStar(unittest.TestCase):
    """Test the NbrStar class."""

    def setUp(self):

        stars = StarData([123, 456, 789], [0.1, 0.2, 0.3],
                         [2.1, 2.2, 2.3], [2.0, 3.0, 4.0],
                         [2.1, 2.1, 4.1], [2.2, 3.2, 4.2],
                         [2.3, 3.3, 4.3], [2.4, 3.4, 4.4],
                         [2.5, 3.5, 4.5])
        stars.setRaInPixel(stars.getRA() * 10)
        stars.setDeclInPixel(stars.getDecl() * 10)

        self.stars = stars
        self.neighboringStar = NbrStar()

    def testGetId(self):

        self._addStar()

        self.assertTrue(123 in self.neighboringStar.getId())
        self.assertTrue(self.neighboringStar.getId()[123], [456])

    def _addStar(self):

        self.neighboringStar.addStar(self.stars, 0, np.array([1]),
                                     FilterType.R)

    def testGetRaDecl(self):

        self._addStar()
        self.assertEqual(self.neighboringStar.getRaDecl(),
                         {456: (0.2, 2.2), 123: (0.1, 2.1)})

    def testGetRaDeclInPixel(self):

        self._addStar()
        self.assertEqual(self.neighboringStar.getRaDeclInPixel(),
                         {456: (2.0, 22.0), 123: (1.0, 21.0)})

    def testGetMag(self):

        self._addStar()
        self.assertEqual(len(self.neighboringStar.getMag(FilterType.R)), 2)
        self.assertEqual(self.neighboringStar.getMag(FilterType.U), {})

    def testAddStarAndGetData(self):

        self._addStar()
        self.assertNotEqual(len(self.neighboringStar.getId()), 0)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
