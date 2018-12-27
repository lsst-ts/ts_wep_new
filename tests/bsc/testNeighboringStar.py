import unittest

from lsst.ts.wep.bsc.StarData import StarData
from lsst.ts.wep.bsc.NeighboringStar import NeighboringStar


class TestNeighboringStar(unittest.TestCase):
    """Test the NeighboringStar class."""
    
    def setUp(self):

        stars = StarData([123, 456, 789], [0.1, 0.2, 0.3],
                         [2.1, 2.2, 2.3], [2.0, 3.0, 4.0], 
                         [2.1, 2.1, 4.1], [2.2, 3.2, 4.2],
                         [2.3, 3.3, 4.3], [2.4, 3.4, 4.4],
                         [2.5, 3.5, 4.5])
        stars.populateRAData([value*10 for value in stars.RA])
        stars.populateDeclData([value*10 for value in stars.Decl])

        self.stars = stars
        self.neighboringStar = NeighboringStar()

    def testAddStar(self):

        self.neighboringStar.addStar(self.stars, 0, 1, "r")

        self.assertTrue(123 in self.neighboringStar.SimobjID)
        self.assertTrue(self.neighboringStar.SimobjID[123], [456])

        self.assertEqual(self.neighboringStar.RaDecl,
                         {456: (0.2, 2.2), 123: (0.1, 2.1)})
        self.assertEqual(self.neighboringStar.RaDeclInPixel, 
                         {456: (2.0, 22.0), 123: (1.0, 21.0)})

        self.assertEqual(len(self.neighboringStar.LSSTMagR), 2)
        self.assertEqual(self.neighboringStar.LSSTMagU, {})


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
