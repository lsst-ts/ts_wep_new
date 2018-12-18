import unittest

from lsst.ts.wep.bsc.StarData import StarData, NeighboringStar


class TestStarData(unittest.TestCase):
    """Test the StarData and NeighboringStar classes."""

    def setUp(self):
        self.stars = StarData([123, 456, 789], [0.1, 0.2, 0.3],
                              [2.1, 2.2, 2.3], [2.0, 3.0, 4.0], 
                              [2.1, 2.1, 4.1], [2.2, 3.2, 4.2],
                              [2.3, 3.3, 4.3], [2.4, 3.4, 4.4],
                              [2.5, 3.5, 4.5])

    def testStarData(self):
        stars = self.stars
        stars.populateDetector("CCD")

        self.assertEqual(stars.SimobjID, [123, 456, 789])
        self.assertEqual(stars.RA, [0.1, 0.2, 0.3])
        self.assertEqual(stars.Decl, [2.1, 2.2, 2.3])
        self.assertEqual(stars.LSSTMagU, [2.0, 3.0, 4.0])
        self.assertEqual(stars.LSSTMagG, [2.1, 2.1, 4.1])
        self.assertEqual(stars.LSSTMagR, [2.2, 3.2, 4.2])
        self.assertEqual(stars.LSSTMagI, [2.3, 3.3, 4.3])
        self.assertEqual(stars.LSSTMagZ, [2.4, 3.4, 4.4])
        self.assertEqual(stars.LSSTMagY, [2.5, 3.5, 4.5])
        self.assertEqual(stars.Detector,"CCD")

    def testCheckCandidateStars(self):
        stars = self.stars

        indexCandidateU = stars.checkCandidateStars("u", 1.9, 2.1)
        indexCandidateG = stars.checkCandidateStars("g", 0, 5)
        indexCandidateR = stars.checkCandidateStars("r", 0, 1)
        indexCandidateI = stars.checkCandidateStars("i", 2.1, 4.0)
        indexCandidateZ = stars.checkCandidateStars("z", 3.0, 5.0)
        indexCandidateY = stars.checkCandidateStars("y", 1.0, 2.0)

        self.assertEqual(indexCandidateU, [0])
        self.assertEqual(indexCandidateG, [0, 1, 2])
        self.assertEqual(indexCandidateR, [])
        self.assertEqual(indexCandidateI, [0, 1])
        self.assertEqual(indexCandidateZ, [1, 2])
        self.assertEqual(indexCandidateY, [])

    def testGetNeighboringStar(self):
        stars = self.stars
        stars.populateRAData([value*10 for value in stars.RA])
        stars.populateDeclData([value*10 for value in stars.Decl])

        neighboringStarU = stars.getNeighboringStar([0], 3, "u", 99)
        neighboringStarG = stars.getNeighboringStar([0, 1], 3, "g", 99)
        neighboringStarR = stars.getNeighboringStar([0], 1, "r", 99)
        neighboringStarI = stars.getNeighboringStar([], 3, "i", 99)
        neighboringStarZ = stars.getNeighboringStar([0, 1], 2, "z", 1)
        neighboringStarY = stars.getNeighboringStar([1], 2, "y", 1)

        self.assertEqual(stars.RAInPixel, [1, 2, 3])
        self.assertEqual(stars.DeclInPixel, [21, 22, 23])

        self.assertEqual(len(neighboringStarU.SimobjID[123]), 2)
        self.assertEqual(len(neighboringStarU.RaDecl), 3)
        self.assertEqual(len(neighboringStarG.SimobjID), 2)
        self.assertEqual(len(neighboringStarR.SimobjID[123]), 0)
        self.assertEqual(neighboringStarI.SimobjID, {})
        self.assertEqual(len(neighboringStarZ.SimobjID), 1)
        self.assertEqual(neighboringStarY.SimobjID, {})


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
