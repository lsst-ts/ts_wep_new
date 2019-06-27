import unittest

from lsst.ts.wep.Utility import CentroidFindType
from lsst.ts.wep.cwfs.CentroidFindFactory import CentroidFindFactory
from lsst.ts.wep.cwfs.CentroidRandomWalk import CentroidRandomWalk
from lsst.ts.wep.cwfs.CentroidOtsu import CentroidOtsu


class TestCentroidFindFactory(unittest.TestCase):
    """Test the CentroidFindFactory class."""

    def testCreateCentroidFindRandomWalk(self):

        centroidFind = CentroidFindFactory.createCentroidFind(
            CentroidFindType.RandomWalk)
        self.assertTrue(isinstance(centroidFind, CentroidRandomWalk))

    def testCreateCentroidFindOtsu(self):

        centroidFind = CentroidFindFactory.createCentroidFind(
            CentroidFindType.Otsu)
        self.assertTrue(isinstance(centroidFind, CentroidOtsu))

    def testCreateCentroidFindWrongType(self):

        self.assertRaises(ValueError, CentroidFindFactory.createCentroidFind,
                          "wrongType")


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
