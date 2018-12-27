import unittest

from lsst.ts.wep.Utility import FilterType
from lsst.ts.wep.bsc.Filter import Filter


class TestFilter(unittest.TestCase):
    """Test the Filter class."""

    def setUp(self):

        self.filter = Filter()

    def testGetFilter(self):

        self.assertEqual(self.filter.getFilter(), FilterType.U)

    def testSetFilter(self):

        filterType = FilterType.G
        self.filter.setFilter(filterType)
        self.assertEqual(self.filter.getFilter(), filterType)

    def testGetMagBoundary(self):

        self.filter.setFilter(FilterType.G)

        lowMagnitude, highMagnitude = self.filter.getMagBoundary() 
        self.assertEqual(lowMagnitude, 9.74)
        self.assertEqual(highMagnitude, 16.17)

    def testGetMagBoundaryWithError(self):

        self.filter.setFilter("r")
        self.assertRaises(ValueError, self.filter.getMagBoundary)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
