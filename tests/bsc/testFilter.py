import unittest

from lsst.ts.wep.bsc.Filter import Filter


class TestFilter(unittest.TestCase):
    """Test the Filter class."""

    def setUp(self):

        self.filter = "u"

    def testFilter(self):

        afilter = Filter()
        afilter.setFilter(self.filter)

        self.assertEqual(afilter.getFilter(), self.filter)
        self.assertEqual(afilter.getMagBoundary(), (7.94, 14.8))


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
