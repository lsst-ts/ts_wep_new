import unittest
import numpy as np

from lsst.ts.wep.ctrlIntf.AstWcsSol import AstWcsSol
from lsst.ts.wep.ctrlIntf.WcsData import WcsData


class TestAstWcsSol(unittest.TestCase):
    """Test the AstWcsSol class."""

    def setUp(self):

        self.astWcsSol = AstWcsSol()

    def testSetWcsData(self):

        wcsCoef = np.arange(4)
        wcsData = WcsData(wcsCoef)
        self.astWcsSol.setWcsData(wcsData)

        self.assertTrue(isinstance(self.astWcsSol.wcsData, WcsData))
        self.assertEqual(self.astWcsSol.wcsData, wcsData)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
