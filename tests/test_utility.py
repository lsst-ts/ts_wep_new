import os
import unittest

from lsst.ts.wep.Utility import abbrevDectectorName, expandDetectorName, \
    mapFilterRefToG, FilterType, getModulePath, getConfigDir


class TestUtility(unittest.TestCase):
    """Test the Utility functions."""

    def testSerializingAndUnserializingFilterType(self):

        filt = FilterType.fromString('y')
        self.assertEqual(filt, FilterType.Y)
        string = filt.toString()
        self.assertEqual(string, 'y')

    def testAbbrevDectectorName(self):

        sciSensorName = "R:2,2 S:1,1"
        self.assertEqual(abbrevDectectorName(sciSensorName), "R22_S11")

        wfsSensorName = "R:4,0 S:0,2,B"
        self.assertEqual(abbrevDectectorName(wfsSensorName), "R40_S02_C1")

        self.assertRaises(ValueError, abbrevDectectorName, "R:4,0 S:0,2,C")

    def testExpandDetectorName(self):

        sciSensorName = "R22_S11"
        self.assertEqual(expandDetectorName(sciSensorName), "R:2,2 S:1,1")

        wfsSensorName = "R40_S02_C1"
        self.assertEqual(expandDetectorName(wfsSensorName), "R:4,0 S:0,2,B")

        self.assertRaises(ValueError, expandDetectorName, "R40_S02_C2")

    def testmapFilterRefToG(self):

        mappedFilterType = mapFilterRefToG(FilterType.REF)
        self.assertEqual(mappedFilterType, FilterType.G)

    def testmapFilterRefToGForFilterU(self):

        mappedFilterType = mapFilterRefToG(FilterType.U)
        self.assertEqual(mappedFilterType, FilterType.U)

    def testGetConfigDir(self):

        ansConfigDir = os.path.join(getModulePath(), "policy")
        self.assertEqual(getConfigDir(), ansConfigDir)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
