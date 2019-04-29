import unittest

from lsst.ts.wep.Utility import FilterType

from lsst.ts.wep.ctrlIntf.WEPCalculation import WEPCalculation
from lsst.ts.wep.ctrlIntf.AstWcsSol import AstWcsSol
from lsst.ts.wep.ctrlIntf.RawExpData import RawExpData
from lsst.ts.wep.ctrlIntf.SensorWavefrontData import SensorWavefrontData


class TestWEPCalculation(unittest.TestCase):
    """Test the WEPCalculation class."""

    def setUp(self):

        self.isrDir = ""
        self.wepCalculation = WEPCalculation(AstWcsSol(), self.isrDir)

    def testGetIsrDir(self):

        isrDir = self.wepCalculation.getIsrDir()
        self.assertEqual(isrDir, self.isrDir)

    def testGetFilter(self):

        filterType = self.wepCalculation.getFilter()
        self.assertEqual(filterType, FilterType.REF)

    def testSetFilter(self):

        filterType = FilterType.R
        self.wepCalculation.setFilter(filterType)

        self.assertEqual(self.wepCalculation.getFilter(), filterType)

    def testGetBoresight(self):

        ra, dec = self.wepCalculation.getBoresight()

        self.assertEqual(ra, 0)
        self.assertEqual(dec, 0)

    def testSetBoresight(self):

        ra = 1.1
        dec = 2.2
        self.wepCalculation.setBoresight(ra, dec)

        self.assertEqual(self.wepCalculation.getBoresight(), (ra, dec))

    def testGetRotAng(self):

        rotAng = self.wepCalculation.getRotAng()
        self.assertEqual(rotAng, 0.0)

    def testSetRotAng(self):

        rotAng = 10.0
        self.wepCalculation.setRotAng(rotAng)

        self.assertEqual(self.wepCalculation.getRotAng(), rotAng)

    def testSetNumOfProc(self):

        self.assertEqual(self.wepCalculation.numOfProc, 1)

        numOfProc = 3
        self.wepCalculation.setNumOfProc(numOfProc)

        self.assertEqual(self.wepCalculation.numOfProc, numOfProc)

    def testSetNumOfProcWithWrongNum(self):

        numOfProc = 0
        self.assertRaises(ValueError, self.wepCalculation.setNumOfProc,
                          numOfProc)

    def testIngestCalibs(self):

        self.assertEqual(self.wepCalculation.calibsDir, "")

        calibsDir = "temp"
        self.wepCalculation.ingestCalibs(calibsDir)

        self.assertEqual(self.wepCalculation.calibsDir, calibsDir)

    def testCalculateWavefrontErrors(self):

        listOfWfErr = self.wepCalculation.calculateWavefrontErrors(RawExpData())
        self.assertTrue(isinstance(listOfWfErr, list))
        self.assertTrue(isinstance(listOfWfErr[0], SensorWavefrontData))


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
