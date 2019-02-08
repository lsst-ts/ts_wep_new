import unittest

from lsst.ts.wep.ctrlIntf.WEPCalculationOfLsstCam import WEPCalculationOfLsstCam


class TestWEPCalculationOfLsstCam(unittest.TestCase):
    """Test the WEPCalculationOfLsstCam class."""

    def setUp(self):
        
        self.wepCalculationOfLsstCam = WEPCalculationOfLsstCam()

    def testIngestRawExp(self):

        self.assertEqual(self.wepCalculationOfLsstCam.visit, 0)
        self.assertEqual(self.wepCalculationOfLsstCam.snap, 0)
        self.assertEqual(self.wepCalculationOfLsstCam.rawExpDir, "")

        visit, snap, rawExpDir = self._ingestRawExp()
        self.assertEqual(self.wepCalculationOfLsstCam.visit, visit)
        self.assertEqual(self.wepCalculationOfLsstCam.snap, snap)
        self.assertEqual(self.wepCalculationOfLsstCam.rawExpDir, rawExpDir)

    def _ingestRawExp(self):

        visit = 99
        snap = 1
        rawExpDir = "temp"
        self.wepCalculationOfLsstCam.ingestRawExp(visit, snap, rawExpDir)

        return visit, snap, rawExpDir

    def testResetRawExpInfo(self):

        self._ingestRawExp()

        self.wepCalculationOfLsstCam._resetRawExpInfo()
        self.assertEqual(self.wepCalculationOfLsstCam.visit, 0)
        self.assertEqual(self.wepCalculationOfLsstCam.snap, 0)
        self.assertEqual(self.wepCalculationOfLsstCam.rawExpDir, "")


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
