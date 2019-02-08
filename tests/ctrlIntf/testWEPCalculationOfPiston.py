import unittest

from lsst.ts.wep.ctrlIntf.WEPCalculationOfPiston import WEPCalculationOfPiston
from lsst.ts.wep.ctrlIntf.AstWcsSol import AstWcsSol


class TestWEPCalculationOfPiston(unittest.TestCase):
    """Test the WEPCalculationOfPiston class."""

    def setUp(self):

        self.wepCalculationOfPiston = WEPCalculationOfPiston(AstWcsSol())

    def testIngestIntraRawExp(self):

        self.assertEqual(self.wepCalculationOfPiston.intraVisit, [])
        self.assertEqual(self.wepCalculationOfPiston.intraSnap, [])
        self.assertEqual(self.wepCalculationOfPiston.intraRawExpDir, [])

        visit, snap, rawExpDir = self._ingestIntraRawExp()
        self.assertEqual(self.wepCalculationOfPiston.intraVisit, [visit])
        self.assertEqual(self.wepCalculationOfPiston.intraSnap, [snap])
        self.assertEqual(self.wepCalculationOfPiston.intraRawExpDir,
                         [rawExpDir])

    def _ingestIntraRawExp(self):

        visit = 99
        snap = 1
        rawExpDir = "temp1"
        self.wepCalculationOfPiston.ingestIntraRawExp(visit, snap, rawExpDir)

        return visit, snap, rawExpDir

    def testIngestExtraRawExp(self):

        self.assertEqual(self.wepCalculationOfPiston.extraVisit, [])
        self.assertEqual(self.wepCalculationOfPiston.extraSnap, [])
        self.assertEqual(self.wepCalculationOfPiston.extraRawExpDir, [])

        visit, snap, rawExpDir = self._ingestExtraRawExp()
        self.assertEqual(self.wepCalculationOfPiston.extraVisit, [visit])
        self.assertEqual(self.wepCalculationOfPiston.extraSnap, [snap])
        self.assertEqual(self.wepCalculationOfPiston.extraRawExpDir,
                         [rawExpDir])

    def _ingestExtraRawExp(self):

        visit = 100
        snap = 2
        rawExpDir = "temp2"
        self.wepCalculationOfPiston.ingestExtraRawExp(visit, snap, rawExpDir)

        return visit, snap, rawExpDir

    def testResetRawExpInfo(self):

        self._ingestIntraRawExp()
        self._ingestExtraRawExp()

        self.wepCalculationOfPiston._resetRawExpInfo()
        self.assertEqual(self.wepCalculationOfPiston.intraVisit, [])
        self.assertEqual(self.wepCalculationOfPiston.intraSnap, [])
        self.assertEqual(self.wepCalculationOfPiston.intraRawExpDir, [])
        self.assertEqual(self.wepCalculationOfPiston.extraVisit, [])
        self.assertEqual(self.wepCalculationOfPiston.extraSnap, [])
        self.assertEqual(self.wepCalculationOfPiston.extraRawExpDir, [])

    def testGetDefocalDisInMm(self):

        defocalDisInMm = self.wepCalculationOfPiston.getDefocalDisInMm()
        self.assertEqual(defocalDisInMm,
                         WEPCalculationOfPiston.DEFOCAL_DIS_IN_MM)

    def testSetDefocalDisInMm(self):

        defocalDisInMm = 3.0
        self.wepCalculationOfPiston.setDefocalDisInMm(defocalDisInMm)

        self.assertEqual(self.wepCalculationOfPiston.getDefocalDisInMm(),
                         defocalDisInMm)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
