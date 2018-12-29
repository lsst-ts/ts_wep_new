import numpy as np
import unittest

from lsst.ts.wep.bsc.WcsSol import WcsSol


class TestWcsSol(unittest.TestCase):
    """Test the ComCam class."""

    def setUp(self):

        self.wcs = WcsSol()

        self.ra = 1.0
        self.dec = 2.0
        rotSkyPos = 0.0
        self.wcs.setObsMetaData(self.ra, self.dec, rotSkyPos)

    def testSetAndGetCamera(self):

        camera = "FaultCamera"
        self.wcs.setCamera(camera)

        self.assertEqual(self.wcs.getCamera(), camera)

    def testSetObservationMetaData(self):

        ra = 10.0
        dec = 20.0
        rotSkyPos = 30.0
        self.wcs.setObsMetaData(ra, dec, rotSkyPos)

        self.assertAlmostEqual(self.wcs._obs.pointingRA, ra)
        self.assertAlmostEqual(self.wcs._obs.pointingDec, dec)
        self.assertAlmostEqual(self.wcs._obs.rotSkyPos, rotSkyPos)

    def testRaDecFromPixelCoordsForSingleChip(self):

        xPix = 2032
        yPix = 2000
        chipName = "R:2,2 S:1,1"
        raByWcs, decByWcs = self.wcs.raDecFromPixelCoords(xPix, yPix, chipName)

        self.assertAlmostEqual(raByWcs, self.ra, places=2)
        self.assertAlmostEqual(decByWcs, self.dec, places=2)

    def testRaDecFromPixelCoordsForChipArray(self):

        xPix = np.array([2032, 2032])
        yPix = np.array([2000, 2000])
        chipName = np.array(["R:2,2 S:1,1", "R:2,2 S:1,1"])
        raByWcs, decByWcs = self.wcs.raDecFromPixelCoords(xPix, yPix, chipName)

        self.assertEqual(len(raByWcs), 2)
        self.assertEqual(len(decByWcs), 2)

        self.assertAlmostEqual(raByWcs[0], self.ra, places=2)
        self.assertAlmostEqual(raByWcs[1], self.ra, places=2)

    def testPixelCoordsFromRaDecWithoutChipName(self):

        xPix, yPix = self.wcs.pixelCoordsFromRaDec(self.ra, self.dec)

        self.assertAlmostEqual(xPix, 2032, places=-1)
        self.assertAlmostEqual(yPix, 1994, places=-1)

    def testPixelCoordsFromRaDecWithChipName(self):

        chipName = "R:2,2 S:1,1"
        xPix, yPix = self.wcs.pixelCoordsFromRaDec(self.ra, self.dec,
                                                   chipName=chipName)

        self.assertAlmostEqual(xPix, 2032, places=-1)
        self.assertAlmostEqual(yPix, 1994, places=-1)

    def testFocalPlaneCoordsFromRaDecWithZeroRot(self):

        self.wcs.setObsMetaData(0, 0, 0)
        xInMm, yInMm = self.wcs.focalPlaneCoordsFromRaDec(0, 0)

        self.assertEqual(xInMm, 0.0)
        self.assertEqual(yInMm, 0.0)

        # 0.2 arcsec = 10 um => 1 um = 0.02 arcsec => 1 mm = 20 arcsec
        # 1 arcsec = 1/3600 degree
        xInMm, yInMm = self.wcs.focalPlaneCoordsFromRaDec(20.0/3600, 0)

        self.assertAlmostEqual(xInMm, 1.0, places=3)
        self.assertAlmostEqual(yInMm, 0.0, places=3)

    def testFocalPlaneCoordsFromRaDecWithNonZeroRot(self):

        self.wcs.setObsMetaData(0, 0, 45)

        # 0.2 arcsec = 10 um => 1 um = 0.02 arcsec => 1 mm = 20 arcsec
        # 1 arcsec = 1/3600 degree
        xInMm, yInMm = self.wcs.focalPlaneCoordsFromRaDec(20.0/3600, 0)

        self.assertAlmostEqual(xInMm, 1/np.sqrt(2), places=3)
        self.assertAlmostEqual(yInMm, -1/np.sqrt(2), places=3)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
