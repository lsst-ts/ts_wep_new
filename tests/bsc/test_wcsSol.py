# This file is part of ts_wep.
#
# Developed for the LSST Telescope and Site Systems.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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

        xPix = 2047.5
        yPix = 2001.5
        chipName = "R22_S11"
        raByWcs, decByWcs = self.wcs.raDecFromPixelCoords(xPix, yPix, chipName)

        self.assertAlmostEqual(raByWcs, self.ra, places=5)
        self.assertAlmostEqual(decByWcs, self.dec, places=5)

    def testRaDecFromPixelCoordsForChipArray(self):

        xPix = np.array([2047.5, 2047.5])
        yPix = np.array([2001.5, 2001.5])
        chipName = np.array(["R22_S11", "R22_S11"])
        raByWcs, decByWcs = self.wcs.raDecFromPixelCoords(xPix, yPix, chipName)

        self.assertEqual(len(raByWcs), 2)
        self.assertEqual(len(decByWcs), 2)

        self.assertAlmostEqual(raByWcs[0], self.ra, places=5)
        self.assertAlmostEqual(raByWcs[1], self.ra, places=5)

        self.assertAlmostEqual(decByWcs[0], self.dec, places=5)
        self.assertAlmostEqual(decByWcs[1], self.dec, places=5)

    def testPixelCoordsFromRaDecWithoutChipName(self):

        xPix, yPix = self.wcs.pixelCoordsFromRaDec(self.ra, self.dec)

        self.assertAlmostEqual(xPix, 2047, places=-1)
        self.assertAlmostEqual(yPix, 2001, places=-1)

    def testPixelCoordsFromRaDecWithChipName(self):

        chipName = "R22_S11"
        xPix, yPix = self.wcs.pixelCoordsFromRaDec(self.ra, self.dec, chipName=chipName)

        self.assertAlmostEqual(xPix, 2047, places=-1)
        self.assertAlmostEqual(yPix, 2001, places=-1)

    def testFocalPlaneCoordsFromRaDecWithZeroRot(self):

        self.wcs.setObsMetaData(0, 0, 0)
        xInMm, yInMm = self.wcs.focalPlaneCoordsFromRaDec(0, 0)

        self.assertEqual(xInMm, 0.0)
        self.assertEqual(yInMm, 0.0)

        # 0.2 arcsec = 10 um => 1 um = 0.02 arcsec => 1 mm = 20 arcsec
        # 1 arcsec = 1/3600 degree
        xInMm, yInMm = self.wcs.focalPlaneCoordsFromRaDec(20.0 / 3600, 0)

        self.assertAlmostEqual(xInMm, 1.0, places=3)
        self.assertAlmostEqual(yInMm, 0.0, places=3)

    def testFocalPlaneCoordsFromRaDecWithNonZeroRot(self):

        self.wcs.setObsMetaData(0, 0, 45)

        # 0.2 arcsec = 10 um => 1 um = 0.02 arcsec => 1 mm = 20 arcsec
        # 1 arcsec = 1/3600 degree
        xInMm, yInMm = self.wcs.focalPlaneCoordsFromRaDec(20.0 / 3600, 0)

        self.assertAlmostEqual(xInMm, 1 / np.sqrt(2), places=3)
        self.assertAlmostEqual(yInMm, -1 / np.sqrt(2), places=3)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
