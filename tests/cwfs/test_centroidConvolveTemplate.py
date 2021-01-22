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

import unittest
import numpy as np

from lsst.ts.wep.cwfs.CentroidConvolveTemplate import CentroidConvolveTemplate


class TestCentroidConvolveTemplate(unittest.TestCase):
    """Test the CentroidConvolveTemplate class."""

    def setUp(self):

        self.centroidConv = CentroidConvolveTemplate()

    def _createData(self):

        radiusInner = 20
        radiusOuter = 40

        # Create two images. One with a single donut and one with two donuts.
        singleDonut = np.zeros((160, 160))
        doubleDonut = np.zeros((160, 160))

        for x in range(160):
            for y in range(160):
                if np.sqrt((80 - x) ** 2 + (80 - y) ** 2) <= radiusOuter:
                    singleDonut[x, y] += 1
                if np.sqrt((80 - x) ** 2 + (80 - y) ** 2) <= radiusInner:
                    singleDonut[x, y] -= 1
                if np.sqrt((50 - x) ** 2 + (80 - y) ** 2) <= radiusOuter:
                    doubleDonut[x, y] += 1
                if np.sqrt((50 - x) ** 2 + (80 - y) ** 2) <= radiusInner:
                    doubleDonut[x, y] -= 1
                if np.sqrt((100 - x) ** 2 + (80 - y) ** 2) <= radiusOuter:
                    doubleDonut[x, y] += 1
                if np.sqrt((100 - x) ** 2 + (80 - y) ** 2) <= radiusInner:
                    doubleDonut[x, y] -= 1
        # Make binary image
        doubleDonut[doubleDonut > 0.5] = 1

        return singleDonut, doubleDonut

    def testGetCenterAndR(self):

        singleDonut, doubleDonut = self._createData()

        # Add noise so the images are not binary
        randState = np.random.RandomState(42)
        doubleDonut += randState.normal(scale=0.1, size=np.shape(doubleDonut))
        singleDonut += randState.normal(scale=0.1, size=np.shape(singleDonut))

        # Use single donut as template and test recovery on double donut
        centX, centY, rad = self.centroidConv.getCenterAndR(
            doubleDonut, singleDonut, 2
        )

        self.assertCountEqual(centX, [50, 100])
        self.assertEqual(centY, [80, 80])
        eff_radius = np.sqrt(40**2 - 20**2)
        self.assertAlmostEqual(rad, eff_radius, delta=0.1)

    def testGetCenterAndRFromImgBinary(self):

        singleDonut, doubleDonut = self._createData()

        # Use single donut as template and test recovery on double donut
        centX, centY, rad = self.centroidConv.getCenterAndRfromImgBinary(
            doubleDonut, singleDonut, 2
        )

        self.assertCountEqual(centX, [50, 100])
        self.assertEqual(centY, [80, 80])
        eff_radius = np.sqrt(40**2 - 20**2)
        self.assertAlmostEqual(rad, eff_radius, delta=0.1)

    def testGetCenterFromTemplateConv(self):

        singleDonut, doubleDonut = self._createData()

        # Test recovery of single donut
        singleCX, singleCY = self.centroidConv.getCenterFromTemplateConv(
            singleDonut, singleDonut, 1
        )
        self.assertEqual(singleCX, [80])
        self.assertEqual(singleCY, [80])

        # Test recovery of two donuts at once
        doubleCX, doubleCY = self.centroidConv.getCenterFromTemplateConv(
            doubleDonut, singleDonut, 2
        )
        self.assertCountEqual(doubleCX, [50, 100])
        self.assertEqual(doubleCY, [80, 80])


if __name__ == "__main__":

    unittest.main()
