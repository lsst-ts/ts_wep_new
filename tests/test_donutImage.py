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

from lsst.ts.wep.DonutImage import DonutImage


class TestDonutImage(unittest.TestCase):
    """Test the donut image class."""

    def setUp(self):

        self.starId = 0
        self.pixelX = 1
        self.pixelY = 2
        self.fieldX = 3
        self.fieldY = 4
        self.donutImg = DonutImage(
            self.starId, self.pixelX, self.pixelY, self.fieldX, self.fieldY
        )

    def testGetStarId(self):

        self.assertEqual(self.donutImg.getStarId(), self.starId)

    def testGetPixelPos(self):

        pixelX, pixelY = self.donutImg.getPixelPos()

        self.assertEqual(pixelX, self.pixelX)
        self.assertEqual(pixelY, self.pixelY)

    def testGetFieldPos(self):

        fieldX, fieldY = self.donutImg.getFieldPos()

        self.assertEqual(fieldX, self.fieldX)
        self.assertEqual(fieldY, self.fieldY)

    def testGetWfErr(self):

        self.assertEqual(len(self.donutImg.getWfErr()), 0)

    def testSetWfErr(self):

        wfErr = np.arange(19)
        self.donutImg.setWfErr(wfErr)

        recordedWfErr = self.donutImg.getWfErr()
        self.assertEqual(np.sum(np.abs(recordedWfErr - wfErr)), 0)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
