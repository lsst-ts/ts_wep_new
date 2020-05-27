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

from lsst.ts.wep.DefocalImage import DefocalImage


class TestDefocalImage(unittest.TestCase):
    """Test the defocal image class."""

    def setUp(self):

        self.intraImg = np.arange(2)
        self.extraImg = np.arange(3)

        self.defocalImg = DefocalImage(self.intraImg, self.extraImg)

    def testGetIntraImg(self):

        intraImg = self.defocalImg.getIntraImg()
        self.assertEqual(np.sum(intraImg), np.sum(self.intraImg))

    def testGetExtraImg(self):

        extraImg = self.defocalImg.getExtraImg()
        self.assertEqual(np.sum(extraImg), np.sum(self.extraImg))

    def testSetImg(self):

        intraImg = np.arange(1)
        extraImg = np.arange(2)
        self.defocalImg.setImg(intraImg=intraImg, extraImg=extraImg)

        newIntraImg = self.defocalImg.getIntraImg()
        newExtraImg = self.defocalImg.getExtraImg()
        self.assertEqual(np.sum(newIntraImg), np.sum(intraImg))
        self.assertEqual(np.sum(newExtraImg), np.sum(extraImg))


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
