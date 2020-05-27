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

from lsst.ts.wep.bsc.LsstCam import LsstCam


class TestLsstCam(unittest.TestCase):
    """Test the LsstCam class."""

    def setUp(self):

        self.cam = LsstCam()
        self.cam.setObsMetaData(0, 0, 0, mjd=59580.0)

    def testGetWfsCcdList(self):

        wfsList = self.cam.getWfsCcdList()
        self.assertEqual(len(wfsList), 8)

    def testGetWavefrontSensor(self):

        wfsData = self.cam.getWavefrontSensor()
        self.assertEqual(len(wfsData), 8)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
