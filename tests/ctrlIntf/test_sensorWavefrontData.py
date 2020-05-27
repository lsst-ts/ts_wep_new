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

from lsst.ts.wep.ctrlIntf.SensorWavefrontData import SensorWavefrontData
from lsst.ts.wep.DonutImage import DonutImage


class TestSensorWavefrontData(unittest.TestCase):
    """Test the SensorWavefrontData class."""

    def setUp(self):

        self.sensorWavefrontData = SensorWavefrontData()

    def testGetMasterDonut(self):

        masterDonut = self.sensorWavefrontData.getMasterDonut()
        self.assertTrue(isinstance(masterDonut, DonutImage))

    def testSetMasterDonut(self):

        masterDonut = self._getDonut()
        self.sensorWavefrontData.setMasterDonut(masterDonut)

        self.assertEqual(self.sensorWavefrontData.getMasterDonut(), masterDonut)

    def _getDonut(self):

        starId = 1
        pixelX = 2.0
        pixelY = 3.0
        fieldX = 1.5
        fieldY = 1.6
        donut = DonutImage(starId, pixelX, pixelY, fieldX, fieldY)

        return donut

    def testGetListOfDonut(self):

        self.assertEqual(self.sensorWavefrontData.getListOfDonut(), [])

        donut = self._getDonut()
        listOfDonut = [donut, donut]
        self.sensorWavefrontData.setListOfDonut(listOfDonut)

        listOfDonutInSensorData = self.sensorWavefrontData.getListOfDonut()
        self.assertEqual(len(listOfDonutInSensorData), len(listOfDonut))


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
