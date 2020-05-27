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

from lsst.ts.wep.ctrlIntf.MapSensorNameAndId import MapSensorNameAndId


class TestMapSensorNameAndId(unittest.TestCase):
    """Test the MapSensorNameAndId class."""

    def setUp(self):

        self.mapping = MapSensorNameAndId()

    def testMapSensorIdToNameWithListInput(self):

        sensorId = [1, 2, 3, 4]
        sensorNameList, numOfsensor = self.mapping.mapSensorIdToName(sensorId)
        self.assertEqual(sensorNameList, ["R00_S21", "R00_S22", "R01_S00", "R01_S01"])
        self.assertEqual(numOfsensor, 4)

        sensorId = [1, -1]
        sensorNameList, numOfsensor = self.mapping.mapSensorIdToName(sensorId)
        self.assertEqual(sensorNameList, ["R00_S21"])
        self.assertEqual(numOfsensor, 1)

        sensorId = []
        sensorNameList, numOfsensor = self.mapping.mapSensorIdToName(sensorId)
        self.assertEqual(sensorNameList, [])
        self.assertEqual(numOfsensor, 0)

    def testMapSensorIdToNameWithIntInput(self):

        sensorId = 1
        sensorNameList, numOfsensor = self.mapping.mapSensorIdToName(sensorId)

        self.assertEqual(sensorNameList, ["R00_S21"])
        self.assertEqual(numOfsensor, 1)

    def testMapSensorNameToIdWithListInput(self):

        sensorName = ["R00_S21", "R00_S22", "R01_S00", "R01_S01"]
        sensorIdList = self.mapping.mapSensorNameToId(sensorName)
        self.assertEqual(sensorIdList, [1, 2, 3, 4])

        sensorName = []
        sensorIdList = self.mapping.mapSensorNameToId(sensorName)
        self.assertEqual(sensorIdList, [])

        incorrectSensorName = ["R00_S21", "R000_S1111"]
        self.assertRaises(
            ValueError, self.mapping.mapSensorNameToId, incorrectSensorName
        )

    def testMapSensorNameToIdWithStrInput(self):

        sensorName = "R00_S21"
        sensorIdList = self.mapping.mapSensorNameToId(sensorName)
        self.assertEqual(sensorIdList, [1])


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
