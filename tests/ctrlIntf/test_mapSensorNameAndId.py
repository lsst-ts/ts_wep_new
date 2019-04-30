import unittest

from lsst.ts.wep.ctrlIntf.MapSensorNameAndId import MapSensorNameAndId


class TestMapSensorNameAndId(unittest.TestCase):
    """Test the MapSensorNameAndId class."""

    def setUp(self):

        self.mapping = MapSensorNameAndId()

    def testMapSensorIdToNameWithListInput(self):

        sensorId = [1, 2, 3, 4]
        sensorNameList, numOfsensor = self.mapping.mapSensorIdToName(
            sensorId)
        self.assertEqual(sensorNameList,
                         ["R00_S21", "R00_S22", "R01_S00", "R01_S01"])
        self.assertEqual(numOfsensor, 4)

        sensorId = [1, -1]
        sensorNameList, numOfsensor = self.mapping.mapSensorIdToName(
            sensorId)
        self.assertEqual(sensorNameList, ["R00_S21"])
        self.assertEqual(numOfsensor, 1)

        sensorId = []
        sensorNameList, numOfsensor = self.mapping.mapSensorIdToName(
            sensorId)
        self.assertEqual(sensorNameList, [])
        self.assertEqual(numOfsensor, 0)

    def testMapSensorIdToNameWithIntInput(self):

        sensorId = 1
        sensorNameList, numOfsensor = self.mapping.mapSensorIdToName(
            sensorId)

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
        self.assertRaises(KeyError, self.mapping.mapSensorNameToId,
                          incorrectSensorName)

    def testMapSensorNameToIdWithStrInput(self):

        sensorName = "R00_S21"
        sensorIdList = self.mapping.mapSensorNameToId(sensorName)
        self.assertEqual(sensorIdList, [1])


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
