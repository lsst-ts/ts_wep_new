import unittest
import numpy as np

from lsst.ts.wep.ctrlIntf.SensorWavefrontData import SensorWavefrontData
from lsst.ts.wep.DonutImage import DonutImage


class TestSensorWavefrontData(unittest.TestCase):
    """Test the SensorWavefrontData class."""

    def setUp(self):

        self.sensorWavefrontData = SensorWavefrontData()

    def testGetSensorId(self):

        self.assertEqual(self.sensorWavefrontData.getSensorId(), 999)

    def testSetSensorId(self):

        sensorId = 10
        self.sensorWavefrontData.setSensorId(sensorId)

        self.assertEqual(self.sensorWavefrontData.getSensorId(), sensorId)

    def testSetSensorIdWithWrongId(self):

        sensorId = -2
        self.assertRaises(ValueError, self.sensorWavefrontData.setSensorId,
                          sensorId)

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

    def testGetAnnularZernikePoly(self):

        wfsErr = self.sensorWavefrontData.getAnnularZernikePoly()
        self.assertEqual(len(wfsErr), SensorWavefrontData.NUM_OF_ZER)

    def testSetAnnularZernikePoly(self):

        wfsErr = np.arange(SensorWavefrontData.NUM_OF_ZER)
        self.sensorWavefrontData.setAnnularZernikePoly(wfsErr)

        wfsErrInSensor = self.sensorWavefrontData.getAnnularZernikePoly()

        delta = np.sum(np.abs(wfsErrInSensor - wfsErr))
        self.assertEqual(delta, 0)

    def testSetAnnularZernikePolyWithWrongLeng(self):

        wfsErr = np.arange(SensorWavefrontData.NUM_OF_ZER + 1)
        self.assertRaises(ValueError,
                          self.sensorWavefrontData.setAnnularZernikePoly,
                          wfsErr)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
