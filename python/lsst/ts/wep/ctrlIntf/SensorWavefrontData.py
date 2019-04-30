import numpy as np

from lsst.ts.wep.DonutImage import DonutImage


class SensorWavefrontData(object):
    """Sensor wavefront data class that has the information of sensor Id, list
    of donut, master donut, and wavefront error.
    """

    NUM_OF_ZER = 19

    def __init__(self):
        """Construct a sensor wavefront data object."""
        super().__init__()

        self.sensorId = 999
        self.listOfDonut = []
        self.masterDonut = DonutImage(0, 0, 0, 0, 0)
        self.effWfErr = np.zeros(self.NUM_OF_ZER)

    def setSensorId(self, sensorId):
        """Set the sensor Id.

        Parameters
        ----------
        sensorId : int
            The Id of the sensor this wavefront error is for.

        Raises
        ------
        ValueError
            sensorId must be >= 0.
        """

        if (sensorId < 0):
            raise ValueError("sensorId must be >= 0.")
        self.sensorId = int(sensorId)

    def getSensorId(self):
        """Get the sensor Id.

        Returns
        -------
        int
            The Id of the sensor this wavefront error is for.
        """

        return self.sensorId

    def setMasterDonut(self, masterDonut):
        """Set the master donut.

        Parameters
        ----------
        masterDonut : DonutImage
            Master donut.
        """

        self.masterDonut = masterDonut

    def getMasterDonut(self):
        """Get the master donut.

        Returns
        -------
        DonutImage
            Master donut.
        """

        return self.masterDonut

    def setListOfDonut(self, listOfDonut):
        """Set the list of donut on the specific sensor.

        Parameters
        ----------
        listOfDonut : list[DonutImage]
            List of donut.
        """

        self.listOfDonut = listOfDonut

    def getListOfDonut(self):
        """Get the list of donut on the specific sensor.

        Returns
        -------
        list[DonutImage]
            List of donut.
        """

        return self.listOfDonut

    def setAnnularZernikePoly(self, annularZernikePoly):
        """Set the effective annular zernike poly.

        Parameters
        ----------
        annularZernikePoly : numpy.ndarray[19] (float)
            The poly describing the effective wavefront error.

        Raises
        ------
        ValueError
            annularZernikePoly must be an array of 19 floats.
        """

        if (len(annularZernikePoly) != self.NUM_OF_ZER):
            raise ValueError("annularZernikePoly must be an array of %d floats."
                             % self.NUM_OF_ZER)
        self.effWfErr = np.array(annularZernikePoly)

    def getAnnularZernikePoly(self):
        """Get the effective annular zernike poly.

        Returns
        -------
        numpy.ndarray[19] (float)
            The poly describing the effective wavefront error.
        """

        return self.effWfErr


if __name__ == "__main__":
    pass
