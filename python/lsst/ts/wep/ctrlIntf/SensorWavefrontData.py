from lsst.ts.wep.ctrlIntf.SensorWavefrontError import SensorWavefrontError
from lsst.ts.wep.DonutImage import DonutImage


class SensorWavefrontData(SensorWavefrontError):
    """Sensor wavefront data class that has the information of sensor Id, list
    of donut, master donut, and wavefront error.
    """

    def __init__(self, numOfZk=19):
        """Construct a sensor wavefront data object.

        Parameters
        ----------
        numOfZk : int, optional
            Number of annular Zernike polynomials. (the default is 19.)
        """

        super(SensorWavefrontData, self).__init__(numOfZk=numOfZk)

        self.listOfDonut = []
        self.masterDonut = DonutImage(0, 0, 0, 0, 0)

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


if __name__ == "__main__":
    pass
