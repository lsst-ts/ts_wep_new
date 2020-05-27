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
