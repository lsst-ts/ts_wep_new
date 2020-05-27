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

from lsst.ts.wep.Utility import CamType

from lsst.ts.wep.ctrlIntf.WEPCalculationOfLsstCam import WEPCalculationOfLsstCam
from lsst.ts.wep.ctrlIntf.WEPCalculationOfLsstFamCam import WEPCalculationOfLsstFamCam
from lsst.ts.wep.ctrlIntf.WEPCalculationOfComCam import WEPCalculationOfComCam


class WEPCalculationFactory(object):
    """Factory for creating the correct WEP calculation based off the camera
    type currently being used."""

    def __init__(self):
        """Construct an WEP calculation factory object."""
        super().__init__()

    @staticmethod
    def getCalculator(camType, isrDir):
        """Get a calculator to process wavefront image.

        Parameters
        ----------
        camType : enum 'CamType'
            The camera type to get the wavefront calculator for.
        isrDir : str
            Instrument signature remocal (ISR) directory. This directory will
            have the input and output that the data butler needs.

        Returns
        -------
        WEPCalculationOfLsstCam, WEPCalculationOfLsstFamCam, or
        WEPCalculationOfComCam
            Concrete child class of WEPCalculation class.

        Raises
        ------
        ValueError
            This camera type is not supported.
        """

        if camType == CamType.LsstCam:
            return WEPCalculationOfLsstCam(isrDir)
        elif camType == CamType.LsstFamCam:
            return WEPCalculationOfLsstFamCam(isrDir)
        elif camType == CamType.ComCam:
            return WEPCalculationOfComCam(isrDir)
        else:
            raise ValueError("This camera type is not supported.")


if __name__ == "__main__":
    pass
