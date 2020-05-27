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

from lsst.ts.wep.DefocalImage import DefocalImage


class DonutImage(DefocalImage):
    def __init__(
        self, starId, pixelX, pixelY, fieldX, fieldY, intraImg=None, extraImg=None
    ):
        """Initialize the donut image class.

        Parameters
        ----------
        starId : int
            Star ID.
        pixelX : float
            Pixel x.
        pixelY : float
            Pixel y.
        fieldX : float
            Field x in degree.
        fieldY : float
            Field y in degree.
        intraImg : numpy.ndarray, optional
            Intra-defocal image. (the default is None.)
        extraImg : numpy.ndarray, optional
            Extra-defocal image. (the default is None.)
        """

        super(DonutImage, self).__init__(intraImg=intraImg, extraImg=extraImg)
        self.starId = int(starId)
        self.pixelX = pixelX
        self.pixelY = pixelY
        self.fieldX = fieldX
        self.fieldY = fieldY

        # Wavefront eror in annular Zk in nm (z4-z22)
        self.zer4UpNm = np.array([])

    def getStarId(self):
        """Get the star Id.

        Returns
        -------
        int
            Star Id.
        """

        return self.starId

    def getPixelPos(self):
        """Get the donut pixel position on sensor.

        Returns
        -------
        float
            Pixel x.
        float
            pixel y.
        """

        return self.pixelX, self.pixelY

    def getFieldPos(self):
        """Get the donut field position in degree.

        Returns
        -------
        float
            Field x in degree.
        float
            Field y in degree.
        """

        return self.fieldX, self.fieldY

    def setWfErr(self, zer4UpNm):
        """Set the wavefront error in annular Zk in nm (z4-z22).

        Parameters
        ----------
        zer4UpNm : numpy.ndarray
            z4 to z22 in nm.
        """

        self.zer4UpNm = zer4UpNm

    def getWfErr(self):
        """Get the wavefront error in annular Zk in nm (z4-z22).

        Returns
        -------
        numpy.ndarray
            z4 to z22 in nm.
        """

        return self.zer4UpNm


if __name__ == "__main__":
    pass
