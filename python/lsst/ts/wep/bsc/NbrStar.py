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

from lsst.ts.wep.Utility import FilterType


class NbrStar(object):
    def __init__(self):
        """Initialize the neighboring star class."""

        self.starId = dict()
        self.raDecl = dict()
        self.raDeclInPixel = dict()

        self.lsstMagU = dict()
        self.lsstMagG = dict()
        self.lsstMagR = dict()
        self.lsstMagI = dict()
        self.lsstMagZ = dict()
        self.lsstMagY = dict()

    def getId(self):
        """Get the star Id.

        Returns
        -------
        dict
            Star Id.
        """

        return self.starId

    def getRaDecl(self):
        """Get the star right ascension (RA) and declination (Decl) in degree.

        Returns
        -------
        dict
            Star RA and Decl in degree.
        """

        return self.raDecl

    def getRaDeclInPixel(self):
        """Get the star right ascension (RA) and declination (Decl) in pixel.

        Returns
        -------
        dict
            Star RA and Decl in pixel.
        """

        return self.raDeclInPixel

    def getMag(self, filterType):
        """Get the star magnitude.

        Parameters
        ----------
        filterType : FilterType
            Filter type.

        Returns
        -------
        dict
            Star magnitude.

        Raises
        ------
        ValueError
            No filter type matches.
        """

        if filterType == FilterType.U:
            return self.lsstMagU
        elif filterType == FilterType.G:
            return self.lsstMagG
        elif filterType == FilterType.R:
            return self.lsstMagR
        elif filterType == FilterType.I:
            return self.lsstMagI
        elif filterType == FilterType.Z:
            return self.lsstMagZ
        elif filterType == FilterType.Y:
            return self.lsstMagY
        else:
            raise ValueError("No filter type matches.")

    def addStar(self, stars, idxCand, idxNbrStar, filterType):
        """Add the information of neightboring stars based on the candidate
        star.

        Parameters
        ----------
        stars : StarData
            Star information.
        idxCand : int
            Index of candidate stars.
        idxNbrStar : numpy.ndarray[int]
            Index of neighboring stars of specific candidate star.
        filterType : str
            Filter type.
        """

        candStarId = stars.getId()
        self.starId[candStarId[idxCand]] = candStarId[idxNbrStar].tolist()

        # Collect coordinates and magnitude of candidate and neighboring stars
        indexStar = np.append(idxNbrStar, idxCand)
        for index in indexStar:
            self.raDecl[candStarId[index]] = (
                stars.getRA()[index],
                stars.getDecl()[index],
            )
            self.raDeclInPixel[candStarId[index]] = (
                stars.getRaInPixel()[index],
                stars.getDeclInPixel()[index],
            )

            starMag = stars.getMag(filterType)[index]
            if filterType == FilterType.U:
                self.lsstMagU[candStarId[index]] = starMag

            elif filterType == FilterType.G:
                self.lsstMagG[candStarId[index]] = starMag

            elif filterType == FilterType.R:
                self.lsstMagR[candStarId[index]] = starMag

            elif filterType == FilterType.I:
                self.lsstMagI[candStarId[index]] = starMag

            elif filterType == FilterType.Z:
                self.lsstMagZ[candStarId[index]] = starMag

            elif filterType == FilterType.Y:
                self.lsstMagY[candStarId[index]] = starMag


if __name__ == "__main__":
    pass
