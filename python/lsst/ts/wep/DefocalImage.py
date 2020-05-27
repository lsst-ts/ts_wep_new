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


class DefocalImage(object):
    def __init__(self, intraImg=None, extraImg=None):
        """Initialize the DefocalImage class.

        Parameters
        ----------
        intraImg : numpy.ndarray, optional
            Intra-defocal image. (the default is None.)
        extraImg : numpy.ndarray, optional
            Extra-defocal image. (the default is None.)
        """

        # Defocal images
        self.intraImg = intraImg
        self.extraImg = extraImg

    def setImg(self, intraImg=None, extraImg=None):
        """Set the image.

        Parameters
        ----------
        intraImg : numpy.ndarray, optional
            Intra-defocal image. (the default is None.)
        extraImg : numpy.ndarray, optional
            Extra-defocal image. (the default is None.)
        """

        self._setValIfNotNone("intraImg", intraImg)
        self._setValIfNotNone("extraImg", extraImg)

    def _setValIfNotNone(self, attrName, val):
        """Set the value to the related class attribute if the value is not
        none.

        Parameters
        ----------
        attrName : str
            Attribute name.
        val : numpy.ndarray
            Assigned value.
        """

        if val is not None:
            setattr(self, attrName, val)

    def getIntraImg(self):
        """Get the intra-defocal image.

        Returns
        -------
        numpy.ndarray
            Intra-defocal image.
        """

        return self.intraImg

    def getExtraImg(self):
        """Get the extra-defocal image.

        Returns
        -------
        numpy.ndarray
            Extra-defocal image.
        """

        return self.extraImg


if __name__ == "__main__":
    pass
