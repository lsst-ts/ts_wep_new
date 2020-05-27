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

from lsst.ts.wep.Utility import FilterType, mapFilterRefToG


class Filter(object):

    # Magnitude boundary for each filter type
    U_LOW_MAG = 7.94
    U_HIGH_MAG = 14.80

    G_LOW_MAG = 9.74
    G_HIGH_MAG = 16.17

    R_LOW_MAG = 9.56
    R_HIGH_MAG = 15.73

    I_LOW_MAG = 9.22
    I_HIGH_MAG = 15.26

    Z_LOW_MAG = 8.83
    Z_HIGH_MAG = 14.68

    Y_LOW_MAG = 8.02
    Y_HIGH_MAG = 13.76

    def __init__(self):
        """Initialize the filter class."""

        self.filter = FilterType.U

    def getFilter(self):
        """Get the filter type.

        Returns
        -------
        FilterType
            Filter type.
        """

        return self.filter

    def setFilter(self, filterType):
        """Set the filter type.

        Parameters
        ----------
        filterType : FilterType
            Filter type.
        """

        self.filter = filterType

    def getMagBoundary(self):
        """Get the boundary of magnitude under the current filter type.

        Returns
        -------
        float
            Lower boundary of magnitude.
        float
            Higher boundary of magnitude.

        Raises
        ------
        ValueError
            No filter type matches.
        """

        lowMagnitude = 0
        highMagnitude = 0

        mappedFilterType = mapFilterRefToG(self.filter)
        if mappedFilterType == FilterType.U:
            lowMagnitude = self.U_LOW_MAG
            highMagnitude = self.U_HIGH_MAG

        elif mappedFilterType == FilterType.G:
            lowMagnitude = self.G_LOW_MAG
            highMagnitude = self.G_HIGH_MAG

        elif mappedFilterType == FilterType.R:
            lowMagnitude = self.R_LOW_MAG
            highMagnitude = self.R_HIGH_MAG

        elif mappedFilterType == FilterType.I:
            lowMagnitude = self.I_LOW_MAG
            highMagnitude = self.I_HIGH_MAG

        elif mappedFilterType == FilterType.Z:
            lowMagnitude = self.Z_LOW_MAG
            highMagnitude = self.Z_HIGH_MAG

        elif mappedFilterType == FilterType.Y:
            lowMagnitude = self.Y_LOW_MAG
            highMagnitude = self.Y_HIGH_MAG

        else:
            raise ValueError("No filter type matches.")

        return lowMagnitude, highMagnitude


if __name__ == "__main__":
    pass
