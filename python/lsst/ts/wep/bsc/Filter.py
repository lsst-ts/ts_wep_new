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

import os

from lsst.ts.wep.Utility import FilterType, mapFilterRefToG, getConfigDir
from lsst.ts.wep.ParamReader import ParamReader


class Filter(object):
    def __init__(self):
        """Initialize the filter class."""

        # Configuration file of the limit of star's magnitude
        pathMagLimitStar = os.path.join(getConfigDir(), "bsc", "magLimitStar.yaml")
        self._fileMagLimitStar = ParamReader(filePath=pathMagLimitStar)

        # Filter type in use
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

        mappedFilterType = mapFilterRefToG(self.filter)
        if mappedFilterType == FilterType.U:
            filterName = "filterU"

        elif mappedFilterType == FilterType.G:
            filterName = "filterG"

        elif mappedFilterType == FilterType.R:
            filterName = "filterR"

        elif mappedFilterType == FilterType.I:
            filterName = "filterI"

        elif mappedFilterType == FilterType.Z:
            filterName = "filterZ"

        elif mappedFilterType == FilterType.Y:
            filterName = "filterY"

        else:
            raise ValueError("No filter type matches.")

        rangeMag = self._fileMagLimitStar.getSetting(filterName)

        return rangeMag["low"], rangeMag["high"]
