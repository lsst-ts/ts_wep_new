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


class RawExpData(object):
    """Raw exposure data class to populate the raw exposure data."""

    def __init__(self):
        """Construct a raw exposure data class."""

        super().__init__()

        self.visit = []
        self.snap = []
        self.rawExpDir = []

    def append(self, visit, snap, rawExpDir):
        """Append the raw exposure data.

        Parameters
        ----------
        visit : int
            Unique visit Id. This value shoule be >=0.
        snap : int
            Snap (0, 1, etc.). This value shoule be >=0.
        rawExpDir : str
            Raw exposure directory in the local disk.
        """

        if self._isNotNegative(visit):
            self.visit.append(int(visit))

        if self._isNotNegative(snap):
            self.snap.append(int(snap))

        self.rawExpDir.append(rawExpDir)

    def _isNotNegative(self, value):
        """Check the input value is not negative.

        Parameters
        ----------
        value : int or float
            Input value.

        Returns
        -------
        bool
            True if the input value is >= 0.

        Raises
        ------
        ValueError
            The input value should be >= 0.
        """

        isNotNegative = False
        if value >= 0:
            isNotNegative = True
        else:
            raise ValueError("The input value should be >= 0.")

        return isNotNegative

    def reset(self):
        """Reset the internal data."""

        self.__init__()

    def getVisit(self):
        """Get the visit.

        Returns
        -------
        list[int]
            Visit.
        """

        return self.visit

    def getSnap(self):
        """Get the snap.

        Returns
        -------
        list[int]
            Snap.
        """

        return self.snap

    def getRawExpDir(self):
        """Get the raw exposure directory.

        Returns
        -------
        list[str]
            Raw exposure directory.
        """

        return self.rawExpDir


if __name__ == "__main__":
    pass
