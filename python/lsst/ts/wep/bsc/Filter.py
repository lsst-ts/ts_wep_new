from lsst.ts.wep.Utility import FilterType


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

        if (self.filter == FilterType.U):
            lowMagnitude = self.U_LOW_MAG
            highMagnitude = self.U_HIGH_MAG

        elif (self.filter == FilterType.G):
            lowMagnitude = self.G_LOW_MAG
            highMagnitude = self.G_HIGH_MAG

        elif (self.filter == FilterType.R):
            lowMagnitude = self.R_LOW_MAG
            highMagnitude = self.R_HIGH_MAG

        elif (self.filter == FilterType.I):
            lowMagnitude = self.I_LOW_MAG
            highMagnitude = self.I_HIGH_MAG

        elif (self.filter == FilterType.Z):
            lowMagnitude = self.Z_LOW_MAG
            highMagnitude = self.Z_HIGH_MAG

        elif (self.filter == FilterType.Y):
            lowMagnitude = self.Y_LOW_MAG
            highMagnitude = self.Y_HIGH_MAG

        else:
            raise ValueError("No filter type matches.")

        return lowMagnitude, highMagnitude


if __name__ == "__main__":
    pass
