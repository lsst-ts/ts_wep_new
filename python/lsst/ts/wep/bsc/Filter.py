from lsst.ts.wep.Utility import FilterType


class Filter(object):

    def __init__(self):

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
            lowMagnitude = 7.94
            highMagnitude = 14.80

        elif (self.filter == FilterType.G):
            lowMagnitude = 9.74
            highMagnitude = 16.17

        elif (self.filter == FilterType.R):
            lowMagnitude = 9.56
            highMagnitude = 15.73

        elif (self.filter == FilterType.I):
            lowMagnitude = 9.22
            highMagnitude = 15.26

        elif (self.filter == FilterType.Z):
            lowMagnitude = 8.83
            highMagnitude = 14.68
            
        elif (self.filter == FilterType.Y):
            lowMagnitude = 8.02
            highMagnitude = 13.76
        else:
            raise ValueError("No filter type matches.")

        return lowMagnitude, highMagnitude


if __name__ == "__main__":
    pass
