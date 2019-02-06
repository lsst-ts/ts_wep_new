class WcsData(object):
    """Contains the world coordinate system (WCS) data of a camera."""

    def __init__(self, wcsCoef):
        """Construct a WCS data class used in AST.

        Parameters
        ----------
        wcsCoef : numpy.ndarray
            WCS coefficients used in the WCS solution.
        """

        self.wcsCoef = wcsCoef

    def setWcsCoef(self, wcsCoef):
        """Set the WCS coefficients.

        Parameters
        ----------
        wcsCoef : numpy.ndarray
            WCS coefficients used in the WCS solution.
        """

        self.wcsCoef = wcsCoef

    def getWcsCoef(self):
        """Get the WCS coefficients.

        Returns
        -------
        numpy.ndarray
            WCS coefficients used in the WCS solution.
        """

        return self.wcsCoef


if __name__ == "__main__":
    pass
