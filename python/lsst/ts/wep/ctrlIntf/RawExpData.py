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

        if (self._isNotNegative(visit)):
            self.visit.append(int(visit))

        if (self._isNotNegative(snap)):
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
        if (value >= 0):
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
