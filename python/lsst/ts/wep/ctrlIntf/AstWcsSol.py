class AstWcsSol(object):
    """AST world coordinate system (WCS) solution provided by DM team."""

    def __init__(self):
        """Construct a Ast WCS Solution class."""

        super().__init__()
        self.wcsData = None

    def setWcsData(self, wcsData):
        """Set the WCS data.

        Parameters
        ----------
        wcsData : WcsData
            WCS data used in the WCS solution.
        """

        self.wcsData = wcsData


if __name__ == "__main__":
    pass
