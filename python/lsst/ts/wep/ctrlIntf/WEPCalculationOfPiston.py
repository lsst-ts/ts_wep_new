from lsst.ts.wep.ctrlIntf.WEPCalculation import WEPCalculation


class WEPCalculationOfPiston(WEPCalculation):
    """The child class of WEPCalculation that gets the defocal images by the
    camera piston."""

    DEFOCAL_DIS_IN_MM = 1.5

    def __init__(self, astWcsSol, isrDir):
        """Construct an WEP calculation of piston object.

        Parameters
        ----------
        astWcsSol : AstWcsSol
            AST world coordinate system (WCS) solution.
        isrDir : str
            Instrument signature remocal (ISR) directory. This directory will
            have the input and output that the data butler needs.
        """
        super(WEPCalculationOfPiston, self).__init__(astWcsSol, isrDir)

        self.defocalDisInMm = self.DEFOCAL_DIS_IN_MM

    def setDefocalDisInMm(self, defocalDisInMm):
        """Set the defocal distance in mm.

        Parameters
        ----------
        float
            Defocal distance in mm.
        """

        self.defocalDisInMm = defocalDisInMm

    def getDefocalDisInMm(self):
        """Set the defocal distance in mm.

        Returns
        -------
        float
            Defocal distance in mm.
        """

        return self.defocalDisInMm


if __name__ == "__main__":
    pass
