from lsst.ts.wep.ctrlIntf.WEPCalculation import WEPCalculation


class WEPCalculationOfPiston(WEPCalculation):
    """The child class of WEPCalculation that gets the defocal images by the
    camera piston."""

    def __init__(self, astWcsSol, camType, isrDir):
        """Construct an WEP calculation of piston object.

        Parameters
        ----------
        astWcsSol : AstWcsSol
            AST world coordinate system (WCS) solution.
        camType : enum 'CamType'
            Camera type.
        isrDir : str
            Instrument signature remocal (ISR) directory. This directory will
            have the input and output that the data butler needs.
        """
        super(WEPCalculationOfPiston, self).__init__(astWcsSol, camType, isrDir)

    def setDefocalDisInMm(self, defocalDisInMm):
        """Set the defocal distance in mm.

        Parameters
        ----------
        float
            Defocal distance in mm.
        """

        wfEsti = self.wepCntlr.getWfEsti()

        inst = wfEsti.getInst()
        inst.setAnnDefocalDisInMm(defocalDisInMm)

    def getDefocalDisInMm(self):
        """Set the defocal distance in mm.

        CWFS: Curvature wavefront sensor.

        Returns
        -------
        float
            Defocal distance in mm used in the cwfs algorithm.
        """

        wfEsti = self.wepCntlr.getWfEsti()

        inst = wfEsti.getInst()
        defocalDisInM = inst.getDefocalDisOffset()

        return defocalDisInM * 1e3


if __name__ == "__main__":
    pass
