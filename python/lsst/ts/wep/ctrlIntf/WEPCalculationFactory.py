from lsst.ts.wep.Utility import CamType

from lsst.ts.wep.ctrlIntf.WEPCalculationOfLsstCam import WEPCalculationOfLsstCam
from lsst.ts.wep.ctrlIntf.WEPCalculationOfLsstFamCam import \
    WEPCalculationOfLsstFamCam
from lsst.ts.wep.ctrlIntf.WEPCalculationOfComCam import WEPCalculationOfComCam


class WEPCalculationFactory(object):
    """Factory for creating the correct WEP calculation based off the camera
    type currently being used."""

    def __init__(self):
        """Construct an WEP calculation factory object."""
        super().__init__()

    @staticmethod
    def getCalculator(camType, isrDir):
        """Get a calculator to process wavefront image.

        Parameters
        ----------
        camType : CamType
            The camera type to get the wavefront calculator for.
        isrDir : str
            Instrument signature remocal (ISR) directory. This directory will
            have the input and output that the data butler needs.

        Returns
        -------
        WEPCalculationOfLsstCam, WEPCalculationOfLsstFamCam, or
        WEPCalculationOfComCam
            Concrete child class of WEPCalculation class.

        Raises
        ------
        ValueError
            This camera type is not supported.
        """

        if (camType == CamType.LsstCam):
            return WEPCalculationOfLsstCam(isrDir)
        elif (camType == CamType.LsstFamCam):
            return WEPCalculationOfLsstFamCam(isrDir)
        elif (camType == CamType.ComCam):
            return WEPCalculationOfComCam(isrDir)
        else:
            raise ValueError("This camera type is not supported.")


if __name__ == "__main__":
    pass
