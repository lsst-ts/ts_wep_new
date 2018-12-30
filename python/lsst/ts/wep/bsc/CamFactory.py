from lsst.ts.wep.Utility import CamType
from lsst.ts.wep.bsc.LsstCam import LsstCam
from lsst.ts.wep.bsc.LsstFamCam import LsstFamCam
from lsst.ts.wep.bsc.ComCam import ComCam


class CamFactory(object):

    @staticmethod
    def createCam(camType):
        """Create the camera object.

        Parameters
        ----------
        camType : CamType
            Camera type.

        Returns
        -------
        LsstCam, LsstFamCam, or ComCam
            Camera object.
        """

        if (camType == CamType.LsstCam):
            return LsstCam()
        elif (camType == CamType.LsstFamCam):
            return LsstFamCam()
        elif (camType == CamType.ComCam):
            return ComCam()


if __name__ == "__main__":
    pass