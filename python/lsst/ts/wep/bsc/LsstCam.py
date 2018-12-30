from lsst.obs.lsstSim import LsstSimMapper
from lsst.afw.cameraGeom import WAVEFRONT

from lsst.ts.wep.bsc.CameraData import CameraData


class LsstCam(CameraData):

    def __init__(self):
        """Initialize the LSST camera class."""

        super(LsstCam, self).__init__(LsstSimMapper().camera)
        self._initDetectors(WAVEFRONT)


if __name__ == "__main__":
    pass
