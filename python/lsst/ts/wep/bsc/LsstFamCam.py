from lsst.obs.lsstSim import LsstSimMapper
from lsst.afw.cameraGeom import SCIENCE

from lsst.ts.wep.bsc.CameraData import CameraData


class LsstFamCam(CameraData):

    def __init__(self):
        """Initialize the LSST full-array mode (FAM) camera class."""

        super(LsstFamCam, self).__init__(LsstSimMapper().camera)
        self._initDetectors(SCIENCE)


if __name__ == "__main__":
    pass
