from lsst.obs.lsstSim import LsstSimMapper
from lsst.afw.cameraGeom.detector.detector import DetectorType

from lsst.ts.wep.bsc.CameraData import CameraData


class LsstFamCam(CameraData):

    def __init__(self):
        """Initialize the LSST full-array mode (FAM) camera class."""

        super(LsstFamCam, self).__init__(LsstSimMapper().camera)
        self._initDetectors(DetectorType.SCIENCE)


if __name__ == "__main__":
    pass
