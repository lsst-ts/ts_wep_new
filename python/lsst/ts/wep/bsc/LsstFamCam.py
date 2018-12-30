from lsst.obs.lsstSim import LsstSimMapper

from lsst.ts.wep.bsc.CameraData import CameraData


class LsstFamCam(CameraData):

    def __init__(self):
        super(LsstFamCam, self).__init__(LsstSimMapper().camera)

        

if __name__ == "__main__":
    pass
