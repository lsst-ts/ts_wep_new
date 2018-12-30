from lsst.obs.lsstSim import LsstSimMapper
from lsst.afw.cameraGeom import SCIENCE

from lsst.ts.wep.bsc.CameraData import CameraData


class ComCam(CameraData):

    def __init__(self):
        """Initialize the commissioning camera class."""

        # The comcam's configuration here is approximated by taking the central 
        # raft of lsst camera.    
        super(ComCam, self).__init__(LsstSimMapper().camera)
        self._initDetectors(SCIENCE)

        # Remove the ccd data that are not belong to ComCam
        detectorList = ["R:2,2 S:0,2", "R:2,2 S:1,2", "R:2,2 S:2,2",
                        "R:2,2 S:0,1", "R:2,2 S:1,1", "R:2,2 S:2,1",
                        "R:2,2 S:0,0", "R:2,2 S:1,0", "R:2,2 S:2,0"]
        self.setWfsCcdList(detectorList)

        wfsCorners = self._rmDictDataNotInList(self.getWfsCorners(),
                                               detectorList)
        self.setWfsCorners(wfsCorners)

        ccdDims = dict()
        for detector in detectorList:
            ccdDims[detector] = self.getCcdDim(detector)
        self.setCcdDims(ccdDims)

    def _rmDictDataNotInList(self, dataDict, keepList):

        # Make a shallow copy
        newDataDict = dict(dataDict)

        # Remove the dictionary data not in the keep list
        for aKey in dataDict.keys():
            if aKey not in keepList:
                newDataDict.pop(aKey)

        return newDataDict


if __name__ == "__main__":
    pass
