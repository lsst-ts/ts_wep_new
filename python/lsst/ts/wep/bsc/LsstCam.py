from lsst.obs.lsstSim import LsstSimMapper

from lsst.ts.wep.bsc.CameraData import CameraData


class LsstCam(CameraData):

    def __init__(self):
        super(LsstCam, self).__init__(LsstSimMapper().camera)

    def getWavefrontSensor(self):
        """
        
        Get the corners of LSST curvature wavefront sensors in (ra, dec) based on the camera_mapper
        list below.
                
        Returns:
            [list] -- (ra, dec) of four corners of each sensor with the name of sensor as a list
        """

        # Camera object
        detectorList = self.getWfsCCdList()
        ra_dec_out = self.getDetectorRaDec(detectorList)

        return ra_dec_out

    def getScineceSensor(self):
        """
        
        Get the corners of LSST science sensors in (ra, dec) based on the camera_mapper list below.
                
        Returns:
            [list] -- (ra, dec) of four corners of each sensor with the name of sensor as a list
        """

        # Camera object
        detectorList = self.getSciCcdList()
        ra_dec_out = self.getDetectorRaDec(detectorList)

        return ra_dec_out


if __name__ == "__main__":
    pass
