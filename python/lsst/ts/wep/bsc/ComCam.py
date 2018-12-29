from lsst.obs.lsstSim import LsstSimMapper

from lsst.ts.wep.bsc.CameraData import CameraData


class ComCam(CameraData):

    def __init__(self):
        # The comcam's configuration here is approximated by taking the central 
        # raft of lsst camera.    
        super(ComCam, self).__init__(LsstSimMapper().camera)

    def getWavefrontSensor(self):
        """
        
        Get the corner sensors of Comcam in (ra, dec) based on the camera_mapper list below.
        The reference is at:
        https://confluence.lsstcorp.org/display/LSWUG/Representation+of+a+Camera
        
        Arguments:
            obs {[metadata]} -- Instantiation of ObservationMetaData that describes the pointing
                                of the telescope.
        
        Returns:
            [list] -- (ra, dec) of four corners of each sensor with the name of sensor as a list
        """

        return self.getSensor("corner")
    
    def getSensor(self, orientation):
        """
        
        Get the sensors of Comcam in (ra, dec) based on the camera list below.
        The reference is at:
        https://confluence.lsstcorp.org/display/LSWUG/Representation+of+a+Camera
        
        Arguments:
            orientation {[string]} -- Orientation of camera to decide which sensor to use.
        
        Returns:
            [list] -- (ra, dec) of four corners of each sensor with the name of sensor as a list.
        """

        # Camera object
        if (orientation == "center"):
            detectorList = ["R:2,2 S:1,1"]
        elif (orientation == "corner"):
            detectorList = ["R:2,2 S:0,2", "R:2,2 S:2,2", "R:2,2 S:0,0", "R:2,2 S:2,0"]
        elif (orientation == "all"):
            detectorList = ["R:2,2 S:0,2", "R:2,2 S:1,2", "R:2,2 S:2,2", "R:2,2 S:0,1", 
                            "R:2,2 S:1,1", "R:2,2 S:2,1", "R:2,2 S:0,0", "R:2,2 S:1,0", 
                            "R:2,2 S:2,0"]

        ra_dec_out = self.getDetectorRaDec(detectorList)

        return ra_dec_out


if __name__ == "__main__":
    pass
