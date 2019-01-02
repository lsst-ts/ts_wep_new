import os
import numpy as np
import unittest

from lsst.ts.wep.SourceSelector import SourceSelector
from lsst.ts.wep.SourceProcessor import SourceProcessor
from lsst.ts.wep.Utility import getModulePath, abbrevDectectorName, \
                                expandDetectorName


class testClass(object):
    # Used only for the test class
    pass


class TestSourceProcessor(unittest.TestCase):
    """Test the source processor class."""

    def setUp(self):

        # Get the path of module
        self.modulePath = getModulePath()

        # CCD focal plane file
        focalPlaneFolder = os.path.join(self.modulePath, "tests", "testData")

        # Set the source processor
        self.sourProc = SourceProcessor()

        # Set the configuration
        self.sourProc.config(sensorName="R00_S22_C0",
                             folderPath2FocalPlane=focalPlaneFolder)

    def testInit(self):

        self.assertEqual(self.sourProc.sensorName, "R00_S22_C0")
        self.assertEqual(len(self.sourProc.sensorDimList), 205)
        self.assertEqual(len(self.sourProc.sensorEulerRot), 205)
        self.assertEqual(len(self.sourProc.sensorFocaPlaneInDeg), 205)
        self.assertEqual(len(self.sourProc.sensorFocaPlaneInUm), 205)

        self.assertEqual(self.sourProc.sensorDimList["R00_S22_C0"],
                         (2000, 4072))
        self.assertEqual(self.sourProc.sensorDimList["R22_S11"],
                         (4000, 4072))
        self.assertEqual(self.sourProc.sensorFocaPlaneInDeg["R22_S11"], (0, 0))
        self.assertNotEqual(self.sourProc.sensorFocaPlaneInDeg["R00_S22_C0"], 
                            self.sourProc.sensorFocaPlaneInDeg["R00_S22_C1"])

    def testConfig(self):

        sensorName = "sensorName"
        self.sourProc.config(sensorName=sensorName)

        self.assertEqual(self.sourProc.sensorName, sensorName)

    def testGetEulerZinDeg(self):

        wfsSensorName = "R40_S02_C1"
        eulerZ = self.sourProc.getEulerZinDeg(wfsSensorName)

        self.assertEqual(eulerZ, 90.004585)

    def testCamXYtoFieldXY(self):

        pixelX = 1000
        pixelY = 2036
        fieldX, fieldY = self.sourProc.camXYtoFieldXY(pixelX, pixelY)

        ansFieldX, ansFieldY = \
            self.sourProc.sensorFocaPlaneInDeg[self.sourProc.sensorName]    
        self.assertEqual(fieldX, ansFieldX)
        self.assertEqual(fieldY, ansFieldY)

    def testCamXYtoFieldXYforWfs(self):

        oxR00S22C0, oyR00S22C0 = self._camXYtoFieldXY("R00_S22_C0", 0, 0)
        oxR00S22C1, oyR00S22C1 = self._camXYtoFieldXY("R00_S22_C1", 0, 0)
        oxR40S02C0, oyR40S02C0 = self._camXYtoFieldXY("R40_S02_C0", 0, 0)
        oxR40S02C1, oyR40S02C1 = self._camXYtoFieldXY("R40_S02_C1", 0, 0)
        oxR44S00C0, oyR44S00C0 = self._camXYtoFieldXY("R44_S00_C0", 0, 0)
        oxR44S00C1, oyR44S00C1 = self._camXYtoFieldXY("R44_S00_C1", 0, 0)
        oxR04S20C0, oyR04S20C0 = self._camXYtoFieldXY("R04_S20_C0", 0, 0)
        oxR04S20C1, oyR04S20C1 = self._camXYtoFieldXY("R04_S20_C1", 0, 0)

        # Compare with the same RXX_SYY
        self.assertEqual(oyR00S22C0, oyR00S22C1)
        self.assertEqual(oxR40S02C0, oxR40S02C1)
        self.assertEqual(oyR44S00C0, oyR44S00C1)
        self.assertEqual(oxR04S20C0, oxR04S20C1)

        # Campare with different RXX_SYY
        self.assertEqual((oxR00S22C0 + oxR44S00C0, oyR00S22C0 + oyR44S00C0),
                         (0, 0))
        self.assertEqual((oxR40S02C1 + oxR04S20C1, oyR40S02C1 + oyR04S20C1),
                         (0, 0))

    def _camXYtoFieldXY(self, sensorName, pixelX, pixelY):

        self.sourProc.config(sensorName=sensorName)
        fieldX, fieldY = self.sourProc.camXYtoFieldXY(pixelX, pixelY)

        return fieldX, fieldY

    def testFocalPlaneXY2CamXY(self):
        pass

    # def testDmXY2CamXY(self):

    #     # Define the database and get the neighboring star map
    #     # Address of local database
    #     dbAdress = os.path.join(self.modulePath, "tests", "testData",
    #                             "bsc.db3")

    #     # Use the focal plane as a reference to double check the DM XY to
    #     # Camera XY
    #     # Boresight (RA, Dec) (unit: degree) (0 <= RA <= 360, -90 <= Dec <= 90)
    #     pointing = (20.0, 30.0)

    #     # Camera rotation
    #     cameraRotation = 0.0

    #     # Active filter type
    #     aFilterType = "u"

    #     # Camera type: "lsst" or "comcam"
    #     cameraType = "lsst"

    #     # Set the camera MJD
    #     cameraMJD = 59580.0

    #     # Camera orientation for ComCam ("center" or "corner" or "all")
    #     # Camera orientation for LSSTcam ("corner" or "all")
    #     orientation = "corner"

    #     # Maximum distance in units of radius one donut must be considered as
    #     # a neighbor.
    #     spacingCoefficient = 2.5

    #     # For the defocus = 1.5 mm, the star's radius is 63 pixel.
    #     starRadiusInPixel = 63

    #     # Collect the data from bright star catalog

    #     # Get the neighboring star map
    #     localDb = SourceSelector()
    #     localDb.configSelector(cameraType=cameraType, dbType="LocalDb",
    #                            aFilter=aFilterType)
    #     localDb.configNbrCriteria(starRadiusInPixel, spacingCoefficient,
    #                               maxNeighboringStar=1)
        
    #     localDb.connect(dbAdress)
    #     neighborStarMapLocal, starMapLocal, wavefrontSensorsLocal = \
    #         localDb.getTargetStar(pointing, cameraRotation,
    #                               orientation=orientation)
    #     localDb.disconnect()

    #     # Collect the data
    #     neighborStarMapLocal = neighborStarMapLocal
    #     sensorList = wavefrontSensorsLocal.keys()

    #     # Test the DM XY to Camera XY directly
    #     self.sourProc.config(sensorName="R22_S11")
    #     self.assertEqual(self.sourProc.dmXY2CamXY(4070, 1000), (3000, 4070))

    #     # Test the Camera XY to DM XY directly
    #     self.assertEqual(self.sourProc.camXY2DmXY(3000, 4070), (4070, 1000))

    #     # Change the DM name to camera team

    #     # Test to get the focal plane position
    #     # When writing the test cases, need to add four corners
    #     camera = LsstSimMapper().camera
    #     obs = ObservationMetaData(pointingRA=pointing[0],
    #                               pointingDec=pointing[1], 
    #                               rotSkyPos=cameraRotation, mjd=cameraMJD)

    #     # Veriry the function of DmXY2CamXY() by (ra, decl) to
    #     # (focal X, focal Y) and then to (cam X, cam Y)
    #     # The focal plane coordinate system is the reference of DM
    #     # and Camera teams
    #     abbrevNameList = ["R40_S02_C0", "R00_S22_C0", "R04_S20_C0",
    #                       "R44_S00_C0", "R40_S02_C1", "R00_S22_C1",
    #                       "R04_S20_C1", "R44_S00_C1"]
    #     for abbrevName in abbrevNameList: 
    #         self.sourProc.config(sensorName=abbrevName)

    #         # Transfrom the abbreviated name to full name
    #         fullName = expandDetectorName(abbrevName)

    #         stars = neighborStarMapLocal[fullName]
    #         for starID in stars.RaDecl.keys():
    #             # Transform star (ra, dec) to focal plane coordinate in mm
    #             focalX, focalY = focalPlaneCoordsFromRaDec(
    #                     stars.RaDecl[starID][0], stars.RaDecl[starID][1], 
    #                     obs_metadata=obs, camera=camera)

    #             # Transform focal plane coordinate in mm to pixel position
    #             # in camera coordinate 
    #             # The input unit is "um" instead of "mm"
    #             camX, camY = self.sourProc.focalPlaneXY2CamXY(focalX*1000,
    #                                                           focalY*1000)

    #             # Transform to camera coordinate directly from the DM
    #             # coordinate
    #             camX1, camY1 = self.sourProc.dmXY2CamXY(
    #                                     stars.RaDeclInPixel[starID][0],
    #                                     stars.RaDeclInPixel[starID][1])

    #             # Do the comparison
    #             delta = np.sqrt( (camX-camX1)**2 + (camY-camY1)**2 )
    #             self.assertLess(delta, 10)

    # def testDeblending(self):

    #     # Donut image folder
    #     imageFolder = os.path.join(self.modulePath, "tests", "testData",
    #                                "testImages")
    #     donutImageFolder = "LSST_C_SN26"

    #     # Give the path to the image folder
    #     imageFolderPath = os.path.join(imageFolder, donutImageFolder)

    #     # Generate the simulated image
    #     defocalDis = 0.25
    #     afilter = "u"
    #     self.sourProc.config(sensorName="R04_S20_C1")

    #     # Create a mocked neighboring star map
    #     neighborStarMap = testClass()
    #     SimobjID = {523572575: [], 
    #                 523572679: [523572671]}
    #     setattr(neighborStarMap, "SimobjID", SimobjID)
    #     LSSTMagU = {523572575: 14.66652, 
    #                 523572671: 16.00000, 
    #                 523572679: 13.25217}
    #     setattr(neighborStarMap, "LSSTMagU", LSSTMagU)
    #     RaDeclInPixel = {523572679: (3966.4462129591157, 1022.9153550029878), 
    #                      523572671: (3968.7766808905071, 1081.0241833910586),
    #                      523572575: (3475.4821263515223, 479.33235991200854)}
    #     setattr(neighborStarMap, "RaDeclInPixel", RaDeclInPixel)

    #     # Simulate the image
    #     ccdImgIntra, ccdImgExtra = self.sourProc.simulateImg(
    #         imageFolderPath, defocalDis, neighborStarMap, afilter,
    #         noiseRatio=0)

    #     # Check the dimension and the image is not zero
    #     self.assertEqual(ccdImgIntra.shape, (4072, 2000))
    #     self.assertNotEqual(np.sum(np.abs(ccdImgIntra)), 0)

    #     # Show the image
    #     # poltExposureImage(ccdImgIntra, name="Intra focal image", scale="log",
    #     #                   cmap=None)
    #     # poltExposureImage(ccdImgIntra, name="Intra focal image",
    #     #                   scale="linear", cmap=None)

    #     # Get the images of one bright star map
    #     starIndex = list(neighborStarMap.SimobjID).index(523572679)
    #     singleSciNeiImg, allStarPosX, allStarPosY, magRatio, offsetX, offsetY = \
    #         self.sourProc.getSingleTargetImage(ccdImgIntra, neighborStarMap,
    #                                            starIndex, afilter)

    #     # Show the image
    #     # poltExposureImage(singleSciNeiImg, name="Single intra focal image",
    #     #                   scale="log", cmap=None)
    #     # poltExposureImage(singleSciNeiImg, name="Single intra focal image",
    #     #                   scale="linear", cmap=None)

    #     # Do the deblending and determine the real position on camera
    #     imgDeblend, realcx, realcy = self.sourProc.doDeblending(
    #             singleSciNeiImg, allStarPosX, allStarPosY, magRatio)

    #     # Show the deblended image
    #     # poltExposureImage(imgDeblend, name="Deblended image", scale="log",
    #     #                   cmap=None)
    #     # poltExposureImage(imgDeblend, name="Deblended image", scale="linear",
    #     #                   cmap=None)

    #     # Get the real camera position x, y after the deblending
    #     realCameraX = realcx + offsetX
    #     realCameraY = realcy + offsetY

    #     # Compared with DM prediction
    #     dmX, dmY = self.sourProc.dmXY2CamXY(RaDeclInPixel[523572679][0],
    #                                         RaDeclInPixel[523572679][1])
    #     delta = np.sqrt( (realCameraX-dmX)**2 + (realCameraY-dmY)**2 )
    #     self.assertLess(delta, 10)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
