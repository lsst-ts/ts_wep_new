import os
import numpy as np
import unittest

from lsst.ts.wep.SourceProcessor import SourceProcessor
from lsst.ts.wep.bsc.NbrStar import NbrStar
from lsst.ts.wep.Utility import getModulePath, FilterType


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

    def testDmXY2CamXY(self):

        self.sourProc.config(sensorName="R22_S11")
        self.assertEqual(self.sourProc.dmXY2CamXY(4070, 1000), (3000, 4070))

    def testCamXY2DmXY(self):

        self.sourProc.config(sensorName="R22_S11")
        self.assertEqual(self.sourProc.camXY2DmXY(3000, 4070), (4070, 1000))

    def testIsVignette(self):

        isVignette = self.sourProc.isVignette(1.76, 0)
        self.assertTrue(isVignette)

        noVignette = self.sourProc.isVignette(0.2, 0.2)
        self.assertFalse(noVignette)

    def testSimulateImg(self):

        ccdImgIntra, ccdImgExtra = self._simulateImg()

        self.assertEqual(ccdImgIntra.shape, (4072, 2000))
        self.assertNotEqual(np.sum(np.abs(ccdImgIntra)), 0)

    def _simulateImg(self):

        imageFolderPath = os.path.join(self.modulePath, "tests", "testData",
                                       "testImages", "LSST_C_SN26")
        defocalDis = 0.25
        nbrStar = self._generateNbrStar()
        ccdImgIntra, ccdImgExtra = self.sourProc.simulateImg(
            imageFolderPath, defocalDis, nbrStar, FilterType.REF,
            noiseRatio=0)

        return ccdImgIntra, ccdImgExtra

    def _generateNbrStar(self):

        nbrStar = NbrStar()
        nbrStar.starId = {523572575: [],
                          523572679: [523572671]}
        nbrStar.lsstMagG = {523572575: 14.66652,
                            523572671: 16.00000,
                            523572679: 13.25217}
        nbrStar.raDeclInPixel = {523572679: (3966.44, 1022.91),
                                 523572671: (3968.77, 1081.02),
                                 523572575: (3475.48, 479.33)}
        return nbrStar

    def testGetSingleTargetImage(self):

        sglSciNeiImg, allStarPosX, allStarPosY, magRatio, offsetX, offsetY = \
            self._getSingleTargetImage()

        self.assertEqual(sglSciNeiImg.shape, (310, 310))
        self.assertAlmostEqual(allStarPosX[0], 126.98)
        self.assertAlmostEqual(allStarPosX[1], 185.09)
        self.assertAlmostEqual(allStarPosY[0], 206.77)
        self.assertAlmostEqual(allStarPosY[1], 204.44)
        self.assertAlmostEqual(magRatio[0], 0.07959174)
        self.assertEqual(magRatio[1], 1)
        self.assertEqual(offsetX, 792.0)
        self.assertEqual(offsetY, 3762.0)

    def _getSingleTargetImage(self):

        nbrStar = self._generateNbrStar()
        ccdImgIntra, ccdImgExtra = self._simulateImg()
        starIndex = list(nbrStar.getId()).index(523572679)
        sglSciNeiImg, allStarPosX, allStarPosY, magRatio, offsetX, offsetY = \
            self.sourProc.getSingleTargetImage(ccdImgIntra, nbrStar,
                                               starIndex, FilterType.REF)

        return sglSciNeiImg, allStarPosX, allStarPosY, magRatio, offsetX, offsetY

    def testDoDeblending(self):

        sglSciNeiImg, allStarPosX, allStarPosY, magRatio, offsetX, offsetY = \
            self._getSingleTargetImage()

        imgDeblend, realcx, realcy = self.sourProc.doDeblending(
            sglSciNeiImg, allStarPosX, allStarPosY, magRatio)

        self.assertEqual(imgDeblend.shape, (310, 310))
        self.assertLess(np.abs(realcx-184.49), 3)
        self.assertLess(np.abs(realcy-205.00), 3)

        # Get the real camera position x, y after the deblending
        realCameraX = realcx + offsetX
        realCameraY = realcy + offsetY

        # Compared with DM prediction
        nbrStar = self._generateNbrStar()
        raDeclInPixel = nbrStar.getRaDeclInPixel()
        camX, camY = self.sourProc.dmXY2CamXY(raDeclInPixel[523572679][0],
                                              raDeclInPixel[523572679][1])
        delta = np.sqrt((realCameraX-camX)**2 + (realCameraY-camY)**2)
        self.assertLess(delta, 10)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
