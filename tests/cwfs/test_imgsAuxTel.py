import os
import numpy as np
import unittest

from lsst.ts.wep.cwfs.Image import Image
from lsst.ts.wep.cwfs.CompensableImage import CompensableImage
from lsst.ts.wep.cwfs.Instrument import Instrument
from lsst.ts.wep.cwfs.Algorithm import Algorithm
from lsst.ts.wep.Utility import getModulePath, getConfigDir, DefocalType, \
    CamType, CentroidFindType


class TestImgsAuxTel(unittest.TestCase):

    def setUp(self):

        testImageDataDir = os.path.join(getModulePath(), "tests", "testData",
                                        "testImages")
        self.testImgDir = os.path.join(testImageDataDir, "auxTel")
        self.validationDir = os.path.join(testImageDataDir, "validation",
                                          "auxTel")

        self.offset = 80
        self.tolMax = 6
        self.tolRms = 2

    def _runWep(self, imgIntraName, imgExtraName, offset, model):

        # Cut the donut image from input files
        centroidFindType = CentroidFindType.Otsu
        imgIntra = Image(centroidFindType=centroidFindType)
        imgExtra = Image(centroidFindType=centroidFindType)

        imgIntraPath = os.path.join(self.testImgDir, imgIntraName)
        imgExtraPath = os.path.join(self.testImgDir, imgExtraName)

        imgIntra.setImg(imageFile=imgIntraPath)

        imgExtra.setImg(imageFile=imgExtraPath)

        xIntra, yIntra, _ = imgIntra.getCenterAndR()
        imgIntraArray = imgIntra.getImg()[
            int(yIntra)-offset:int(yIntra)+offset,
            int(xIntra)-offset:int(xIntra)+offset]

        xExtra, yExtra, _ = imgExtra.getCenterAndR()
        imgExtraArray = imgExtra.getImg()[
            int(yExtra)-offset:int(yExtra)+offset,
            int(xExtra)-offset:int(xExtra)+offset]

        # Set the images
        fieldXY = (0, 0)
        imgCompIntra = CompensableImage(centroidFindType=centroidFindType)
        imgCompIntra.setImg(fieldXY, DefocalType.Intra, image=imgIntraArray)

        imgCompExtra = CompensableImage(centroidFindType=centroidFindType)
        imgCompExtra.setImg(fieldXY, DefocalType.Extra, image=imgExtraArray)

        # Calculate the wavefront error

        # Set the instrument
        instDir = os.path.join(getConfigDir(), "cwfs", "instData")
        instAuxTel = Instrument(instDir)
        instAuxTel.config(CamType.AuxTel, imgCompIntra.getImgSizeInPix(),
                          announcedDefocalDisInMm=0.8)

        # Set the algorithm
        algoFolderPath = os.path.join(getConfigDir(), "cwfs", "algo")
        algoAuxTel = Algorithm(algoFolderPath)
        algoAuxTel.config("exp", instAuxTel)
        algoAuxTel.runIt(imgCompIntra, imgCompExtra, model)

        return algoAuxTel.getZer4UpInNm()

    def testCase1paraxial(self):

        imgIntraName, imgExtraName = self._getImgsCase1()
        zer4UpNm = self._runWep(imgIntraName, imgExtraName, self.offset,
                                "paraxial")

        ans = self._getDataVerify("case1_auxTel_paraxial.txt")
        self._compareDifferenceWithTol(zer4UpNm, ans)

    def _getImgsCase1(self):

        imgIntraName = "1579925613-16Pup_intra-0-1.fits"
        imgExtraName = "1579925662-16Pup_extra-0-1.fits"

        return imgIntraName, imgExtraName

    def _getDataVerify(self, fileName):

        filePath = os.path.join(self.validationDir, fileName)

        return np.loadtxt(filePath)

    def _compareDifferenceWithTol(self, zer4UpNm, ans):

        diffMax = np.max(np.abs(zer4UpNm - ans))
        self.assertLess(diffMax, self.tolMax)

        diffRms = np.sqrt(np.sum(np.abs(zer4UpNm - ans) ** 2) / 19)
        self.assertLess(diffRms, self.tolRms)

    def testCase1onaxis(self):

        imgIntraName, imgExtraName = self._getImgsCase1()
        zer4UpNm = self._runWep(imgIntraName, imgExtraName, self.offset,
                                "onAxis")
        ans = self._getDataVerify("case1_auxTel_onaxis.txt")
        self._compareDifferenceWithTol(zer4UpNm, ans)

    def testCase2paraxial(self):

        imgIntraName, imgExtraName = self._getImgsCase2()
        zer4UpNm = self._runWep(imgIntraName, imgExtraName, self.offset,
                                "paraxial")

        ans = self._getDataVerify("case2_auxTel_paraxial.txt")
        self._compareDifferenceWithTol(zer4UpNm, ans)

    def _getImgsCase2(self):

        imgIntraName = "1579925833-16Pup_intra-0-1.fits"
        imgExtraName = "1579925882-16Pup_extra-0-1.fits"

        return imgIntraName, imgExtraName

    def testCase2onaxis(self):

        imgIntraName, imgExtraName = self._getImgsCase2()
        zer4UpNm = self._runWep(imgIntraName, imgExtraName, self.offset,
                                "onAxis")

        ans = self._getDataVerify("case2_auxTel_onaxis.txt")
        self._compareDifferenceWithTol(zer4UpNm, ans)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
