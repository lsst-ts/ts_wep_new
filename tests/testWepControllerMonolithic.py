import os
import shutil
import unittest

from lsst.ts.wep.CamDataCollector import CamDataCollector
from lsst.ts.wep.CamIsrWrapper import CamIsrWrapper
from lsst.ts.wep.SourceProcessor import SourceProcessor
from lsst.ts.wep.SourceSelector import SourceSelector
from lsst.ts.wep.WfEstimator import WfEstimator
from lsst.ts.wep.WepController import WepController

from lsst.ts.wep.Utility import getModulePath, FilterType, CamType, BscDbType,\
                                runProgram


class TestWepControllerMonolithic(unittest.TestCase):
    """Test the WepController class."""

    def setUp(self):

        self.modulePath = getModulePath()
        self.dataDir = os.path.join(self.modulePath, "tests", "tmp")
        self.isrDir = os.path.join(self.dataDir, "input")

        self._makeDir(self.isrDir)

        # Configurate the WEP components
        dataCollector = CamDataCollector(self.isrDir)
        isrWrapper = CamIsrWrapper(self.isrDir)
        sourSelc = SourceSelector(CamType.ComCam, BscDbType.LocalDbForStarFile)
        sourProc = self._configSourceProcessor()
        wfsEsti = self._configWfEstimator()

        # Instantiate the WEP controller
        self.wepCntlr = WepController(dataCollector, isrWrapper, sourSelc,
                                      sourProc, wfsEsti)

    def _makeDir(self, directory):

        if (not os.path.exists(directory)):
            os.makedirs(directory)

    def _configSourceProcessor(self):

        folderPath2FocalPlane = os.path.join(self.modulePath, "tests",
                                             "testData")
        sourProc = SourceProcessor(folderPath2FocalPlane=folderPath2FocalPlane)

        return sourProc

    def _configWfEstimator(self):

        instruFolderPath = os.path.join(self.modulePath, "configData", "cwfs",
                                        "instruData")
        algoFolderPath = os.path.join(self.modulePath, "configData", "cwfs",
                                      "algo")
        wfsEsti = WfEstimator(instruFolderPath, algoFolderPath)

        return wfsEsti

    def testSteps(self):
        """Do the test based on the steps defined in the child class."""

        for name, step in self._steps():
            try:
                step()
            except Exception as e:
                self.fail("{} failed ({}: {})".format(step, type(e), e))

    def _steps(self):
        """Sort the order of test steps.

        Yields
        ------
        str
            Step function name.
        func
            Step function generator.
        """

        for name in sorted(dir(self)):
            if name.startswith("step"):
                yield name, getattr(self, name)

    @unittest.skip
    def tearDown(self):

        shutil.rmtree(self.dataDir)

    def step1_genCalibsAndIngest(self):

        # Generate the fake flat images
        fakeFlatDir = os.path.join(self.dataDir, "fake_flats")
        self._makeDir(fakeFlatDir)

        detector = "R22_S11 R22_S10"
        self._genFakeFlat(fakeFlatDir, detector)

        # Generate the PhoSim mapper
        self.wepCntlr.dataCollector.genPhoSimMapper()

        # Do the ingestion
        calibFiles = os.path.join(fakeFlatDir, "*")
        self.wepCntlr.dataCollector.ingestCalibs(calibFiles)

    def _genFakeFlat(self, fakeFlatDir, detector):
        
        currWorkDir = os.getcwd()

        os.chdir(fakeFlatDir)
        self._makeFakeFlat(detector)
        os.chdir(currWorkDir)

    def _makeFakeFlat(self, detector):

        command = "makeGainImages.py"
        argstring = "--detector_list %s" % detector
        runProgram(command, argstring=argstring)

    def step2_ingestExp(self):

        print("Step 2.")

    def step3_doIsr(self):
        pass

    #     # Instantiate the WEP controller
    #     self.wepCntlr = WEPController()

    # def testCornerWfsFunction(self):

    #     # Test to get the list of corner wavefront sensors
    #     wfsList = self.wepCntlr.getWfsList()
    #     self.assertEqual(len(wfsList), 8)

    #     # Instintiate the components
    #     sourSelc = SourceSelector()
    #     dataCollector = WFDataCollector()
    #     sourProc = SourceProcessor()

    #     instruFolderPath = os.path.join(self.modulePath, "algoData", "cwfs", "instruData")
    #     algoFolderPath = os.path.join(self.modulePath, "algoData", "cwfs", "algo")
    #     wfsEsti = WFEstimator(instruFolderPath, algoFolderPath)

    #     # Configurate the source selector
    #     cameraType = "lsst"
    #     dbType = "LocalDb"
    #     aFilter = "g"
    #     cameraMJD = 59580.0

    #     sourSelc.configSelector(cameraType=cameraType, dbType=dbType, aFilter=aFilter, 
    #                             cameraMJD=cameraMJD)

    #     # Set the criteria of neighboring stars
    #     starRadiusInPixel = 63
    #     spacingCoefficient = 2.5
    #     sourSelc.configNbrCriteria(starRadiusInPixel, spacingCoefficient)

    #     # Configurate the WFS data collector
    #     # Data butler does not support the corner WFS at this moment.
    #     pathOfRawData = os.path.join(self.modulePath, "test", "phosimOutput")
    #     destinationPath = butlerInputs = butlerOutputs = os.path.join(self.modulePath, "test")
    #     dataCollector.config(pathOfRawData=pathOfRawData, destinationPath=destinationPath)

    #     # Configurate the source processor
    #     focalPlaneFolder = os.path.join(self.modulePath, "test")
    #     sourProc.config(donutRadiusInPixel=starRadiusInPixel, folderPath2FocalPlane=focalPlaneFolder, 
    #                     pixel2Arcsec=0.2)

    #     # Configurate the wavefront estimator
    #     defocalDisInMm = None
        
    #     # Size of donut in pixel for corner WFS
    #     sizeInPix = 120
    #     wfsEsti.config(solver="exp", instName=cameraType, opticalModel="offAxis", 
    #                     defocalDisInMm=defocalDisInMm, sizeInPix=sizeInPix)

    #     # Configurate the WEP controller
    #     self.wepCntlr.config(sourSelc=sourSelc, dataCollector=dataCollector, 
    #                         sourProc=sourProc, wfsEsti=wfsEsti)

    #     # Test the configuration
    #     self.assertTrue(isinstance(self.wepCntlr.wfsEsti, WFEstimator))

    #     # Get the target stars by file

    #     # Set the database address
    #     dbAdress = os.path.join(self.modulePath, "test", "bsc.db3")

    #     # Do the query
    #     pointing = (0,0)
    #     cameraRotation = 0.0
    #     skyInfoFilePath = os.path.join(self.modulePath, "test", "phosimOutput", "realWfs", "output", 
    #                                    "skyWfsInfo.txt")

    #     camOrientation = "corner"
    #     neighborStarMap, starMap, wavefrontSensors = self.wepCntlr.getTargetStarByFile(dbAdress, 
    #                                                     skyInfoFilePath, pointing, cameraRotation, 
    #                                             orientation=camOrientation, tableName="TempTable")
    #     self.assertEqual(len(starMap), 8)

    #     starData = starMap["R:0,0 S:2,2,A"]
    #     self.assertEqual(len(starData.SimobjID), 2)

    #     # Get the available sensor name list
    #     sensorNameList = list(starMap.keys())

    #     # Get the eimage
    #     wfsDir = os.path.join("realWfs", "output")
    #     wfsImgMap = self.wepCntlr.getPostISRDefocalImgMap(sensorNameList, wfsDir=wfsDir)

    #     wfsImg = wfsImgMap["R:0,0 S:2,2,A"]
    #     self.assertNotEqual(np.sum(wfsImg.intraImg), None)
    #     self.assertEqual(wfsImg.extraImg, None)

    #     # Get the donut map
    #     donutMap = self.wepCntlr.getDonutMap(neighborStarMap, wfsImgMap, aFilter, 
    #                                         doDeblending=False, sglDonutOnly=True)
        
    #     donutList = donutMap["R:0,0 S:2,2,A"]
    #     self.assertEqual(len(donutList), 2)

    #     donutImg = donutList[0]
    #     self.assertNotEqual(np.sum(donutImg.intraImg), None)
    #     self.assertEqual(donutImg.extraImg, None)
    #     self.assertEqual(donutImg.starId, 6)
    #     self.assertEqual(int(donutImg.pixelX), 506)
    #     self.assertEqual(int(donutImg.pixelY), 1008)

    #     # Calculate the wavefront error for the individual donut
    #     partDonutMap = dict()
    #     partDonutMap["R:0,0 S:2,2,A"] = donutMap["R:0,0 S:2,2,A"]
    #     partDonutMap["R:0,0 S:2,2,B"] = donutMap["R:0,0 S:2,2,B"]
        
    #     partDonutMap = self.wepCntlr.calcWfErr(partDonutMap)
        
    #     donutList = partDonutMap["R:0,0 S:2,2,A"]
    #     donutImg = donutList[0]
    #     self.assertEqual(len(donutImg.zer4UpNm), 19)

    #     # Test the weighting ratio
    #     weightingRatio = self.wepCntlr.calcWeiRatio(donutList)
    #     self.assertEqual(np.sum(weightingRatio), 1)
    #     self.assertEqual(weightingRatio[0], 0.5)

    #     # Test to calculate the average wavefront error
    #     avgErr = self.wepCntlr.calcSglAvgWfErr(donutList)
    #     ans = donutList[0].zer4UpNm*weightingRatio[0] + donutList[1].zer4UpNm*weightingRatio[1]
    #     self.assertEqual(np.sum(avgErr), np.sum(ans))


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
