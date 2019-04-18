import os
import numpy as np
from astropy.io import fits
import shutil
import unittest

from lsst.ts.wep.cwfs.Tool import ZernikeAnnularFit
from lsst.ts.wep.CamDataCollector import CamDataCollector
from lsst.ts.wep.CamIsrWrapper import CamIsrWrapper
from lsst.ts.wep.SourceProcessor import SourceProcessor
from lsst.ts.wep.SourceSelector import SourceSelector
from lsst.ts.wep.WfEstimator import WfEstimator
from lsst.ts.wep.WepController import WepController

from lsst.ts.wep.Utility import getModulePath, FilterType, CamType, BscDbType,\
    runProgram, getConfigDir


class TestWepControllerMonolithic(unittest.TestCase):
    """Test the WepController class."""

    def setUp(self):

        self.modulePath = getModulePath()
        self.dataDir = os.path.join(self.modulePath, "tests", "tmp")
        self.isrDir = os.path.join(self.dataDir, "input")
        self.opdDir = os.path.join(self.modulePath, "tests", "testData",
                                   "opdOutput", "9005000")

        self._makeDir(self.isrDir)

        # Configurate the WEP components
        dataCollector = CamDataCollector(self.isrDir)
        isrWrapper = CamIsrWrapper(self.isrDir)
        sourSelc = self._configSourceSelector()
        sourProc = self._configSourceProcessor()
        wfEsti = self._configWfEstimator()

        # Instantiate the WEP controller
        self.wepCntlr = WepController(dataCollector, isrWrapper, sourSelc,
                                      sourProc, wfEsti)

        # Intemediate data used in the test
        self.filter = FilterType.REF

        self.neighborStarMap = dict()
        self.starMap = dict()
        self.wavefrontSensors = dict()

        self.wfsImgMap = dict()
        self.donutMap = dict()
        self.masterDonutMap = dict()

    def _makeDir(self, directory):

        if (not os.path.exists(directory)):
            os.makedirs(directory)

    def _configSourceSelector(self):

        sourSelc = SourceSelector(CamType.ComCam, BscDbType.LocalDbForStarFile)

        # Set the criteria of neighboring stars
        starRadiusInPixel = 63
        spacingCoefficient = 2.5
        maxNeighboringStar = 1
        sourSelc.configNbrCriteria(starRadiusInPixel, spacingCoefficient,
                                   maxNeighboringStar=maxNeighboringStar)

        # Connest the database
        dbAdress = os.path.join(self.modulePath, "tests", "testData",
                                "bsc.db3")
        sourSelc.connect(dbAdress)

        return sourSelc

    def _configSourceProcessor(self):

        folderPath2FocalPlane = os.path.join(self.modulePath, "tests",
                                             "testData")
        sourProc = SourceProcessor(folderPath2FocalPlane=folderPath2FocalPlane)

        return sourProc

    def _configWfEstimator(self):

        cwfsConfigDir = os.path.join(getConfigDir(), "cwfs")
        instDir = os.path.join(cwfsConfigDir, "instData")
        algoDir = os.path.join(cwfsConfigDir, "algo")
        wfEsti = WfEstimator(instDir, algoDir)

        # Use the comcam to calculate the LSST central raft image
        # with 1.5 mm defocal distance
        wfEsti.config(solver="exp", instName="comcam",
                      opticalModel="offAxis", defocalDisInMm=1.5,
                      sizeInPix=160, debugLevel=0)

        return wfEsti

    def tearDown(self):

        self.wepCntlr.getSourSelc().disconnect()

        shutil.rmtree(self.dataDir)

    def testGetDataCollector(self):

        self.assertTrue(isinstance(self.wepCntlr.getDataCollector(),
                                   CamDataCollector))

    def testGetIsrWrapper(self):

        self.assertTrue(isinstance(self.wepCntlr.getIsrWrapper(),
                                   CamIsrWrapper))

    def testGetSourSelc(self):

        self.assertTrue(isinstance(self.wepCntlr.getSourSelc(),
                                   SourceSelector))

    def testGetSourProc(self):

        self.assertTrue(isinstance(self.wepCntlr.getSourProc(),
                                   SourceProcessor))

    def testGetWfEsti(self):

        self.assertTrue(isinstance(self.wepCntlr.getWfEsti(),
                                   WfEstimator))

    def testGetButlerWrapper(self):

        self.assertEqual(self.wepCntlr.getButlerWrapper(), None)

    def testMonolithicSteps(self):
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

    def step1_genCalibsAndIngest(self):

        # Generate the fake flat images
        fakeFlatDir = os.path.join(self.dataDir, "fake_flats")
        self._makeDir(fakeFlatDir)

        detector = "R22_S11 R22_S10"
        self._genFakeFlat(fakeFlatDir, detector)

        # Generate the PhoSim mapper
        self.wepCntlr.getDataCollector().genPhoSimMapper()

        # Do the ingestion
        calibFiles = os.path.join(fakeFlatDir, "*")
        self.wepCntlr.getDataCollector().ingestCalibs(calibFiles)

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

        intraImgFiles = os.path.join(getModulePath(), "tests", "testData",
                                     "phosimOutput", "realComCam",
                                     "repackagedFiles", "intra", "*.fits")
        extraImgFiles = os.path.join(getModulePath(), "tests", "testData",
                                     "phosimOutput", "realComCam",
                                     "repackagedFiles", "extra", "*.fits")

        self.wepCntlr.getDataCollector().ingestImages(intraImgFiles)
        self.wepCntlr.getDataCollector().ingestImages(extraImgFiles)

    def step3_doIsr(self):

        fileName = "isr_config.py"
        self.wepCntlr.getIsrWrapper().config(doFlat=True, fileName=fileName)

        rerunName = "run1"
        self.wepCntlr.getIsrWrapper().doISR(self.isrDir, rerunName=rerunName)

    def step4_setButlerInputsPath(self):

        inputs = os.path.join(self.isrDir, "rerun", "run1")
        self.wepCntlr.setPostIsrCcdInputs(inputs)

    def step5_getTargetStar(self):

        # Set the observation meta data
        ra = 0.0
        dec = 0.0
        rotSkyPos = 0.0
        self.wepCntlr.getSourSelc().setObsMetaData(ra, dec, rotSkyPos)

        # Set the filter
        self.wepCntlr.getSourSelc().setFilter(self.filter)

        # Get the target star by file
        skyFilePath = os.path.join(self.modulePath, "tests", "testData",
                                   "phosimOutput", "realComCam",
                                   "skyComCamInfo.txt")
        neighborStarMap, starMap, wavefrontSensors = \
            self.wepCntlr.getSourSelc().getTargetStarByFile(skyFilePath, offset=0)

        # Assign the data for the following steps to use
        self.neighborStarMap = neighborStarMap
        self.starMap = starMap
        self.wavefrontSensors = wavefrontSensors

        # Do the assertion
        self.assertEqual(len(neighborStarMap), 2)
        self.assertEqual(len(starMap), 2)
        self.assertEqual(len(wavefrontSensors), 2)

    def step6_getPostIsrDefocalImgMap(self):

        sensorNameList = list(self.wavefrontSensors)

        intraObsId = 9005001
        extraObsId = 9005000
        obsIdList = [intraObsId, extraObsId]

        wfsImgMap = self.wepCntlr.getPostIsrImgMapByPistonDefocal(
            sensorNameList, obsIdList)

        # Assign the data for the following steps to use
        self.wfsImgMap = wfsImgMap

        # Do the assertion
        self.assertEqual(len(wfsImgMap), 2)

    def step7_getDonutMap(self):

        self.donutMap = self.wepCntlr.getDonutMap(
            self.neighborStarMap, self.wfsImgMap, self.filter,
            doDeblending=False)

        # Do the assertion
        for sensor, donutList in self.donutMap.items():
            self.assertEqual(len(donutList), 2)

    def step8a_calcWfErr(self):

        self.donutMap = self.wepCntlr.calcWfErr(self.donutMap)

        # Do the assertion
        for sensor, donutList in self.donutMap.items():
            for donut in donutList:
                wfErr = donut.getWfErr()
                self.assertEqual(wfErr.argmax(), 2)
                self.assertGreater(wfErr.max(), 100)

        # Compare with OPD
        donutList = []
        donutList.extend(self.donutMap["R:2,2 S:1,1"])
        donutList.extend(self.donutMap["R:2,2 S:1,0"])

        for aId in range(4):
            wfErr = donutList[aId].getWfErr()
            zkOfOpd = self._getZkInNmFromOpd(aId)[3:]
            delta = np.abs(wfErr[7]-zkOfOpd[7])
            self.assertLess(delta, 5)

    def _getZkInNmFromOpd(self, opdId):

        # Get the OPD data
        opdFileName = "opd_9005000_%d.fits.gz" % opdId
        opdFitsFile = os.path.join(self.opdDir, opdFileName)
        opd = fits.getdata(opdFitsFile)

        # Get x-, y-coordinate in the OPD image
        opdSize = opd.shape[0]
        opdGrid1d = np.linspace(-1, 1, opdSize)
        opdx, opdy = np.meshgrid(opdGrid1d, opdGrid1d)

        # Fit the OPD map with Zk
        idx = (opd != 0)
        znTerms = 22
        obscuration = 0.61
        zk = ZernikeAnnularFit(opd[idx], opdx[idx], opdy[idx], znTerms,
                               obscuration)

        # Return the unit in nm (the unit in OPD is um)
        return zk * 1e3

    def step8b_calcAvgWfErrOnSglCcd(self):

        avgErrMap = dict()
        for sensor, donutList in self.donutMap.items():
            avgErr = self.wepCntlr.calcAvgWfErrOnSglCcd(donutList)
            avgErrMap[sensor] = avgErr

            # Do the assertion
            self.assertEqual(avgErr.argmax(), 2)
            self.assertGreater(avgErr.max(), 100)

        # Compare with the central OPD
        wfErrOnR22S11 = avgErrMap["R:2,2 S:1,1"]
        zkOfOpdOnR22S11 = self._getZkInNmFromOpd(4)[3:]
        deltaOnR22S11 = np.abs(wfErrOnR22S11[7] - zkOfOpdOnR22S11[7])
        self.assertLess(deltaOnR22S11, 5)

        wfErrOnR22S10 = avgErrMap["R:2,2 S:1,0"]
        zkOfOpdOnR22S10 = self._getZkInNmFromOpd(5)[3:]
        deltaOnR22S10 = np.abs(wfErrOnR22S10[7] - zkOfOpdOnR22S10[7])
        self.assertLess(deltaOnR22S10, 5)

    def step9a_genMasterDonut(self):

        self.masterDonutMap = self.wepCntlr.genMasterDonut(
            self.donutMap, zcCol=np.zeros(22))

        # Do the assertion
        self.assertEqual(len(self.masterDonutMap), 2)
        for sensor, masterDonutList in self.masterDonutMap.items():
            self.assertEqual(len(masterDonutList), 1)

    def step9b_calcWfErrOfMasterDonut(self):

        self.masterDonutMap = self.wepCntlr.calcWfErr(self.masterDonutMap)

        masterDonutList = self.masterDonutMap["R:2,2 S:1,1"]
        wfErr = masterDonutList[0].getWfErr()

        # Do the assertion
        self.assertEqual(wfErr.argmax(), 2)
        self.assertGreater(wfErr.max(), 100)

        # Compare with the central OPD
        zkOfOpdOnR22S11 = self._getZkInNmFromOpd(4)[3:]
        deltaOnR22S11 = np.abs(wfErr[7] - zkOfOpdOnR22S11[7])
        self.assertLess(deltaOnR22S11, 30)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
