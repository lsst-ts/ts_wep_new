import os

import astropy
import lsst.utils.tests
import numpy as np
from lsst.daf.butler import Butler
from lsst.ts.wep.task.fitDonutRadiusTask import (
    FitDonutRadiusTask,
    FitDonutRadiusTaskConfig,
)
from lsst.ts.wep.utils import (
    getModulePath,
    runProgram,
    writeCleanUpRepoCmd,
    writePipetaskCmd,
)


class TestFitDonutRadiusTaskScienceSensor(lsst.utils.tests.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Generate donut stamps needed for task.
        """

        moduleDir = getModulePath()
        cls.testDataDir = os.path.join(moduleDir, "tests", "testData")
        testPipelineConfigDir = os.path.join(cls.testDataDir, "pipelineConfigs")
        cls.repoDir = os.path.join(cls.testDataDir, "gen3TestRepo")
        cls.runNameScience = "run1"
        cls.baseRunNameScience = "pretest_run_science"
        cls.runNameCwfs = "run2"
        cls.baseRunNameCwfs = "pretest_run_cwfs"

        # Check that run doesn't already exist due to previous improper cleanup
        butler = Butler(cls.repoDir)
        registry = butler.registry
        collectionsList = list(registry.queryCollections())
        for runName in [cls.runNameScience, cls.runNameCwfs]:
            if runName in collectionsList:
                cleanUpCmd = writeCleanUpRepoCmd(cls.repoDir, runName)
                runProgram(cleanUpCmd)

        collections = "refcats/gen2,LSSTCam/calib,LSSTCam/raw/all"
        instrument = "lsst.obs.lsst.LsstCam"
        cls.cameraName = "LSSTCam"

        # First run science sensors
        pipelineYaml = os.path.join(
            testPipelineConfigDir, "testFitDonutRadiusFamPipeline.yaml"
        )

        if "pretest_run_science" in collectionsList:
            pipelineYaml += "#fitDonutRadiusTask"
            collections += ",pretest_run_science"
            cls.baseRunNameScience = "pretest_run_science"
        pipeCmd = writePipetaskCmd(
            cls.repoDir,
            cls.runNameScience,
            instrument,
            collections,
            pipelineYaml=pipelineYaml,
        )
        # The only available science detectors here are 93 and 94
        pipeCmd += (
            ' -d "exposure IN (4021123106001, 4021123106002) and detector.id = 94" '
        )
        runProgram(pipeCmd)

        # Run CWFS sensors
        collections = "refcats/gen2,LSSTCam/calib,LSSTCam/raw/all"
        pipelineYaml = os.path.join(
            testPipelineConfigDir, "testFitDonutRadiusCwfsPipeline.yaml"
        )

        if "pretest_run_cwfs" in collectionsList:
            pipelineYaml += "#fitDonutRadiusTask"
            collections += ",pretest_run_cwfs"
            cls.baseRunNameCwfs = "pretest_run_cwfs"
        pipeCmd = writePipetaskCmd(
            cls.repoDir,
            cls.runNameCwfs,
            instrument,
            collections,
            pipelineYaml=pipelineYaml,
        )
        # Use only one corner sensor pair
        pipeCmd += ' -d "exposure = 4021123106000"'
        runProgram(pipeCmd)

    def setUp(self):
        self.config = FitDonutRadiusTaskConfig()
        self.task = FitDonutRadiusTask(config=self.config, name="Base Task")

        self.butler = Butler(self.repoDir)
        self.registry = self.butler.registry

        self.dataIdExtraScience = {
            "instrument": "LSSTCam",
            "detector": 94,
            "exposure": 4021123106001,
            "visit": 4021123106001,
        }
        self.dataIdExtraCwfs = {
            "instrument": "LSSTCam",
            "detector": 191,
            "exposure": 4021123106000,
            "visit": 4021123106000,
        }

    def testValidateConfigs(self):
        # Test default value
        self.OrigTask = FitDonutRadiusTask(config=self.config, name="Orig Task")
        self.assertEqual(self.OrigTask.config.widthMultiplier, 0.8)
        self.assertEqual(self.OrigTask.config.filterSigma, 3)
        self.assertEqual(self.OrigTask.config.minPeakWidth, 5)
        self.assertEqual(self.OrigTask.config.minPeakHeight, 0.3)

        # Test changing configs
        self.config.widthMultiplier = 2.00
        self.config.filterSigma = 4
        self.config.minPeakWidth = 7
        self.config.minPeakHeight = 0.8
        self.ModifiedTask = FitDonutRadiusTask(config=self.config, name="Mod Task")
        self.assertEqual(self.ModifiedTask.config.widthMultiplier, 2.00)
        self.assertEqual(self.ModifiedTask.config.filterSigma, 4)
        self.assertEqual(self.ModifiedTask.config.minPeakWidth, 7)
        self.assertEqual(self.ModifiedTask.config.minPeakHeight, 0.8)

    def testTaskRunScienceSensor(self):

        donutStampsExtra = self.butler.get(
            "donutStampsExtra",
            dataId=self.dataIdExtraScience,
            collections=[self.baseRunNameScience],
        )
        donutStampsIntra = self.butler.get(
            "donutStampsIntra",
            dataId=self.dataIdExtraScience,
            collections=[self.baseRunNameScience],
        )

        taskOut = self.task.run(donutStampsExtra, donutStampsIntra)
        self.assertEqual(type(taskOut.donutRadiiTable), astropy.table.table.QTable)

        self.assertEqual(len(taskOut.donutRadiiTable), 6)

        # test that the expected columns are present
        expected_columns = [
            "VISIT",
            "DFC_TYPE",
            "DET_NAME",
            "DFC_DIST",
            "RADIUS",
            "X_LEFT_EDGE",
            "X_RIGHT_EDGE",
            "FAIL_FLAG",
        ]
        self.assertLessEqual(
            set(expected_columns), set(taskOut.donutRadiiTable.colnames)
        )
        # test that correct detector names are present
        self.assertEqual(
            set(np.unique(taskOut.donutRadiiTable["DET_NAME"].value)),
            set(["R22_S11"]),
        )
        # test that correct visits ids are present
        self.assertEqual(
            set(np.unique(taskOut.donutRadiiTable["VISIT"].value)),
            set([4021123106001, 4021123106002]),
        )
        # test that the mean radius is correct
        self.assertFloatsAlmostEqual(
            69.69552203692189,
            np.mean(taskOut.donutRadiiTable["RADIUS"]),
            rtol=1e-1,
            atol=1e-1,
        )

    def testPipelineRunScienceSensor(self):
        donutRadiiTable = self.butler.get(
            "donutRadiiTable",
            dataId=self.dataIdExtraScience,
            collections=[self.runNameScience],
        )
        # as above, check that the table was made automatically
        self.assertEqual(type(donutRadiiTable), astropy.table.table.QTable)

        self.assertEqual(len(donutRadiiTable), 6)

    def testTaskRunCwfs(self):

        donutStampsExtra = self.butler.get(
            "donutStampsExtra",
            dataId=self.dataIdExtraCwfs,
            collections=[self.baseRunNameCwfs],
        )
        donutStampsIntra = self.butler.get(
            "donutStampsIntra",
            dataId=self.dataIdExtraCwfs,
            collections=[self.baseRunNameCwfs],
        )

        taskOut = self.task.run(donutStampsExtra, donutStampsIntra)
        self.assertEqual(type(taskOut.donutRadiiTable), astropy.table.table.QTable)

        self.assertEqual(len(taskOut.donutRadiiTable), 4)

        expected_columns = [
            "VISIT",
            "DFC_TYPE",
            "DET_NAME",
            "DFC_DIST",
            "RADIUS",
            "X_LEFT_EDGE",
            "X_RIGHT_EDGE",
            "FAIL_FLAG",
        ]
        self.assertLessEqual(
            set(expected_columns), set(taskOut.donutRadiiTable.colnames)
        )
        # test that correct detector names are present
        self.assertEqual(
            set(np.unique(taskOut.donutRadiiTable["DET_NAME"].value)),
            set(
                [
                    "R00_SW0",
                    "R00_SW1",
                ]
            ),
        )
        # test that correct visit id is present
        self.assertEqual(
            set(
                np.unique(
                    np.asarray(taskOut.donutRadiiTable["VISIT"].value).astype(int)
                )
            ),
            set([4021123106000]),
        )
        self.assertFloatsAlmostEqual(
            66.38605539619834,
            np.mean(taskOut.donutRadiiTable["RADIUS"]),
            rtol=1e-1,
            atol=1e-1,
        )

    def testPipelineRunCwfs(self):
        donutRadiiTable = self.butler.get(
            "donutRadiiTable",
            dataId=self.dataIdExtraCwfs,
            collections=[self.runNameCwfs],
        )
        # check that the table was made automatically
        self.assertEqual(type(donutRadiiTable), astropy.table.table.QTable)
        self.assertEqual(len(donutRadiiTable), 4)

    @classmethod
    def tearDownClass(cls):
        for runName in [cls.runNameScience, cls.runNameCwfs]:
            cleanUpCmd = writeCleanUpRepoCmd(cls.repoDir, runName)
            runProgram(cleanUpCmd)
