# -*- coding: utf-8 -*-

# This script is for the unit test.

import os, time, unittest
import os
import numpy as np
from runWEP import runWEP
from lsst.ts.wep.Utility import getModulePath

class TestFile(object): 
    
    def __init__(self, ImageFolderName, ImageName, FieldXY, UseAlgorithm, Orientation, ValidationPath):
        
        self.imageFolderName = ImageFolderName
        self.imageName = ImageName
        self.fieldXY = [FieldXY, FieldXY]
        self.useAlgorithm = UseAlgorithm
        self.orientation = Orientation
        
        # Get the refFilePath
        refFileName = ImageFolderName + "_" + ImageName + UseAlgorithm + ".txt"
        self.refFilePath = os.path.join(ValidationPath,refFileName)
        
class DataWEP(object):
    
    def __init__(self, InstruFolder, AlgoFolderPath, InstruName, ImageFolder, TestFile):
    
        self.instruFolder = InstruFolder
        self.algoFolderPath = AlgoFolderPath
        self.instruName = InstruName
        self.useAlgorithm = TestFile.useAlgorithm
        self.imageFolderPath = os.path.join(ImageFolder,TestFile.imageFolderName)
        self.intra_image_name = TestFile.imageName + "intra.txt"
        self.extra_image_name = TestFile.imageName + "extra.txt"
        self.fieldXY = TestFile.fieldXY
        self.orientation = TestFile.orientation
        self.refFilePath = TestFile.refFilePath

class TestWEP(unittest.TestCase):    
    
    def setUp(self):

        # Get the path of module
        modulePath = getModulePath()
        
        # Restart time
        self.startTime = time.time()
        self.difference = 0
        self.validationDir = os.path.join(modulePath, "test", "testImages", "validation")
    
    def tearDown(self):
        
        # Calculate the time of test case
        t = time.time() - self.startTime
        print("%s: %.3f s. Differece is %.3f." % (self.id(), t, self.difference))
    
    def generateTestCase(self, ImageFolderName, ImageName ,FieldXY, UseAlgorithm, 
                         Orientation, ValidationPath):

        case_test = TestFile(ImageFolderName, ImageName, FieldXY, UseAlgorithm, Orientation, 
                             ValidationPath)  
    
        case = DataWEP(InstruFolder, AlgoFolderPath, InstruName, ImageFolderPath, 
                       case_test)
    
        return case
    
    def compareCalculation(self, DataWEP, tor):
                
        # Run WEP to get Zk
        zer4UpNm = runWEP(DataWEP.instruFolder, DataWEP.algoFolderPath, DataWEP.instruName, 
                          DataWEP.useAlgorithm, DataWEP.imageFolderPath, DataWEP.intra_image_name, 
                          DataWEP.extra_image_name, DataWEP.fieldXY, DataWEP.orientation)
        
        # Load the reference data
        refZer4UpNm = np.loadtxt(DataWEP.refFilePath)

        # Compare the result
        difference = np.sum((zer4UpNm-refZer4UpNm)**2)
        
        if difference <= tor:
            result = "true"
        else:
            result = "false"

        self.difference = difference

        return result
    
    def testCase1(self):
        
        case = self.generateTestCase("LSST_NE_SN25", "z11_0.25_", [1.185, 1.185], "exp", "offAxis", 
                                     self.validationDir)        
        result = self.compareCalculation(case, tor)           
        self.assertEqual(result, "true")


    def testCase2(self):

        case = self.generateTestCase("LSST_NE_SN25", "z11_0.25_", [1.185, 1.185], "fft", "offAxis", 
                                     self.validationDir)        
        result = self.compareCalculation(case, tor)
        self.assertEqual(result, "true")

    def testCase3(self):

        case = self.generateTestCase("F1.23_1mm_v61", "z7_0.25_", [0, 0], "fft", "paraxial", 
                                     self.validationDir)        
        result = self.compareCalculation(case, tor)
        self.assertEqual(result, "true")

    def testCase4(self):

        case = self.generateTestCase("LSST_C_SN26", "z7_0.25_", [0, 0], "fft", "onAxis", 
                                     self.validationDir)        
        result = self.compareCalculation(case, tor)
        self.assertEqual(result, "true")

    def testCase5(self):

        case = self.generateTestCase("LSST_C_SN26", "z7_0.25_", [0, 0], "exp", "onAxis", 
                                     self.validationDir)        
        result = self.compareCalculation(case, tor)
        self.assertEqual(result, "true")


if __name__ == "__main__": 

    # Get the path of module
    modulePath = getModulePath()
    
    # Information of test
    InstruFolder = os.path.join(modulePath, "algoData", "cwfs", "instruData")
    AlgoFolderPath = os.path.join(modulePath, "algoData", "cwfs", "algo")
    InstruName = "lsst"
    ImageFolderPath = os.path.join(modulePath, "test", "testImages")

    # Set the tolerance
    tor = 3
    
    # Run the test
    suite = unittest.TestLoader().loadTestsFromTestCase(TestWEP)
    unittest.TextTestRunner(verbosity=0).run(suite)
    
    