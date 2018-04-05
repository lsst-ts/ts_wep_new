import numpy as np
from scipy.spatial.distance import cdist

from lsst.ts.wep.bsc.Filter import Filter

import unittest

class StarData(object):

    def __init__(self, simobjid, ra, decl, lsstMagU, lsstMagG, lsstMagR, lsstMagI, 
                lsstMagZ, lsstMagY):
        
        self.Detector = "" 
      
        self.SimobjID = simobjid

        self.RA = ra
        self.RAInPixel = []
        self.Decl = decl
        self.DeclInPixel = []
        
        self.LSSTMagU = lsstMagU 
        self.LSSTMagG = lsstMagG
        self.LSSTMagR = lsstMagR
        self.LSSTMagI = lsstMagI
        self.LSSTMagZ = lsstMagZ
        self.LSSTMagY = lsstMagY

        # Get the filter type defined in the Filter class
        afilter = Filter()
        self.FilterU = afilter.FilterU
        self.FilterG = afilter.FilterG
        self.FilterR = afilter.FilterR
        self.FilterI = afilter.FilterI
        self.FilterZ = afilter.FilterZ
        self.FilterY = afilter.FilterY
        
    def populateDetector(self, detector):
        """
        
        Populates the detector information for this set of stars.
        
        Arguments:
            detector {[string]} -- The name of the detector.
        """
        self.Detector = detector
        
    def populateRAData(self, raInPixel):
        """
        
        Populates the RA pixel and mm data for this set of stars.
        
        Arguments:
            raInPixel {[float]} -- The ra pixel coordinate each star falls on the detector.
        """
        self.RAInPixel = raInPixel
        
    def populateDeclData(self, declInPixel):
        """
        
        Populates the Decl pixel and mm data for this set of stars.
        
        Arguments:
            declInPixel {[float]} -- The decl pixel coordinate each star falls on the detector.
        """
        self.DeclInPixel = declInPixel

    def checkCandidateStars(self, cameraFilter, lowMagnitude, highMagnitude):
        """
        
        Determine the candidate stars based on the magnitude.
        
        Arguments:
            cameraFilter {[string]} -- Filter type of camera: u, g, r, i, z, y.
            lowMagnitude {[float]} -- Lower boundary of magnitude.
            highMagnitude {[float]} -- Higher boundary of magnitude.
        
        Returns:
            [int] -- index of candidate stars
        """

        indexCandidate = []
        if (self.RA):
            if (cameraFilter == self.FilterU):
                indexCandidate = [index for index in range(len(self.RA)) 
                    if self.LSSTMagU[index] >= lowMagnitude and self.LSSTMagU[index] <= highMagnitude]
            
            elif (cameraFilter == self.FilterG):
                indexCandidate = [index for index in range(len(self.RA)) 
                    if self.LSSTMagG[index] >= lowMagnitude and self.LSSTMagG[index] <= highMagnitude]
            
            elif (cameraFilter == self.FilterR):
                indexCandidate = [index for index in range(len(self.RA)) 
                    if self.LSSTMagR[index] >= lowMagnitude and self.LSSTMagR[index] <= highMagnitude]
                
            elif (cameraFilter == self.FilterI):
                indexCandidate = [index for index in range(len(self.RA)) 
                    if self.LSSTMagI[index] >= lowMagnitude and self.LSSTMagI[index] <= highMagnitude]
                
            elif (cameraFilter == self.FilterZ):
                indexCandidate = [index for index in range(len(self.RA)) 
                    if self.LSSTMagZ[index] >= lowMagnitude and self.LSSTMagZ[index] <= highMagnitude]
                
            elif (cameraFilter == self.FilterY):
                indexCandidate = [index for index in range(len(self.RA)) 
                    if self.LSSTMagY[index] >= lowMagnitude and self.LSSTMagY[index] <= highMagnitude]

        return indexCandidate

    def getNeighboringStar(self, indexCandidate, maxDistance, cameraFilter, maxNeighboringStar):
        """
        
        Get the neighboring stars of candidate stars based on specific max distance.
        
        Arguments:
            indexCandidate {[int]} -- Index of candidate star in "stars" data.
            maxDistance {[float]} -- Maximum distance in pixel.
            cameraFilter {[string]} -- Filter type of camera: u, g, r, i, z, y.

        Returns:
            [NeighboringStar] -- Information of neighboring stars.
        """

        # Calculate the distance in pixel between candidate stars and all stars
        if (indexCandidate):
            allStarXY = np.array([self.RAInPixel, self.DeclInPixel]).transpose()
            candidateStarXY = allStarXY[np.array(indexCandidate),:]
            starDistances = cdist(candidateStarXY,allStarXY)

        # Decide the number of neighboring stars
        neighboringStar = NeighboringStar()
        for ii in range(len(indexCandidate)):
            indexNeighboringStar = np.where(starDistances[ii,:] < maxDistance)[0]
            
            # Delete candidate star itself
            indexNeighboringStar = np.delete(indexNeighboringStar, 
                                             np.where(indexNeighboringStar==indexCandidate[ii]))

            # Remove the candidate star if there is the neighboring star brighter than itself
            if (cameraFilter == self.FilterU):
                magSelf = self.LSSTMagU[indexCandidate[ii]]
                magNeighboringStar = np.array(self.LSSTMagU)[[indexNeighboringStar]]

            elif (cameraFilter == self.FilterG):
                magSelf = self.LSSTMagG[indexCandidate[ii]]
                magNeighboringStar = np.array(self.LSSTMagG)[[indexNeighboringStar]]

            elif (cameraFilter == self.FilterR):
                magSelf = self.LSSTMagR[indexCandidate[ii]]
                magNeighboringStar = np.array(self.LSSTMagR)[[indexNeighboringStar]]

            elif (cameraFilter == self.FilterI):
                magSelf = self.LSSTMagI[indexCandidate[ii]]
                magNeighboringStar = np.array(self.LSSTMagI)[[indexNeighboringStar]]

            elif (cameraFilter == self.FilterZ):
                magSelf = self.LSSTMagZ[indexCandidate[ii]]
                magNeighboringStar = np.array(self.LSSTMagZ)[[indexNeighboringStar]]

            elif (cameraFilter == self.FilterY):
                magSelf = self.LSSTMagY[indexCandidate[ii]]
                magNeighboringStar = np.array(self.LSSTMagY)[[indexNeighboringStar]]

            if ((np.where(magNeighboringStar < magSelf)[0]).size != 0):
                brighterNeighbor = True
            else:
                brighterNeighbor = False 

            # Restrict the maximum number of neighboring stars
            if (indexNeighboringStar.size > maxNeighboringStar):
                highNeighboringStar = True
            else:
                highNeighboringStar = False

            # Record the information of neighboring stars
            if (brighterNeighbor == False and highNeighboringStar == False):
                neighboringStar.addStar(self, indexCandidate[ii], indexNeighboringStar, cameraFilter)

        return neighboringStar

class NeighboringStar(object):

    def __init__(self):

        self.SimobjID = {} 
        self.RaDecl = {}
        self.RaDeclInPixel = {}
        
        self.LSSTMagU = {}
        self.LSSTMagG = {}
        self.LSSTMagR = {}
        self.LSSTMagI = {}
        self.LSSTMagZ = {}
        self.LSSTMagY = {}

    def addStar(self, stars, indexCandidate, indexNeighboringStar, cameraFilter):
        """
        
        Add the information of neightboring stars based on the candidate star.
        
        Arguments:
            stars {[StarData]} -- Star information.
            indexCandidate {[int]} -- Index of candidate stars.
            indexNeighboringStar {[int]} -- Index of neighboring stars of specific candidate star.
            cameraFilter {[string]} -- Filter type of camera: u, g, r, i, z, y.
        """

        self.SimobjID[stars.SimobjID[indexCandidate]] = np.array(stars.SimobjID)[[indexNeighboringStar]].tolist()

        # Collect coordinates and magnitude of candidate and neighboring stars 
        indexStar = np.append(indexNeighboringStar, indexCandidate)
        for index in indexStar:
            self.RaDecl[stars.SimobjID[index]] = (stars.RA[index], stars.Decl[index])
            self.RaDeclInPixel[stars.SimobjID[index]] = (stars.RAInPixel[index], stars.DeclInPixel[index])

            if (cameraFilter == stars.FilterU):
                self.LSSTMagU[stars.SimobjID[index]] = stars.LSSTMagU[index]

            elif (cameraFilter == stars.FilterG):
                self.LSSTMagG[stars.SimobjID[index]] = stars.LSSTMagG[index]

            elif (cameraFilter == stars.FilterR):
                self.LSSTMagR[stars.SimobjID[index]] = stars.LSSTMagR[index]

            elif (cameraFilter == stars.FilterI):
                self.LSSTMagI[stars.SimobjID[index]] = stars.LSSTMagI[index]

            elif (cameraFilter == stars.FilterZ):
                self.LSSTMagZ[stars.SimobjID[index]] = stars.LSSTMagZ[index]

            elif (cameraFilter == stars.FilterY):
                self.LSSTMagY[stars.SimobjID[index]] = stars.LSSTMagY[index]

class StarDataTest(unittest.TestCase):
    """
    Test the function of StarData and NeighboringStar.
    """

    def setUp(self):
        self.stars = StarData([123, 456, 789], [0.1, 0.2, 0.3], [2.1, 2.2, 2.3], [2.0, 3.0, 4.0], 
                              [2.1, 2.1, 4.1], [2.2, 3.2, 4.2], [2.3, 3.3, 4.3], [2.4, 3.4, 4.4], 
                              [2.5, 3.5, 4.5])

    def testStarData(self):
        stars = self.stars
        stars.populateDetector("CCD")

        self.assertEqual(stars.SimobjID, [123, 456, 789])
        self.assertEqual(stars.RA, [0.1, 0.2, 0.3])
        self.assertEqual(stars.Decl, [2.1, 2.2, 2.3])
        self.assertEqual(stars.LSSTMagU, [2.0, 3.0, 4.0])
        self.assertEqual(stars.LSSTMagG, [2.1, 2.1, 4.1])
        self.assertEqual(stars.LSSTMagR, [2.2, 3.2, 4.2])
        self.assertEqual(stars.LSSTMagI, [2.3, 3.3, 4.3])
        self.assertEqual(stars.LSSTMagZ, [2.4, 3.4, 4.4])
        self.assertEqual(stars.LSSTMagY, [2.5, 3.5, 4.5])
        self.assertEqual(stars.Detector,"CCD")

    def testCheckCandidateStars(self):
        stars = self.stars

        indexCandidateU = stars.checkCandidateStars("u", 1.9, 2.1)
        indexCandidateG = stars.checkCandidateStars("g", 0, 5)
        indexCandidateR = stars.checkCandidateStars("r", 0, 1)
        indexCandidateI = stars.checkCandidateStars("i", 2.1, 4.0)
        indexCandidateZ = stars.checkCandidateStars("z", 3.0, 5.0)
        indexCandidateY = stars.checkCandidateStars("y", 1.0, 2.0)

        self.assertEqual(indexCandidateU, [0])
        self.assertEqual(indexCandidateG, [0, 1, 2])
        self.assertEqual(indexCandidateR, [])
        self.assertEqual(indexCandidateI, [0, 1])
        self.assertEqual(indexCandidateZ, [1, 2])
        self.assertEqual(indexCandidateY, [])

    def testGetNeighboringStar(self):
        stars = self.stars
        stars.populateRAData([value*10 for value in stars.RA])
        stars.populateDeclData([value*10 for value in stars.Decl])

        neighboringStarU = stars.getNeighboringStar([0], 3, "u", 99)
        neighboringStarG = stars.getNeighboringStar([0, 1], 3, "g", 99)
        neighboringStarR = stars.getNeighboringStar([0], 1, "r", 99)
        neighboringStarI = stars.getNeighboringStar([], 3, "i", 99)
        neighboringStarZ = stars.getNeighboringStar([0, 1], 2, "z", 1)
        neighboringStarY = stars.getNeighboringStar([1], 2, "y", 1)

        self.assertEqual(stars.RAInPixel, [1, 2, 3])
        self.assertEqual(stars.DeclInPixel, [21, 22, 23])

        self.assertEqual(len(neighboringStarU.SimobjID[123]), 2)
        self.assertEqual(len(neighboringStarU.RaDecl), 3)
        self.assertEqual(len(neighboringStarG.SimobjID), 2)
        self.assertEqual(len(neighboringStarR.SimobjID[123]), 0)
        self.assertEqual(neighboringStarI.SimobjID, {})
        self.assertEqual(len(neighboringStarZ.SimobjID), 1)
        self.assertEqual(neighboringStarY.SimobjID, {})

if __name__ == "__main__":

    # Do the unit test
    unittest.main()