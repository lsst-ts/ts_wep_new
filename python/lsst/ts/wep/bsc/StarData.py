import numpy as np
from scipy.spatial.distance import cdist

from lsst.ts.wep.bsc.NeighboringStar import NeighboringStar
from lsst.ts.wep.bsc.Filter import Filter


class StarData(object):

    def __init__(self, simobjid, ra, decl, lsstMagU, lsstMagG, lsstMagR,
                 lsstMagI, lsstMagZ, lsstMagY):
        
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


if __name__ == "__main__":
    pass
