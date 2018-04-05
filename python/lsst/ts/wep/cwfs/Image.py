#!/usr/bin/env python
##
# @package cwfs
# @file image.py
##
# @authors: Bo Xin & Chuck Claver
# @       Large Synoptic Survey Telescope
##
# Refactored by Te-Wei Tsai at June, 2017
##

# getCenterAndR() is partly based on the EF wavefront sensing software
# by Laplacian Optics

import os, unittest
import numpy as np
from astropy.io import fits

from scipy.stats import entropy
from scipy.ndimage.measurements import center_of_mass
from lsst.ts.wep.Utility import getModulePath

class Image(object):
    
    def __init__(self):
        """
        
        Image class for wavefront estimation.
        """

        # Image parameters
        self.image = None
        self.name = None

    def setImg(self, image=None, imageFile=None):
        """
        
        Set the wavefront image.
                
        Keyword Arguments:
            image {[ndarray]} -- Array of image. (default: {None})
            imageFile {[str]} -- Path of image file. (default: {None})
        """

        # Read the file if there is no input image
        if (image is not None):
            self.image = image
        else:
            if (imageFile is not None):
                self.image = self.__readImgFile(imageFile)
                self.name = imageFile

    def __readImgFile(self, imageFile):
        """
        
        Read the donut image.
        
        Arguments:
            imageFile {[path]} -- Path of image file.
        
        Returns:
            [ndarray] -- Image data.
        
        Raises:
            IOError -- IO error if the file type is not ".txt" or ".fits".
        """
 
        image = None

        # Check the format of image 
        if (os.path.isfile(imageFile)):
            if (imageFile.endswith((".fits", ".fits.gz"))):
                image = fits.getdata(imageFile)
            else:
                image = np.loadtxt(imageFile)
                # This assumes this "txt" file is in the format
                # I[0,0]   I[0,1]
                # I[1,0]   I[1,1]
                image = image[::-1, :]
    
        if (image is None):
            raise IOError("Unrecognised file type for %s" % imageFile)
        
        return image

    def getCenterAndR_ef(self, image=None, randNumFilePath=None, histogram_len=256, checkEntropy=False, 
                            entroThres=3.5, debugLevel=0):
        """
        
        Centroid finding code based on northcott_ef_bundle/ef/ef/efimageFunc.cc.
        This is the modified version of getCenterAndR_ef.m by Bo Xin at 6/25/14.
        This is further modified by Te-Wei Tsai for the deblending use at 6/8/17.
        
        Keyword Arguments:
            image {[ndarray]} -- Image to do the analysis (default: {None}).
            randNumFilePath {[str]} -- Random table file path. If not None, read this table instead of 
                                      using numpy random number function (default: {None}).
            histogram_len {[int]} -- Nuber of bins in histogram (default: {256}).
            checkEntropy {[bool]} -- Check the entropy of figure intensity to decide the image quality 
                                    (default: {False}).
            entroThres {[float]} -- Threshold of entropy check (default: {1.5}).

            debugLevel {[int]} -- Show the information under the running. If the value is higher, 
                                  the information shows more. It can be 0, 1, 2, or 3. (default: {0})

        Returns:
            [float] -- Centroid x.
            [float] -- Centroid y.
            [float] -- Effective weighting radius.
            [ndarray] -- Binary image of bright star.
        """

        # Parameters of circle
        realcx = []
        realcy = []
        realR = []

        # Binary image of bright star
        imgBinary = []

        # Parameters to decide the signal of bright star
        slide = int(0.1*histogram_len)       
        stepsize = int(0.06*histogram_len)
        nwalk = int(1.56*histogram_len)

        # Copy the image
        if (image is not None):
            tempImage = image
        else:
            tempImage = self.image.copy()

        # Reshape the image to 1D array
        array1d = tempImage.flatten()
    
        # Generate the histogram of intensity
        phist, cen = np.histogram(array1d, bins=histogram_len)

        # Check the entropy of intensity distribution
        if (checkEntropy):

            # Square the distribution to magnify the difference, and calculate the entropy 
            figureEntropy = entropy(phist**2)
            
            if (figureEntropy > entroThres) or (figureEntropy == 0):
                print("Entropy is %f > %f. Can not differentiate the star signal." % (figureEntropy, 
                                                                                        entroThres))
                return realcx, realcy, realR, imgBinary

        # Parameters for random walk search
        start = int(histogram_len/2.1)
        end = slide + 25  # Go back 
        startidx = range(start, end, -15)

        foundvalley = False

        # Use a fixed rand table for the evaluation
        if (randNumFilePath is not None):
            iRand = 0
            myRand = np.loadtxt(randNumFilePath)
            myRand = np.tile(np.reshape(myRand, (1000, 1)), (10, 1))
            
        for istartPoint in range(len(startidx)):
            minind = startidx[istartPoint]

            # Check the condition of start index
            if ((minind <= 0) or (max(phist[minind - 1:]) == 0)):
                continue
            minval = phist[minind - 1]

            # Do the randon walk search
            for ii in range(nwalk + 1):

                # if (minind <= slide):
                if (minind >= slide):
                    
                    # Find the index of bin that the count is not zero
                    while (minval == 0):
                        minind = minind - 1
                        minval = phist[int(minind - 1)]

                    # Generate the thermal fluctuation based on the random table
                    # to give a random walk/ step with a random thermal fluctuation.
                    if (randNumFilePath is not None):
                        ind = np.round(stepsize * (2 * myRand[iRand, 0] - 1))
                        iRand += 1
                        thermal = 1 + 0.5*myRand[iRand, 0] * \
                            np.exp(1.0 * ii / (nwalk * 0.3))
                        iRand += 1
                    else:
                        ind = np.round(stepsize * (2 * np.random.rand() - 1))
                        thermal = 1 + 0.5*np.random.rand()*np.exp(1.0*ii/(nwalk*0.3))
    
                    # Check the index of bin is whithin the range of histogram
                    if ((minind + ind < 1) or (minind + ind > (histogram_len))):
                        continue

                    # Look for the minimum point
                    if (phist[int(minind + ind - 1)] < (minval * thermal)):

                        # Add the panality to go to the high intensity position
                        if (ind>0):
                            ind = int(ind/3)
                        else:
                            ind = int(ind/2)

                        # Update the value of minind
                        minval = phist[int(minind + ind - 1)]
                        minind = minind + ind

                else:
                    break

            # Find the signal of bright star in histogram
            if (minind >= slide):
                foundvalley = True
                break

        # Try to close the second peak
        while (minind >= slide) and (foundvalley == True):
            if np.abs(phist[int(minind-5)]-phist[int(minind)]) < 4*np.median(phist[len(phist)-20:]):
                minind = minind - 1
            else:
                print("Stop search. Final minind in histogram is %f." % minind)
                break

        # If no valley (signal) is found for the noise, use the value at start index 
        # of histogram to be the threshold.
        if (not foundvalley):
            minind = start
            # Show the value of minind
            if (debugLevel >=3):
                print("Valley is not found. Use minind = %f." % minind)

        # Get the threshold value of bright star
        pval = cen[int(minind)]

        # Get the binary image
        imgBinary = tempImage.copy()
        imgBinary[tempImage > max(0, pval - 1e-8)] = 1
        imgBinary[tempImage < pval] = 0

        # Calculate the weighting radius
        realR = np.sqrt(np.sum(imgBinary) / np.pi)
    
        # Calculate the center of mass
        realcy, realcx = center_of_mass(imgBinary)

        return realcx, realcy, realR, imgBinary

    def updateImage(self, image):
        """
        
        Update the image of donut.
        
        Arguments:
            image {[float]} -- Donut image.
        """

        # Update the image
        if (self.image is not None):
            self.image = image
        else:
            print("The attribute:image is None. Use setImg() instead.")

    def getSNR(self):
        """
        
        Get the signal to noise ratio of donut.
        
        Returns:
            [float] -- Signal to noise ratio.
        """

        # Get the signal binary image
        realcx, realcy, realR, imgBinary = self.getCenterAndR_ef()

        # Get the background binary img
        bgBinary = 1-imgBinary

        # Get the donut image signal and calculate the intensity
        signal = self.image*imgBinary
        signal = np.mean(signal[signal!=0])

        # Get the backgrond signal
        bg = self.image*bgBinary
        bg = bg[bg!=0]

        # Calculate the noise
        noise = np.std(bg-np.mean(bg))

        # Calculate SNR
        snr = signal/noise

        return snr

class ImageTest(unittest.TestCase):
    """
    Test functions in Image.
    """
    def setUp(self):

        # Get the path of module
        modulePath = getModulePath()

        # Define the algorithm folder 
        algoFolderPath = os.path.join(modulePath, "algoData", "cwfs", "algo")
                
        # Define the image folder and image names
        # Image data -- Don't know the final image format. 
        # It is noted that image.readFile inuts is based on the txt file  
        imageFolderPath = os.path.join(modulePath, "test", "testImages", "LSST_NE_SN25")
        imgName = "z11_0.25_intra.txt"
                
        # Image files Path  
        imgFile = os.path.join(imageFolderPath, imgName)

        # There is the difference between intra and extra images
        self.img = Image()
        self.img.setImg(imageFile=imgFile)

    def testZeroImg(self):

        # Creat a zero image
        zeroImg = Image()
        zeroImg.setImg(image=np.zeros([4, 4]))
        self.assertEqual(np.sum(zeroImg.image), 0)

        # Update Image
        zeroImg.updateImage(np.ones([4, 4]))
        self.assertEqual(np.sum(zeroImg.image), 16)

        realcx, realcy, realR, imgBinary = zeroImg.getCenterAndR_ef(randNumFilePath=None, 
                                                                    checkEntropy=True)
        self.assertEqual(realcx, [])

        # update to the random image
        zeroImg.updateImage(np.random.rand(100,100))
        realcx, realcy, realR, imgBinary = zeroImg.getCenterAndR_ef(randNumFilePath=None, 
                                                                    checkEntropy=True)
        self.assertEqual(realcx, [])

    def testImg(self):

        realcx, realcy, realR, imgBinary = self.img.getCenterAndR_ef(randNumFilePath=None, 
                                                                     checkEntropy=True)
        self.assertEqual(int(realcx), 61)
        self.assertEqual(int(realcy), 61)
        self.assertGreater(int(realR), 35)

        # Calculate the S/N
        # Add the noise to the image
        noisedImg = self.img.image + np.random.random(self.img.image.shape)*0.1
        self.img.setImg(image=noisedImg)
        snr = self.img.getSNR()

        self.assertGreater(snr, 15)

if __name__ == "__main__":

    # Do the unit test
    unittest.main()