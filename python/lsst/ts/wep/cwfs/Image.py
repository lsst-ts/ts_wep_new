import os
import numpy as np
from astropy.io import fits

from scipy.stats import entropy
from scipy.ndimage.measurements import center_of_mass


class Image(object):

    def __init__(self):
        """Image class for wavefront estimation."""

        # Image parameters
        self.image = None
        self.name = None

    def getImg(self):
        """Get the image.

        Returns
        -------
        numpy.ndarray
            Get the image.
        """

        return self.image

    def setImg(self, image=None, imageFile=None):
        """Set the wavefront image.

        Parameters
        ----------
        image : numpy.ndarray, optional
            Array of image. (the default is None.)
        imageFile : str, optional
            Path of image file. (the default is None.)
        """

        # Read the file if there is no input image
        if (image is not None):
            self.image = image
        else:
            if (imageFile is not None):
                self.image = self._readImgFile(imageFile)
                self.name = imageFile

    def _readImgFile(self, imageFile):
        """Read the donut image.

        Parameters
        ----------
        imageFile : str
            Path of image file.

        Returns
        -------
        numpy.ndarray
            Image data.

        Raises
        ------
        IOError
            IO error if the file type is not ".txt" or ".fits".
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

    def getCenterAndR_ef(self, image=None, randNumFilePath=None,
                         histogram_len=256, checkEntropy=False,
                         entroThres=3.5, debugLevel=0):
        """Get the centroid data by the random walk model.

        Parameters
        ----------
        image : numpy.ndarray, optional
            Image to do the analysis. (the default is None.)
        randNumFilePath : str, optional
            Random table file path. If not None, read this table instead of
            using numpy random number function. (the default is None.)
        histogram_len : int, optional
            Nuber of bins in histogram. (the default is 256.)
        checkEntropy : bool, optional
            Check the entropy of figure intensity to decide the image quality.
            (the default is False.)
        entroThres : float, optional
            Threshold of entropy check. (the default is 3.5.)
        debugLevel : int, optional
            Show the information under the running. If the value is higher,
            the information shows more. It can be 0, 1, 2, or 3. (the default
            is 0.)

        Returns
        -------
        float
            Centroid x.
        float
            Centroid y.
        float
            Effective weighting radius.
        numpy.ndarray[int]
            Binary image of bright star.
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

            # Square the distribution to magnify the difference, and calculate
            # the entropy
            figureEntropy = entropy(phist**2)

            if (figureEntropy > entroThres) or (figureEntropy == 0):
                print("Entropy is %f > %f. Can not differentiate the star signal." % (
                    figureEntropy, entroThres))
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

                    # Generate the thermal fluctuation based on the random
                    # table to give a random walk/ step with a random thermal
                    # fluctuation.
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
                        if (ind > 0):
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
        while (minind >= slide) and (foundvalley is True):
            if np.abs(phist[int(minind-5)]-phist[int(minind)]) < 4*np.median(phist[len(phist)-20:]):
                minind = minind - 1
            else:
                print("Stop search. Final minind in histogram is %f." % minind)
                break

        # If no valley (signal) is found for the noise, use the value at start
        # index of histogram to be the threshold.
        if (not foundvalley):
            minind = start
            # Show the value of minind
            if (debugLevel >= 3):
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
        """Update the image of donut.

        Parameters
        ----------
        image : numpy.ndarray
            Donut image.
        """

        # Update the image
        if (self.image is not None):
            self.image = image
        else:
            print("The attribute:image is None. Use setImg() instead.")

    def getSNR(self):
        """Get the signal to noise ratio of donut.

        Returns
        -------
        float
            Signal to noise ratio.
        """

        # Get the signal binary image
        realcx, realcy, realR, imgBinary = self.getCenterAndR_ef()

        # Get the background binary img
        bgBinary = 1-imgBinary

        # Get the donut image signal and calculate the intensity
        signal = self.image*imgBinary
        signal = np.mean(signal[signal != 0])

        # Get the backgrond signal
        bg = self.image*bgBinary
        bg = bg[bg != 0]

        # Calculate the noise
        noise = np.std(bg-np.mean(bg))

        # Calculate SNR
        snr = signal/noise

        return snr


if __name__ == "__main__":
    pass
