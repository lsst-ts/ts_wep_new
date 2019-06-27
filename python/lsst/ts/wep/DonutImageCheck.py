import numpy as np
from scipy.stats import entropy


class DonutImageCheck(object):

    def __init__(self, numOfBins=256, entroThres=3.5):
        """Donut image check class to judge the donut image is effective or
        not.

        Parameters
        ----------
        numOfBins : int, optional
            Number of bins in the histogram. (the default is 256.)
        entroThres : float, optional
            Threshold of entropy (the default is 3.5.)
        """

        # Number of bins in the histogram
        self.numOfBins = int(numOfBins)

        # Threshold of entropy
        self.entroThres = entroThres

    def isEffDonut(self, donutImg):
        """Is effective donut image or not.

        Parameters
        ----------
        donutImg : numpy.ndarray
            Donut image.

        Returns
        -------
        bool
            True if the donut image is effective.
        """

        array1d = donutImg.flatten()
        hist = np.histogram(array1d, bins=self.numOfBins)[0]

        # Square the distribution to magnify the difference in entropy
        imgEntropy = entropy(hist**2)
        if (imgEntropy < self.entroThres) and (imgEntropy != 0):
            return True
        else:
            return False


if __name__ == "__main__":
    pass
