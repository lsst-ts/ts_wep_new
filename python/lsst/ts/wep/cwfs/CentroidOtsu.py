from skimage.filters import threshold_otsu

from lsst.ts.wep.cwfs.CentroidDefault import CentroidDefault


class CentroidOtsu(CentroidDefault):

    def __init__(self):
        """CentroidDefault child class to get the centroid of donut by the
        Otsu's method."""

        # Number of bins in the histogram
        self.numOfBins = 256

    def getImgBinary(self, imgDonut):
        """Get the binary image.

        Parameters
        ----------
        imgDonut : numpy.ndarray
            Donut image to do the analysis.

        Returns
        -------
        numpy.ndarray [int]
            Binary image of donut.
        """

        threshold = threshold_otsu(imgDonut, nbins=self.numOfBins)
        imgBinary = (imgDonut > threshold).astype(int)

        return imgBinary
