import numpy as np

from lsst.ts.wep.cwfs.CentroidRandomWalk import CentroidRandomWalk
from scipy.signal import correlate
from sklearn.cluster import KMeans


class CentroidConvolveTemplate(CentroidRandomWalk):

    """CentroidDefault child class to get the centroid of donut by
    convolution."""

    def getCenterFromTemplateConv(self, imageBinary, templateImgBinary,
                                  n_donuts):
        """Get the centers of the donuts in the binary image.

        Parameters
        ----------
        imageBinary: numpy.ndarray
            Binary image of postage stamp

        templateImgBinary: numpy.ndarray
            Binary image of template donut

        Returns
        -------
        nx, ny: float
            x, y pixel coordinates for donut centroid
        """

        temp_convolve = correlate(imageBinary, templateImgBinary, mode='same', method='fft')
        ranked_convolve = np.argsort(temp_convolve.flatten())[::-1]
        cutoff = len(np.where(temp_convolve.flatten() > 0.95*np.max(temp_convolve))[0])
        ranked_convolve = ranked_convolve[:cutoff]
        nx, ny = np.unravel_index(ranked_convolve, np.shape(imageBinary))

        kmeans = KMeans(n_clusters=n_donuts).fit(np.array([nx, ny]).T)
        labels = kmeans.labels_

        cent_x = []
        cent_y = []

        for label_num in range(n_donuts):
            nx_label, ny_label = np.unravel_index(ranked_convolve[labels == label_num][0],
                                                  np.shape(imageBinary))
            cent_x.append(nx_label)
            cent_y.append(ny_label)

        return cent_x, cent_y
