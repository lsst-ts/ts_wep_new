# This file is part of ts_wep.
#
# Developed for the LSST Telescope and Site Systems.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np

from lsst.ts.wep.cwfs.CentroidDefault import CentroidDefault
from lsst.ts.wep.cwfs.CentroidRandomWalk import CentroidRandomWalk
from scipy.signal import correlate
from sklearn.cluster import KMeans


class CentroidConvolveTemplate(CentroidDefault):

    """CentroidDefault child class to get the centroid of donut by
    convolution."""

    def getCenterAndR(self, imgDonut, templateDonut, nDonuts,
                      peakThreshold=0.95):
        """Get the centroid data and effective weighting radius.

        Parameters
        ----------
        imgDonut : numpy.ndarray
            Donut image.
        templateDonut : numpy.ndarray
            Template image for a single donut.
        nDonuts : int
            Number of donuts in the image
        peakThreshold : float
            This value is a specifies a number between 0 and 1 that is
            the fraction of the highest pixel value in the convolved image.
            The code then sets all pixels with a value below this to 0 before
            running the K-means algorithm to find peaks that represent possible
            donut locations. (The default is 0.95)

        Returns
        -------
        float
            Centroid x.
        float
            Centroid y.
        float
            Effective weighting radius.
        """

        self._centRandomWalk = CentroidRandomWalk()
        imgBinary = self._centRandomWalk.getImgBinary(imgDonut)
        templateBinary = self._centRandomWalk.getImgBinary(templateDonut)

        return self.getCenterAndRfromImgBinary(imgBinary, templateBinary,
                                               nDonuts,
                                               peakThreshold=peakThreshold)

    def getCenterAndRfromImgBinary(self, imgBinary, templateBinary, nDonuts,
                                   peakThreshold=0.95):
        """Get the centroid data and effective weighting radius.

        Parameters
        ----------
        imgBinary : numpy.ndarray
            Binary image of donut.
        templateBinary : numpy.ndarray
            Binary image of template for a single donut.
        nDonuts : int
            Number of donuts in the image
        peakThreshold : float
            This value is a specifies a number between 0 and 1 that is
            the fraction of the highest pixel value in the convolved image.
            The code then sets all pixels with a value below this to 0 before
            running the K-means algorithm to find peaks that represent possible
            donut locations. (The default is 0.95)

        Returns
        -------
        float
            Centroid x.
        float
            Centroid y.
        float
            Effective weighting radius.
        """

        x, y = self.getCenterFromTemplateConv(imgBinary, templateBinary,
                                              nDonuts,
                                              peakThreshold=peakThreshold)
        radius = np.sqrt(np.sum(templateBinary) / np.pi)

        return x, y, radius

    def getCenterFromTemplateConv(self, imageBinary, templateImgBinary,
                                  nDonuts, peakThreshold=0.95):
        """
        Get the centers of the donuts by convolving a binary template image
        with the binary image of the donut or donuts.

        Peaks will appear as bright spots in the convolved image. Since we
        use binary images the brightness of the stars does not matter and
        the peaks of any stars in the image should have about the same
        brightness if the template is correct.

        Parameters
        ----------
        imageBinary: numpy.ndarray
            Binary image of postage stamp.
        templateImgBinary: numpy.ndarray
            Binary image of template donut.
        nDonuts: int
            Number of donuts there should be in the binary image.
        peakThreshold: float
            This value is a specifies a number between 0 and 1 that is
            the fraction of the highest pixel value in the convolved image.
            The code then sets all pixels with a value below this to 0 before
            running the K-means algorithm to find peaks that represent possible
            donut locations. (The default is 0.95)

        Returns
        -------
        list
            X pixel coordinates for donut centroid.
        list
            Y pixel coordinates for donut centroid.
        """

        # We set the mode to be "same" because we need to return the same
        # size image to the code.
        tempConvolve = correlate(
            imageBinary, templateImgBinary, mode="same"
        )

        # Then we rank the pixel values keeping only those above
        # some fraction of the highest value.
        rankedConvolve = np.argsort(tempConvolve.flatten())[::-1]
        cutoff = len(
            np.where(
                tempConvolve.flatten() > peakThreshold * np.max(tempConvolve)
            )[0]
        )
        rankedConvolve = rankedConvolve[:cutoff]
        nx, ny = np.unravel_index(rankedConvolve, np.shape(imageBinary))

        # Then to find peaks in the image we use K-Means with the
        # specified number of donuts
        kmeans = KMeans(n_clusters=nDonuts).fit(np.array([nx, ny]).T)
        labels = kmeans.labels_

        centX = []
        centY = []

        # Then in each cluster we take the brightest pixel as the centroid
        for labelNum in range(nDonuts):
            nxLabel, nyLabel = np.unravel_index(
                rankedConvolve[labels == labelNum][0], np.shape(imageBinary)
            )
            centX.append(nxLabel)
            centY.append(nyLabel)

        return centX, centY
