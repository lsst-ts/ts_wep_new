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
from copy import deepcopy
from lsst.ts.wep.cwfs.CentroidDefault import CentroidDefault
from lsst.ts.wep.cwfs.CentroidRandomWalk import CentroidRandomWalk
from scipy.signal import correlate
from sklearn.cluster import KMeans


class CentroidConvolveTemplate(CentroidDefault):

    """CentroidDefault child class to get the centroid of donut by
    convolution."""

    def getCenterAndR(
        self, imgDonut, templateDonut=None, nDonuts=1, peakThreshold=0.95
    ):
        """Get the centroid data and effective weighting radius.

        Parameters
        ----------
        imgDonut : numpy.ndarray
            Donut image.
        templateDonut : None or numpy.ndarray
            Template image for a single donut. If set to None
            then the image will be convolved with itself. (The Default is None)
        nDonuts : int
            Number of donuts in the image. (The Default is 1)
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

        if templateDonut is None:
            templateDonut = deepcopy(imgDonut)

        self._centRandomWalk = CentroidRandomWalk()
        imgBinary = self._centRandomWalk.getImgBinary(imgDonut)
        templateBinary = self._centRandomWalk.getImgBinary(templateDonut)

        return self.getCenterAndRfromImgBinary(
            imgBinary,
            templateBinary=templateBinary,
            nDonuts=nDonuts,
            peakThreshold=peakThreshold,
        )

    def getCenterAndRfromImgBinary(
        self, imgBinary, templateBinary=None, nDonuts=1, peakThreshold=0.95
    ):
        """Get the centroid data and effective weighting radius.

        Parameters
        ----------
        imgBinary : numpy.ndarray
            Binary image of donut.
        templateBinary : None or numpy.ndarray
            Binary image of template for a single donut. If set to None
            then the image will be convolved with itself. (The Default is None)
        nDonuts : int
            Number of donuts in the image. (The Default is 1)
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

        if templateBinary is None:
            templateBinary = deepcopy(imgBinary)

        x, y = self.getCenterFromTemplateConv(
            imgBinary,
            templateImgBinary=templateBinary,
            nDonuts=nDonuts,
            peakThreshold=peakThreshold,
        )
        radius = np.sqrt(np.sum(templateBinary) / np.pi)

        return x, y, radius

    def getCenterFromTemplateConv(
        self, imageBinary, templateImgBinary=None, nDonuts=1, peakThreshold=0.95
    ):
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
        templateImgBinary: None or numpy.ndarray
            Binary image of template donut. If set to None then the image
            is convolved with itself. (The default is None)
        nDonuts: int
            Number of donuts there should be in the binary image. (The default
            is 1)
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

        if templateImgBinary is None:
            templateImgBinary = deepcopy(imageBinary)

        # We set the mode to be "same" because we need to return the same
        # size image to the code.
        tempConvolve = correlate(imageBinary, templateImgBinary, mode="same")

        # Then we rank the pixel values keeping only those above
        # some fraction of the highest value.
        rankedConvolve = np.argsort(tempConvolve.flatten())[::-1]
        cutoff = len(
            np.where(tempConvolve.flatten() > peakThreshold * np.max(tempConvolve))[0]
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
