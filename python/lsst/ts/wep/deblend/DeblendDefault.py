import numpy as np

from lsst.ts.wep.Utility import CentroidFindType
from lsst.ts.wep.cwfs.CentroidFindFactory import CentroidFindFactory


class DeblendDefault(object):
    """Default deblend class."""

    def generateMultiDonut(self, template, spaceCoef, magRatio, theta):
        """Gemerate multiple donut images.

        Only one neightboring star will be generated for test, which is the
        baseline of LSST.

        Parameters
        ----------
        template : numpy.ndarray
            Template donut image.
        spaceCoef : float
            Spacing coefficient to decide the distance between donuts.
        magRatio : float
            Magnitude ratio of new donut compared with the original one.
        theta : float
            Theta angle of generated neighboring star in degree.

        Returns
        -------
        numpy.ndarray
            Image of donuts.
        numpy.ndarray
            Image of bright star.
        numpy.ndarray
            Image of neighboring star.
        float
            Position x of neighboring star.
        float
            Position y of neighboring star.

        Raises
        ------
        ValueError
            spaceCoef should be greater than zero.
        ValueError
            magRatio should be postive and less than 1.
        """

        # Check the inputs
        if spaceCoef <= 0:
            raise ValueError("spaceCoef should be greater than zero.")
        elif (magRatio < 0) or (magRatio > 1):
            raise ValueError("magRatio should be postive and less than 1.")

        # Get the center and radius of self-donut
        centroidFind = CentroidFindFactory.createCentroidFind(
            CentroidFindType.Otsu)
        imgBinary = centroidFind.getImgBinary(template)
        selfX, selfY, selfR = centroidFind.getCenterAndRfromImgBinary(
            imgBinary)

        # Get the position of new donut based on spaceCoef and theta
        thetaInRad = np.deg2rad(theta)
        newX = selfX + spaceCoef * selfR * np.cos(thetaInRad)
        newY = selfY + spaceCoef * selfR * np.sin(thetaInRad)

        # Calculate the frame size and shift the center of donuts
        lengthX = max(selfX, newX) - min(selfX, newX) + 5 * selfR
        lengthY = max(selfY, newY) - min(selfY, newY) + 5 * selfR
        length = int(max(lengthX, lengthY))

        # Enforce the length to be even for the symmetry
        if (length % 2 == 1):
            length += 1

        shiftX = length / 2 - (selfX + newX) / 2
        shiftY = length / 2 - (selfY + newY) / 2

        # Get the new coordinate
        selfX += shiftX
        selfY += shiftY

        newX += shiftX
        newY += shiftY

        # Generate image of multiple donuts
        imageMain = np.zeros([length, length])
        imageNeighbor = np.zeros([length, length])

        # Get the shifted main donut image
        m, n = template.shape
        imageMain[int(selfY-m/2):int(selfY+m/2), int(selfX-n/2):int(selfX+n/2)] += template

        # Get the shifted neighboring donut image
        imageNeighbor[int(newY-m/2):int(newY+m/2), int(newX-n/2):int(newX+n/2)] += magRatio * template

        # Get the synthesized multi-donut image
        image = imageMain + imageNeighbor

        return image, imageMain, imageNeighbor, newX, newY

    def deblendDonut(self, imgToDeblend, iniGuessXY, **kwargs):
        """Deblend the donut image.

        Parameters
        ----------
        imgToDeblend : numpy.ndarray
            Image to deblend.
        iniGuessXY : list[tuple]
            The list contains the initial guess of (x, y) positions of
            neighboring stars as [star 1, star 2, etc.].
        **kwargs : dict[str, any]
            Dictionary of input argument: new value for that input argument.

        Returns
        -------
        numpy.ndarray
            Deblended donut image.
        float
            Position x of donut in pixel.
        float
            Position y of donut in pixel.

        Raises
        ------
        NotImplementedError
            Child class should implement this.
        """
        raise NotImplementedError("Child class should implement this.")


if __name__ == "__main__":
    pass
