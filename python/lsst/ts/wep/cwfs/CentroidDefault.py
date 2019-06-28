import numpy as np
from scipy.ndimage import center_of_mass


class CentroidDefault(object):
    """Default Centroid class."""

    def getCenterAndR(self, imgDonut):
        """Get the centroid data and effective weighting radius.

        Parameters
        ----------
        imgDonut : numpy.ndarray
            Donut image.

        Returns
        -------
        float
            Centroid x.
        float
            Centroid y.
        float
            Effective weighting radius.
        """

        imgBinary = self.getImgBinary(imgDonut)

        return self.getCenterAndRfromImgBinary(imgBinary)

    def getCenterAndRfromImgBinary(self, imgBinary):
        """Get the centroid data and effective weighting radius from the binary
        image.

        Parameters
        ----------
        imgBinary : numpy.ndarray [int]
            Binary image of donut.

        Returns
        -------
        float
            Centroid x.
        float
            Centroid y.
        float
            Effective weighting radius.
        """

        y, x = center_of_mass(imgBinary)
        radius = np.sqrt(np.sum(imgBinary) / np.pi)

        return x, y, radius

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

        Raises
        ------
        NotImplementedError
            Child class should implement this.
        """

        raise NotImplementedError("Child class should implement this.")


if __name__ == "__main__":
    pass
