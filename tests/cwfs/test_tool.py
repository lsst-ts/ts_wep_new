import numpy as np
from scipy.integrate import nquad
import unittest

from lsst.ts.wep.cwfs.Tool import ZernikeAnnularEval, padArray, extractArray


class TestTool(unittest.TestCase):
    """Test the fuctions in Tool."""

    def setUp(self):

        # Create the mesh of x, y-coordinate
        point = 400
        ratio = 0.9
        xx, yy = self._genGridXy(point, ratio)

        # Define the attributes
        self.xx = xx
        self.yy = yy

    def _genGridXy(self, point, ratio):

        yy, xx = np.mgrid[-(point / 2 - 0.5):(point / 2 + 0.5),
                          -(point / 2 - 0.5):(point / 2 + 0.5)]

        xx = xx / (point * ratio / 2)
        yy = yy / (point * ratio / 2)

        return xx, yy

    def testZernikeAnnularEval(self):

        # Obscuration
        e = 0.61

        # Calculate the radius
        dd = np.sqrt(self.xx**2 + self.yy**2)

        # Define the invalid range
        idx = (dd > 1) | (dd < e)

        # Create the Zernike terms
        Z = np.zeros(22)

        # Generate the map of Z12
        Z[11] = 1

        # Calculate the map of Zernike polynomial
        Zmap = ZernikeAnnularEval(Z, self.xx, self.yy, e)
        Zmap[idx] = np.nan

        # Put the elements to be 0 in the invalid region
        Zmap[np.isnan(Zmap)] = 0

        # Check the normalization for Z1 - Z28
        e = 0.61
        ansValue = np.pi*(1-e**2)
        for ii in range(28):
            Z = np.zeros(28)
            Z[ii] = 1

            normalization = nquad(self._genNormalizedFunc,
                                  [[e, 1], [0, 2*np.pi]], args=(Z, e))[0]

            self.assertAlmostEqual(normalization, ansValue)

        # Check the orthogonality for Z1 - Z28
        for jj in range(28):
            Z1 = np.zeros(28)
            Z1[jj] = 1
            for ii in range(28):
                if (ii != jj):
                    Z2 = np.zeros(28)
                    Z2[ii] = 1

                    orthogonality = nquad(self._genOrthogonalFunc,
                                          [[e, 1], [0, 2*np.pi]],
                                          args=(Z1, Z2, e))[0]

                    self.assertAlmostEqual(orthogonality, 0)

    def _genNormalizedFunc(self, r, theta, Z, e):

        func = r * ZernikeAnnularEval(Z, r*np.cos(theta), r*np.sin(theta), e)**2

        return func

    def _genOrthogonalFunc(self, r, theta, Z1, Z2, e):

        func = r * ZernikeAnnularEval(Z1, r*np.cos(theta), r*np.sin(theta), e) * \
            ZernikeAnnularEval(Z2, r*np.cos(theta), r*np.sin(theta), e)

        return func

    def testPadArray(self):

        imgDim = 10
        padPixelSize = 20

        img, imgPadded = self._padRandomImg(imgDim, padPixelSize)

        self.assertEqual(imgPadded.shape[0], imgDim + padPixelSize)

    def _padRandomImg(self, imgDim, padPixelSize):

        img = np.random.rand(imgDim, imgDim)
        imgPadded = padArray(img, imgDim + padPixelSize)

        return img, imgPadded

    def testExtractArray(self):

        imgDim = 10
        padPixelSize = 20
        img, imgPadded = self._padRandomImg(imgDim, padPixelSize)

        imgExtracted = extractArray(imgPadded, imgDim)

        self.assertEqual(imgExtracted.shape[0], imgDim)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
