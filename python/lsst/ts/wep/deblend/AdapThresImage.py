import os, unittest
import numpy as np

from skimage.filters import threshold_local
from scipy.ndimage.measurements import center_of_mass

from lsst.ts.wep.cwfs.Image import Image
from lsst.ts.wep.Utility import getModulePath

class AdapThresImage(Image):

    def generateMultiDonut(self, spaceCoef, magRatio, theta):
        """
        
        Gemerate multiple donut images. Only one neightboring star will be generated for test,
        which is the baseline of LSST.
        
        Arguments:
            spaceCoef {[float]} -- Spacing coefficient to decide the distance between donuts.
            magRatio {[float]} -- Magnitude ratio of new donut compared with the original one.
            theta {[float]} -- Theta angle of generated neighboring star.
        
        Returns:
            [float] -- Image of donuts.
            [float] -- Neighboring donut (x, y) position.
        """

        # Check the inputs
        if (spaceCoef < 0):
            print("spaceCoef should be greater than zero.")
            return -1
        elif (magRatio < 0 or magRatio > 1):
            print("magRatio should be postive and less than 1.")
            return -1

        # Get the center and radius of self-donut
        selfX, selfY, selfR, imgBinary = self.getCenterAndR_ef(checkEntropy = True)

        # Get the position of new donut based on spaceCoef and theta
        newX = selfX + spaceCoef*selfR*np.cos(theta)
        newY = selfY + spaceCoef*selfR*np.sin(theta)

        # Calculate the frame size and shift the center of donuts
        lengthX = max(selfX, newX) - min(selfX, newX) + 5*selfR
        lengthY = max(selfY, newY) - min(selfY, newY) + 5*selfR
        length = int(max(lengthX, lengthY))

        # Enforce the length to be even for the symmetry 
        if (length%2 == 1):
            length += 1

        shiftX = length/2.0 - (selfX + newX)/2
        shiftY = length/2.0 - (selfY + newY)/2

        # Get the new coordinate
        selfX += shiftX
        selfY += shiftY

        newX += shiftX
        newY += shiftY

        # Generate image of multiple donuts
        imageMain = np.zeros([length, length])
        imageNeighbor = np.zeros([length, length])
        m, n = self.image.shape

        # Get the shifted main donut image
        imageMain[int(selfY-m/2):int(selfY+m/2), int(selfX-n/2):int(selfX+n/2)] += self.image

        # Get the shifted neighboring donut image
        imageNeighbor[int(newY-m/2):int(newY+m/2), int(newX-n/2):int(newX+n/2)] += magRatio*self.image

        # Get the synthesized multi-donut image
        image = imageMain + imageNeighbor

        return image, imageMain, imageNeighbor, newX, newY

    def getCenterAndR_adap(self, blockSize=33):
        """
        
        Calculate the weighting center and radius of circle based on the adapative 
        threshold. 

        Arguments:
            blockSize {[int]} -- Block size for adaptive threshold. This value should 
                                 be odd.

        Returns:
            [float] -- Values of weighting center (realcx, realcy) and radius (realR).
        """
        
        # Adaptive threshold
        delta = 1
        times = 0
        while (delta > 1e-2) and (times < 10):
            img = self.image.copy()
            imgBinary = (img > threshold_local(img, blockSize)).astype(float)

            # Calculate the weighting radius
            realR = np.sqrt(np.sum(imgBinary) / np.pi)

            # Calculte the nearest odd number of radius for the blockSize
            if (int(realR)%2 == 0):
                oddRearR = int(realR+1)
            else:
                oddRearR = int(realR)

            # Critera check of while loop
            delta = abs(blockSize - oddRearR)
            times += 1

            # New value of blockSize
            blockSize = oddRearR
        
        # Calculate the center of mass
        realcy, realcx = center_of_mass(imgBinary)

        # The values of (realcx, realcy, realR) will be (nan, nan, 0.0) for the invalid image.
        if (not np.isfinite([realcx, realcy]).any()):
            print("Can not fit donut to circle.")

        return realcx, realcy, realR, imgBinary

class AdapThresImageTest(unittest.TestCase):
    """
    
    Test the function of AdapThresImage.
    """
    
    def setUp(self):

        # Get the path of module
        modulePath = getModulePath()

        # Image file
        imageFolderPath = os.path.join(modulePath, "test", "testImages", "LSST_NE_SN25")
        imageName = 'z11_0.25_intra.txt'
        imageFile = os.path.join(imageFolderPath, imageName)

        # Set the image file
        self.adapImage = AdapThresImage()
        self.adapImage.setImg(imageFile = imageFile)

        self.zeroImage = AdapThresImage()
        self.zeroImage.setImg(image = np.zeros([120,120]))

        self.randImage = AdapThresImage()
        self.randImage.setImg(image = np.random.rand(120,120))

    def testFunc(self):

        # Calculate the centroid
        realcx, realcy, realR, imgBinary = self.adapImage.getCenterAndR_adap()
        self.assertEqual([round(realcx), round(realcy), round(realR)], [61, 61, 38])
        
        realcx, realcy, realR, imgBinary = self.adapImage.getCenterAndR_ef(checkEntropy=True)
        self.assertEqual([round(realcx), round(realcy), round(realR)], [61, 61, 38])

        realcx, realcy, realR, imgBinary = self.zeroImage.getCenterAndR_ef(checkEntropy=True)
        self.assertEqual(realcx, [])

        realcx, realcy, realR, imgBinary = self.randImage.getCenterAndR_ef(checkEntropy=True)
        self.assertEqual(realcx, [])

if __name__ == "__main__":

    # Do the unit test
    unittest.main()
