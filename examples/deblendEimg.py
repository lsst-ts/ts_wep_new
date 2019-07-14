import os
import numpy as np

from lsst.ts.wep.deblend.AdapThresImage import AdapThresImage
from lsst.ts.wep.deblend.BlendedImageDecorator import BlendedImageDecorator
from lsst.ts.wep.PlotUtil import plotImage
from lsst.ts.wep.Utility import getModulePath


if __name__ == "__main__":

    # Get the path of module
    modulePath = getModulePath()

    # Define the image folder and image names
    # Image data -- Don't know the final image format.
    # It is noted that image.readFile inuts is based on the txt file
    imageFolderPath = os.path.join(modulePath, "tests", "testData",
                                   "testImages", "LSST_NE_SN25")
    imageName = 'z11_0.25_intra.txt'
    imageFile = os.path.join(imageFolderPath, imageName)

    imageIntra = AdapThresImage()
    imageIntra.setImg(imageFile=imageFile)

    # Coefficient of distance between donuts
    # (Space coefficient should be >= 1.2)
    spaceCoef = np.random.rand()*1.5 + 1.2
    print("Random space coeffcient is %f." % spaceCoef)

    # Ratio of magnitude between donuts (If the magnitudes of stars differ by
    # 5, the brightness differs by 100.)
    # (Magnitude difference shoulbe be >= 1.)
    magDiff = np.random.rand()*3 + 1.0
    magRatio = 1/100**(magDiff/5.0)
    print("Random magnitude difference is %f and the ratio is %f." % (
        magDiff, magRatio))

    # Theta angle of generated neighboring star
    theta = np.random.rand()*2*np.pi
    print("Random theta in radian is %f." % theta)

    # Generate the blended image
    image, imageMain, imageNeighbor, neighborX, neighborY = \
        imageIntra.generateMultiDonut(spaceCoef, magRatio, theta)

    # Add a constant system error
    sysError = np.random.rand()
    image += sysError
    print("Random system error is %f." % sysError)

    # Final blended image
    # blendImage = BlendedImage.BlendedImage(image=image, atype="intra")
    blendImage = BlendedImageDecorator()
    blendImage.setImg(image=image)

    # Generate a random error to simulate the uncertainty of position of
    # neighboring star
    tempX = np.random.randint(-5, high=6)
    tempY = np.random.randint(-5, high=6)
    print("Random shift (X, Y) of neighboring star is (%f, %f)." % (
        tempX, tempY))

    # Add the random error to position of neighboring star
    neighborX += tempX
    neighborY += tempY

    # Generate a random error to simulate the uncertainty of magnitude ratio of
    # neighboring star
    if (np.random.rand() > 0.5):
        tmpMag = np.random.rand()*0.3
    else:
        tmpMag = -np.random.rand()*0.3

    print("Random shift of magnitude ratio of neighboring star is %f." % tmpMag)

    # Do the deblending
    imgDeblend, realcx, realcy = blendImage.deblendDonut((neighborX, neighborY))

    # Do the comaprison
    delta = np.sum(np.abs(imageMain-imgDeblend))
    diffRatio = delta/np.sum(np.abs(imageMain))
    print("Difference is %f." % delta)
    print("Ratio of difference is %f." % diffRatio)
    print("Pixel x, y of donut is (%f, %f)." % (realcx, realcy))

    # Show the images
    plotImage(image, title="Blended Image", saveFilePath="temp1.png")
    plotImage(imgDeblend, title="Deblended Image", saveFilePath="temp2.png")
