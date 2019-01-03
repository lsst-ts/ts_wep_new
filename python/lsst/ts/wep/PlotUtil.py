import numpy as np

from lsst.ts.wep.ButlerWrapper import ButlerWrapper
from matplotlib.colors import LogNorm, SymLogNorm

import matplotlib
# Must be before importing matplotlib.pyplot or pylab!
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def poltExposureImage(exposure, name="", scale="log", cmap="gray", vmin=None, vmax=None, saveFilePath=None):
    """
    
    Plot the exposure image.
    
    Arguments:
        exposure {[exposure]} -- Data butler of exposure image.
    
    Keyword Arguments:
        name {[str]} -- Image title name. (default: {""})
        scale {[str]} -- Scale of image map (log or linear). (default: {"log"})
        cmap {[str]} -- Color map definition. (default: {"gray"})
        vmin {[float]} -- Mininum value to show. This normalizes the luminance data. (default: {None})
        vmax {[float]} -- Maximum value to show. This normalizes the luminance data. (default: {None})
        saveFilePath {[str]} -- Save image to file path. (default: {None})
    """
    # Get the image data
    img = ButlerWrapper.getImageData(exposure)

    # Change the scale if needed
    if scale not in ("linear", "log"):
        print("No %s scale to choose. Only 'linear' and 'log' scales are allowed." % scale)
        return

    # Decide the norm in imshow for the ploting
    if (scale == "linear"):
        plotNorm = None
    elif (scale == "log"):
        if (img.min()) < 0:
            plotNorm = SymLogNorm(linthresh=0.03)
        else:
            plotNorm = LogNorm()
    
    # Plot the image
    plt.figure()
    plt.imshow(img, cmap=cmap, origin="lower", norm=plotNorm, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.title(name)

    if (saveFilePath is not None):
        plt.savefig(saveFilePath, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def plotHist(exposure, name="", numOfBin=1000, log=False, saveFilePath=None):
    """
    
    Plot the histogram.
    
    Arguments:
        exposure {[exposure]} -- Data butler of exposure image.
    
    Keyword Arguments:
        name {string} -- Image title name (default: {""}).
        numOfBin {int} -- Number of bins (default: {1000}).
        log {bool} -- The histogram axis will be set to a log scale if log=True 
                      (default: {False}).
        saveFilePath {[str]} -- Save image to file path. (default: {None})
    """

    # Get the image data
    img = ButlerWrapper.getImageData(exposure)

    # Plot the histogram
    plt.figure()
    plt.hist(img.flatten(), bins=int(numOfBin), log=log)
    plt.title(name)

    if (saveFilePath is not None):
        plt.savefig(saveFilePath, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def plotDonutImg(donutMap, saveToDir=None, dpi=None):
    """
    
    Plot the donut image.
    
    Arguments:
        donutMap {[dict]} --  Donut image map.
    
    Keyword Arguments:
        saveToDir {[str]} -- Directory to save the images. (default: {None})
        dpi {[int]} -- The resolution in dots per inch. (default: {None})
    """

    intraType = "intra"
    extraType = "extra"

    for sensorName, donutList in donutMap.items():
        # Generate the image name
        imgTitle = abbrevDectectorName(sensorName) + "_DonutImg"

        # Collect all images and titles
        intraImgList = []
        extraImgList = []
        
        intraTitleList = []
        extraTitleList = []

        intraPixelXYList = []
        extraPixelXYList = []
        
        # Collect intra- and extra-focal donut images
        for donutImg in donutList:
            for ii in range(2):

                # Assign the image (0: intra, 1: extra)
                if (ii == 0):
                    img = donutImg.intraImg
                else:
                    img = donutImg.extraImg

                if (img is not None):

                    pixelXy = (donutImg.pixelX, donutImg.pixelY)

                    if (ii == 0):
                        intraImgList, intraTitleList, intraPixelXYList = _collectDonutImgList(
                                                intraImgList, intraTitleList, intraPixelXYList, 
                                                img, donutImg.starId, intraType, pixelXy)
                    else:
                        extraImgList, extraTitleList, extraPixelXYList = _collectDonutImgList(
                                                extraImgList, extraTitleList, extraPixelXYList, 
                                                img, donutImg.starId, extraType, pixelXy)

        # Decide the figure grid shape
        numOfRow = np.max([len(intraImgList), len(extraImgList)])

        if (len(intraImgList) == 0) or (len(extraImgList) == 0):
            numOfCol = 1
        else:
            numOfCol = 2

        gridShape = (numOfRow, numOfCol)

        # Plot the donut figure
        plt.figure()

        # Plot the intra-focal donut
        locOfCol = 0
        for ii in range(len(intraImgList)):
            _subPlot(plt, gridShape, (ii, locOfCol), intraImgList[ii], intraTitleList[ii], intraPixelXYList[ii])

        # Plot the extra-focal donut

        # Update the location of column if necessary
        if (numOfCol == 2):
            locOfCol = 1
        
        for ii in range(len(extraImgList)):
            _subPlot(plt, gridShape, (ii, locOfCol), extraImgList[ii], extraTitleList[ii], extraPixelXYList[ii])

        # Adjust the space between xlabel and title for neighboring sub-figures
        plt.tight_layout()

        # Save the file or not
        if (saveToDir is not None):
            # Generate the filepath
            imgType = ".png"
            imgFilePath = os.path.join(saveToDir, imgTitle+imgType)
            plt.savefig(imgFilePath, bbox_inches="tight", dpi=dpi)
            plt.close()
        else:
            plt.show()


def _subPlot(plt, gridShape, loc, img, aTitle, pixelXy):
    """
    
    Do the subplot of figure.
    
    Arguments:
        plt {[pyplot]} -- Plotting framework.
        gridShape {[tuple]} -- Shape of grid.
        loc {[tuple]} -- Location of subplot in grid.
        img {[ndarray]} -- Image of donut.
        aTitle {[str]} -- Title of subplot.
        pixelXy {[tuple]} -- Chip position of donut in (x, y).
    """

    # Chip position of donut
    pixelXy = np.round(pixelXy)
    pixelPos = "Pixel XY: (%d, %d)" % (pixelXy[0], pixelXy[1])

    # Decide the position of subplot
    ax = plt.subplot2grid(gridShape, loc)

    # Show the figure
    axPlot = ax.imshow(img, origin="lower")
    
    # Set the title
    ax.set_title(aTitle)

    # Set the x lavel
    ax.set_xlabel(pixelPos)

    # Set the colorbar
    plt.colorbar(axPlot, ax=ax)


def _collectDonutImgList(imgList, titleList, pixelXyList, img, starId, aType, pixelXy):
    """
    
    Collect the donut data in list.
    
    Arguments:
        imgList {[list]} -- List of image.
        titleList {[list]} -- List of title.
        pixelXyList {[list]} -- List of pixel XY.
        img {[ndarray]} -- Donut image.
        starId {[int]} -- Star Id.
        aType {[str]} -- Type of donut.
        pixelXy {[tuple]} -- Pixel position in (x, y).
    
    Returns:
        [list] -- List of image.
        [list] -- List of title.
        [list] -- List of pixel XY.
    """

    # Get the title
    aTitle = "_".join([str(starId), aType])    

    # Append the list
    imgList.append(img)
    titleList.append(aTitle)
    pixelXyList.append(pixelXy)

    return imgList, titleList, pixelXyList


if __name__ == "__main__":
    pass
