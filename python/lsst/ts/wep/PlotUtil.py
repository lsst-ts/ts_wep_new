import os
import numpy as np

from lsst.ts.wep.ButlerWrapper import ButlerWrapper
from lsst.ts.wep.Utility import abbrevDectectorName
from matplotlib.colors import LogNorm, SymLogNorm

import matplotlib.pyplot as plt
plt.switch_backend('Agg')


def poltExposureImage(exposure, name="", scale="log", cmap="gray", vmin=None,
                      vmax=None, saveFilePath=None):
    """Plot the exposure image.

    Parameters
    ----------
    exposure : Exposure
        Data butler of exposure image.
    name : str, optional
        Image title name. (the default is "".)
    scale : str, optional
        Scale of image map (log or linear). (the default is "log".)
    cmap : str, optional
        Color map definition. (the default is "gray".)
    vmin : float, optional
        Mininum value to show. This normalizes the luminance data. (the default
        is None.)
    vmax : float, optional
        Maximum value to show. This normalizes the luminance data. (the default
        is None.)
    saveFilePath : str, optional
        Save image to file path. (the default is None.)
    """

    # Get the image data
    img = ButlerWrapper.getImageData(exposure)

    # Change the scale if needed
    if scale not in ("linear", "log"):
        print("No %s scale to choose. Only 'linear' and 'log' scales are allowed." % scale)
        return

    # Decide the norm in imshow for the ploting
    if scale == "linear":
        plotNorm = None
    elif scale == "log":
        if (img.min()) < 0:
            plotNorm = SymLogNorm(linthresh=0.03)
        else:
            plotNorm = LogNorm()

    # Plot the image
    plt.figure()
    plt.imshow(img, cmap=cmap, origin="lower", norm=plotNorm, vmin=vmin,
               vmax=vmax)
    plt.colorbar()
    plt.title(name)

    if saveFilePath is None:
        plt.show()
    else:
        plt.savefig(saveFilePath, bbox_inches="tight")
        plt.close()


def plotHist(exposure, name="", numOfBin=1000, log=False, saveFilePath=None):
    """Plot the histogram.

    Parameters
    ----------
    exposure : Exposure
        Data butler of exposure image.
    name : str, optional
        Image title name. (the default is "".)
    numOfBin : int, optional
        Number of bins. (the default is 1000.)
    log : bool, optional
        The histogram axis will be set to a log scale if log=True. (the default
        is False.)
    saveFilePath : str, optional
        Save image to file path. (the default is None.)
    """

    # Get the image data
    img = ButlerWrapper.getImageData(exposure)

    # Plot the histogram
    plt.figure()
    plt.hist(img.flatten(), bins=int(numOfBin), log=log)
    plt.title(name)

    if saveFilePath is None:
        plt.show()
    else:
        plt.savefig(saveFilePath, bbox_inches="tight")
        plt.close()


def plotDonutImg(donutMap, saveToDir=None, dpi=None):
    """Plot the donut image.

    Parameters
    ----------
    donutMap : dict
        Donut image map.
    saveToDir : str, optional
        Directory to save the images. (the default is None.)
    dpi : int, optional
        The resolution in dots per inch. (the default is None.)
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
                if ii == 0:
                    img = donutImg.intraImg
                else:
                    img = donutImg.extraImg

                if img is not None:

                    pixelXy = (donutImg.pixelX, donutImg.pixelY)

                    if ii == 0:
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
            _subPlot(plt, gridShape, (ii, locOfCol), intraImgList[ii],
                     intraTitleList[ii], intraPixelXYList[ii])

        # Plot the extra-focal donut

        # Update the location of column if necessary
        if numOfCol == 2:
            locOfCol = 1

        for ii in range(len(extraImgList)):
            _subPlot(plt, gridShape, (ii, locOfCol), extraImgList[ii],
                     extraTitleList[ii], extraPixelXYList[ii])

        # Adjust the space between xlabel and title for neighboring sub-figures
        plt.tight_layout()

        # Save the file or not
        if saveToDir is not None:
            # Generate the filepath
            imgType = ".png"
            imgFilePath = os.path.join(saveToDir, imgTitle+imgType)
            plt.savefig(imgFilePath, bbox_inches="tight", dpi=dpi)
            plt.close()
        else:
            plt.show()


def _subPlot(plt, gridShape, loc, img, aTitle, pixelXy):
    """Do the subplot of figure.

    Parameters
    ----------
    plt : pyplot
        Plotting framework.
    gridShape : tuple
        Shape of grid.
    loc : tuple
        Location of subplot in grid.
    img : numpy.ndarray
        Image of donut.
    aTitle : str
        Title of subplot.
    pixelXy : tuple
        Chip position of donut in (x, y).
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


def _collectDonutImgList(imgList, titleList, pixelXyList, img, starId, aType,
                         pixelXy):
    """Collect the donut data in list.

    Parameters
    ----------
    imgList : list
        List of image.
    titleList : list
        List of title.
    pixelXyList : list
        List of pixel XY.
    img : numpy.ndarray
        Donut image.
    starId : int
        Star Id.
    aType : str
        Type of donut.
    pixelXy : tuple
        Pixel position in (x, y).

    Returns
    -------
    list
        List of image.
    list
        List of title.
    list
        List of pixel XY.
    """

    # Get the title
    aTitle = "_".join([str(starId), aType])

    # Append the list
    imgList.append(img)
    titleList.append(aTitle)
    pixelXyList.append(pixelXy)

    return imgList, titleList, pixelXyList


def plotImage(imageDonut, title=None, mask=None, show=True, fitParameters=[],
              saveFilePath=None):
    """Show the wavefront (donut) image.

    Parameters
    ----------
    imageDonut : numpy.ndarray[float]
        Input donut Image.
    title : str, optional
        Title of image. (the default is None.)
    mask : numpy.ndarray[int], optional
        Mask. (the default is None.)
    show : bool, optional
        Show the figure or not. (the default is True.)
    fitParameters : list, optional
        Fitting parameter of circle (center position and radius) (the default
        is [].)
    saveFilePath : str, optional
        File path to save the image. (the default is None.)
    """

    if mask is not None:
        image = np.where(mask == 0, np.nan, imageDonut)
    else:
        image = imageDonut

    plt.figure()

    # Plot the fitted circle
    if fitParameters:
        if len(fitParameters) == 3:
            theta = np.linspace(0, 2*np.pi, 101)
            x = fitParameters[0] + fitParameters[2]*np.cos(theta)
            y = fitParameters[1] + fitParameters[2]*np.sin(theta)
            plt.plot(x, y, "b")
        else:
            print("fitParameters should have 3 elements.")

    # Plot the wavefront image
    plt.imshow(image, origin="lower")
    plt.colorbar()
    if title:
        plt.title(title)

    if saveFilePath is not None:
        plt.savefig(saveFilePath, bbox_inches="tight")
        plt.close()

    if show:
        plt.show()


def plotZernike(zkIdx, zk, unit, saveFilePath=None):
    """Plot the Zernike polynomials (zk).

    Parameters
    ----------
    zkIdx : list[int] or numpy.array[int]
        Index of zk.
    zk : numpy.array
        Zernike polynomials.
    unit : str
        Unit of Zernike polynomials.
    saveFilePath : str, optional
        File path to save the image. (the default is None.)
    """

    plt.plot(zkIdx, zk, marker="o", color="r", markersize=10)

    plt.xlabel("Zernike Index")
    plt.ylabel("Zernike coefficient (%s)" % unit)
    plt.grid()

    if saveFilePath is None:
        plt.show()
    else:
        plt.savefig(saveFilePath, bbox_inches="tight")
        plt.close()
