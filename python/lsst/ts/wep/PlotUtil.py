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


if __name__ == "__main__":
    pass
