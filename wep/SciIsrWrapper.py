import os, unittest
import numpy as np

import matplotlib
# Must be before importing matplotlib.pyplot or pylab!
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm

from lsst.daf.persistence import Butler

from wep.SciWFDataCollector import runProgram

class SciIsrWrapper(object):

    def __init__(self):
        """
        
        Initialize the SciIsrWrapper class.        
        """

        self.pathData = None
        self.outputPath = None
        self.butler = None

    def configWrapper(self, inputs=None, outputs=None):
        """
        
        Configurate the ISR wrapper.

        Keyword Arguments:
            inputs {[str]} -- Input directory. (default: {None})
            outputs {[str]} -- Output directory. (default: {None})
        """
        
        self.pathData = inputs
        self.outputPath = outputs

    def configBulter(self, inputs, outputs=None):
        """
        
        Configurate the data butler.

        Arguments:
              inputs {[RepositoryArg or string]} -- Can be a single item or a list. Provides arguments 
                      to load an existing repository (or repositories). String is assumed to be a URI 
                      and is used as the cfgRoot (URI to the location of the cfg file). (Local file system 
                      URI does not have to start with 'file://' and in this way can be a relative path).
        
        Keyword Arguments:                      
              outputs {[RepositoryArg or string]} -- Can be a single item or a list. Provides arguments to 
                      load one or more existing repositories or create new ones. String is assumed to be a 
                      URI and as used as the repository root. (default: {None})
        """

        self.butler = Butler(inputs=inputs, outputs=outputs)

    def doISR(self, visit, snap, raft, sensor):
        """
        
        Do the instrument signature removal (ISR). 
        
        Arguments:
            visit {[int]} -- Visit time.
            snap {int} -- Snap time (0 or 1) means first/ second exposure.
            raft {[str]} -- Raft name.
            sensor {[str]} -- Sensor name.
        """
        
        if (self.pathData is not None) and (self.outputPath is not None):

            # Process the phoSim CCD images
            command = "processSimCcd.py"

            # Set data Id
            dataId = "visit=%d raft=%s sensor=%s snap=%d" % (visit, raft, sensor, snap)

            # Set the input and output directories
            argstring = "%s --id %s --output %s" % (self.pathData, dataId, self.outputPath)

            # Do the ISR and assemble CCD image
            runProgram(command, argstring=argstring)

        else:
            print("Can not do ISR because input/output path is not defined.")

    def getButlerData(self, visit, snap, raft, sensor, channel=None, atype=None):
        """
        
        Get the butler data.
        
        Arguments:
            visit {[int]} -- Visit time.
            snap {int} -- Snap time (0 or 1) means first/ second exposure.
            raft {[string]} -- Raft name.
            sensor {[string]} -- Sensor name.

        Keyword Arguments:
            channel {[string]} -- Channel name (default: {"None"}).
            atype {[string]} -- Type of arrangement: eimage, raw, bias, flat, dark, postISRCCD 
                                (default: {"None"}).
        
        Returns:
            [butler] -- Butler data.
        """

        # Define the data ID and get the raw amplifer image
        if atype in ("eimage", "postISRCCD"):
            dataId = dict(visit=int(visit), snap=snap, raft=raft, sensor=sensor, immediate=True)
        else:
            dataId = dict(visit=int(visit), snap=snap, raft=raft, sensor=sensor, channel=channel, 
                            immediate=True)
        butlerData = self.butler.get(atype, dataId=dataId)

        return butlerData

def getImageData(image):
    """
    
    Get the image data in numpy array.
    
    Arguments:
        image {[obj]} -- Image data. It can be numpy array or Exposure image.
    
    Returns:
        [float] -- Image data in numpy array.
    """

    # Get the numpy array data based on the input object type
    if isinstance(image, np.ndarray):
        data = image
    elif hasattr(image, "getMaskedImage"):
        data = image.getMaskedImage().getImage().getArray() 
    elif hasattr(image, "getImage"):
        data = image.getImage().getArray() 
    else:
        data = image.getArray()

    # Return the data in numpy array
    return data

def poltExposureImage(exposure, name="", scale="log", cmap="gray", vmin=None, vmax=None):
    """
    
    Plot the exposure image.
    
    Arguments:
        exposure {[exposure]} -- Data butler of exposure image.
    
    Keyword Arguments:
        name {string} -- Image title name (default: {""}).
        scale {string} -- Scale of image map (log or linear) (default: {"log"}).
        cmap {string} -- Color map definition (default: {"gray"}).
        vmin {[float]} -- Mininum value to show. This normalizes the luminance data (default: {None}).
        vmax {[float]} -- Maximum value to show. This normalizes the luminance data (default: {None}).      
    """
    # Get the image data
    img = getImageData(exposure)

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
    plt.show()

def plotHist(exposure, name="", numOfBin=1000, log=False):
    """
    
    Plot the histogram.
    
    Arguments:
        exposure {[exposure]} -- Data butler of exposure image.
    
    Keyword Arguments:
        name {string} -- Image title name (default: {""}).
        numOfBin {int} -- Number of bins (default: {1000}).
        log {bool} -- The histogram axis will be set to a log scale if log=True 
                      (default: {False}).
    """

    # Get the image data
    img = getImageData(exposure)

    # Plot the histogram
    plt.figure()
    plt.hist(img.flatten(), bins=int(numOfBin), log=log)
    plt.title(name)
    plt.show()

class SciIsrWrapperTest(unittest.TestCase):


    """ 
    Test the function of SciIsrWrapper.
    """

    def setUp(self):

        # Path of data folder
        dataFolderPath = "../test"
        self.pathData = dataFolderPath
        self.outputPath = dataFolderPath

    def testFunction(self):
        
        # Instantiate the isr wrapper
        isrWrapper = SciIsrWrapper()

        # Configurate the path
        isrWrapper.configWrapper(inputs=self.pathData, outputs=self.outputPath)

        self.assertEqual(isrWrapper.pathData, self.pathData)
        self.assertEqual(isrWrapper.outputPath, self.outputPath)

        # Test to get the butler data
        obsId = 99999999
        snap = 0
        raft = "2,2"
        sensor = "1,1"
        channel = "1,4"
        isrWrapper.configBulter(self.pathData)
        butlerDataRaw = isrWrapper.getButlerData(obsId, snap, raft, sensor, channel=channel, atype="raw")
        self.assertEqual(butlerDataRaw.getMetadata().get("GAIN"), 1.83546)

        # Test to get the image data
        data = getImageData(butlerDataRaw)
        self.assertEqual(data.shape, (2001, 513))

if __name__ == "__main__":

    # Do the unit test
    unittest.main()
    
    # # Path
    # pathData = "/home/ttsai/Document/phosimObsData/input"
    # outputPath = "/home/ttsai/Document/phosimObsData/output"

    # # Instantiate the isr wrapper
    # isrWrapper = SciIsrWrapper()

    # # Configurate the path
    # isrWrapper.configWrapper(inputs=pathData, outputs=outputPath)

    # # Do the isr
    # visit = 9006000
    # snap = 0
    # raft = "2,2"
    # sensor = "1,1"
    # # isrWrapper.doISR(visit, snap, raft, sensor)
    # isrWrapper.configBulter(outputPath)
    # butlerData = isrWrapper.getButlerData(visit, snap, raft, sensor, atype="postISRCCD")
