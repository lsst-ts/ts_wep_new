import os, unittest
import lsst.daf.persistence as dafPersistence
from lsst.ts.wep.Utility import getModulePath

class SliceFacade(object):

    def __init__(self, dataButlerObject):

        self.dataButlerObject = dataButlerObject

    def getDimensions(self, key=None):
        """

        Get the dimension of exposure image.

        Arguments:
            key {[string]} -- Name of amplifier. (default: {None})

        Returns:
            [float] -- Dimension of exposure image.
        """

        return self.__getItem(key).getDimensions()

    def getBoresightRotAngle(self, key=None):
        """

        Get rotation angle at boresight at middle of exposure. The meaning of rotation
        angle depends on rotType. For example, if rotType is SKY the angle is the position
        angle of the focal plane +Y with respect to North.

        Arguments:
            key {[string]} -- Name of amplifier. (default: {None})

        Returns:
            [float] -- Rotation angle at boresight.
        """

        return self.__getItem(key).getInfo().getVisitInfo().getBoresightRotAngle()

    def getRotType(self, key=None):
        """

        Get rotation type of boresightRotAngle

        Arguments:
            key {[string]} -- Name of amplifier. (default: {None})

        Returns:
            [string] -- Get the rotation type.
        """

        return self.__getItem(key).getInfo().getVisitInfo().getRotType()

    def getBoresightAzAlt(self, key=None):
        """

        Get the boresight of telescope in Az and Alt.

        Arguments:
            key {[string]} -- Name of amplifier. (default: {None})

        Returns:
            [float] -- Boresight of telescope in Az and Alt.
        """

        return self.__getItem(key).getInfo().getVisitInfo().getBoresightAzAlt()

    def getBoresightRaDec(self, key=None):
        """

        Get the boresight of telescope in Ra and Dec.

        Arguments:
            key {[string]} -- Name of amplifier. (default: {None})

        Returns:
            [float] -- Boresight of telescope in Ra and Dec.
        """

        return self.__getItem(key).getInfo().getVisitInfo().getBoresightRaDec()

    def getGain(self, key=None):
        """

        Get the gain value of amplifier.

        Arguments:
            key {[string]} -- Name of amplifier. (default: {None})

        Returns:
            [float] -- Gain of amplifier.
        """

        return self.__getItem(key).getMetadata().get("GAIN")

    def __getItem(self, key=None):
        """
        
        Get the specific exposure.
        
        Keyword Arguments:
            key {[string]} -- Name of amplifier. (default: {None})
        
        Returns:
            [butler] -- Exposure butler object.
        """

        if key is not None:
            return self.dataButlerObject[key]
        else:
            return self.dataButlerObject

class SliceFacadeTest(unittest.TestCase):

    """
    Test the function of SliceFacade.
    """

    def setUp(self):

        # Get the path of module
        modulePath = getModulePath()

        # Path of data folder
        dataFolderPath = os.path.join(modulePath, "test")
        self.dataFolderPath = dataFolderPath

    def testFunction(self):

        # Constuct the butler
        butler = dafPersistence.Butler(inputs=self.dataFolderPath)

        # Get the amplifier slice data
        obsId = 99999999
        snap = 0
        raft = "2,2"
        sensor = "1,1"
        channel = "1,4"

        dataId = dict(visit=obsId, snap=snap, raft=raft, sensor=sensor, channel=channel)
        exposure = butler.get("raw", dataId=dataId)
        ampSlice = SliceFacade(exposure)
        ampTwoSlice = SliceFacade(dict([("S1", exposure), ("S2", exposure)]))

        # Test Slice functions
        self.assertEqual(ampSlice.getDimensions()[0], 513)
        self.assertEqual(ampSlice.getDimensions()[1], 2001)

        self.assertEqual(ampSlice.getGain(), 1.83546)

        self.assertEqual(ampSlice.getBoresightAzAlt()[0], 0)

        self.assertEqual(ampSlice.getBoresightRotAngle(), 0)

        self.assertEqual(ampSlice.getBoresightRaDec()[0], 0)
        self.assertEqual(ampSlice.getBoresightRaDec()[1], 0)

        # Test Slice functions for a dictionary
        self.assertEqual(ampTwoSlice.getDimensions("S1")[0], 513)
        self.assertEqual(ampTwoSlice.getDimensions("S2")[1], 2001)

if __name__ == "__main__":

    # Do the unit test
    unittest.main()