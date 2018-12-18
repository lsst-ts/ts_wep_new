import os
import unittest

from lsst.daf.persistence import Butler

from lsst.ts.wep.isr.SliceFacade import SliceFacade
from lsst.ts.wep.Utility import getModulePath


class SliceFacadeTest(unittest.TestCase):
    """Test the SliceFacade class."""

    def setUp(self):

        # Get the path of module
        modulePath = getModulePath()

        # Path of data folder
        dataFolderPath = os.path.join(modulePath, "tests", "testData")
        self.dataFolderPath = dataFolderPath

    def testFunction(self):

        # Constuct the butler
        butler = Butler(inputs=self.dataFolderPath)

        # Get the amplifier slice data
        obsId = 99999999
        snap = 0
        raft = "2,2"
        sensor = "1,1"
        channel = "1,4"

        dataId = dict(visit=obsId, snap=snap, raft=raft, sensor=sensor,
                      channel=channel)
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
