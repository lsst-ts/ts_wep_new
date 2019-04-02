import os
import unittest

from lsst.ts.wep.cwfs.Instrument import Instrument
from lsst.ts.wep.Utility import getModulePath


class TestInstrument(unittest.TestCase):
    """Test the Instrument class."""

    def setUp(self):

        # Get the path of module
        modulePath = getModulePath()

        # Define the instrument folder
        self.instruFolder = os.path.join(modulePath, "configData", "cwfs",
                                         "instruData")

        # Define the instrument name
        self.instruName = "lsst"

    def testInstrument(self):
        inst = Instrument(self.instruFolder)
        inst.config(self.instruName, 120)
        self.assertEqual(inst.parameter["sensorSamples"], 120)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
