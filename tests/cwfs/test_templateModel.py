# This file is part of ts_wep.
#
# Developed for the LSST Telescope and Site Systems.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import unittest
import os
import numpy as np

from lsst.ts.wep.cwfs.TemplateModel import TemplateModel
from lsst.ts.wep.Utility import CamType, getConfigDir
from lsst.ts.wep.cwfs.Instrument import Instrument
from lsst.ts.wep.cwfs.CompensableImage import CompensableImage


class TestTemplateModel(unittest.TestCase):
    """Test the TemplateModel class."""

    def setUp(self):

        self.templateMaker = TemplateModel()

    def _createTestTemplate(self, imageSize, defocalState):
        """Create a test template located at center of focal plane."""

        # Load Instrument parameters
        instDir = os.path.join(getConfigDir(), "cwfs", "instData")
        dimOfDonutImgOnSensor = imageSize
        inst = Instrument(instDir)
        inst.config(CamType.LsstCam, dimOfDonutImgOnSensor)

        # Create image for mask
        img = CompensableImage()
        img.defocalType = defocalState

        # define position of donut at center of current sensor in degrees
        img.fieldX, img.fieldY = 0.0, 0.0
        img.makeMask(inst, "offAxis", 0, 1)

        return img.cMask

    def testMakeTemplate(self):

        testTemplate = self._createTestTemplate(160, "extra")

        # Generate a test template on the center chip
        templateArray = self.templateMaker.makeTemplate("R22_S11", "extra", 160)

        self.assertEqual(np.max(templateArray), 1.0)
        self.assertEqual(np.shape(templateArray), (160, 160))
        np.testing.assert_array_equal(templateArray, testTemplate)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
