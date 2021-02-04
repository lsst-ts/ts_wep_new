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
import numpy as np

from lsst.ts.wep.cwfs.DonutTemplateModel import DonutTemplateModel
from lsst.ts.wep.Utility import DefocalType


class TestTemplateModel(unittest.TestCase):
    """Test the TemplateModel class."""

    def setUp(self):

        self.templateMaker = DonutTemplateModel()

    def testMakeTemplate(self):

        # Generate a test template on the center chip
        imageSize = 160
        templateArray = self.templateMaker.makeTemplate(
            "R22_S11", DefocalType.Extra, imageSize
        )

        self.assertTrue(isinstance(templateArray, np.ndarray))
        self.assertEqual(templateArray.dtype, int)
        self.assertEqual(np.max(templateArray), 1)

        # Center of donut should have hole in it
        self.assertEqual(templateArray[int(imageSize / 2), int(imageSize / 2)], 0)

        # Donut at center of focal plane should be symmetrical
        np.testing.assert_array_equal(
            templateArray[: int(imageSize / 2)],
            templateArray[-1 : -1 * (int(imageSize / 2) + 1) : -1],
        )

    def testTemplateSize(self):

        # Generate a test template on the center chip
        imageSize = 160
        templateArray = self.templateMaker.makeTemplate(
            "R22_S11", DefocalType.Extra, imageSize
        )

        smallTemplate = self.templateMaker.makeTemplate(
            "R22_S11", DefocalType.Extra, imageSize - 20
        )

        self.assertEqual(np.shape(templateArray), (imageSize, imageSize))
        self.assertEqual(np.shape(smallTemplate), (imageSize - 20, imageSize - 20))
        np.testing.assert_array_equal(templateArray[10:-10, 10:-10], smallTemplate)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
