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

import os
import shutil
import unittest

from lsst.ts.wep.SourceSelector import SourceSelector
from lsst.ts.wep.Utility import getModulePath, FilterType, CamType, BscDbType
from lsst.ts.wep.bsc.PlotStarFunc import plotRaDecl, plotStarInPixelOnDetector


class TestPlotStarFunc(unittest.TestCase):
    """Test the functions of PlotStarFunc."""

    def setUp(self):

        self.modulePath = getModulePath()

        # Query the database to get the star maps
        self.sourSelc = SourceSelector(CamType.ComCam, BscDbType.LocalDb)
        self.sourSelc.setObsMetaData(0.0, 63, 0.0)
        self.sourSelc.configNbrCriteria(63.0, 2.5, maxNeighboringStar=99)
        self.sourSelc.setFilter(FilterType.U)

        dbAdress = os.path.join(self.modulePath, "tests", "testData", "bsc.db3")
        self.sourSelc.connect(dbAdress)
        neighborStarMap, starMap, wavefrontSensors = self.sourSelc.getTargetStar(
            offset=0
        )
        self.sourSelc.disconnect()

        # Collect the data for the test of plot functions
        self.wavefrontSensors = wavefrontSensors
        self.starMap = starMap
        self.neighborStarMap = neighborStarMap

        # Save image files directory
        self.dataDir = os.path.join(self.modulePath, "tests", "tmp")
        self._makeDir(self.dataDir)

    def _makeDir(self, directory):

        if not os.path.exists(directory):
            os.makedirs(directory)

    def tearDown(self):

        shutil.rmtree(self.dataDir)

    def testPlotRaDecl(self):

        saveFilePath = os.path.join(self.dataDir, "starInRaDecl.png")
        plotRaDecl(
            self.wavefrontSensors,
            self.starMap,
            self.neighborStarMap,
            self.sourSelc.db.STD_DEV_SPLIT,
            saveFilePath=saveFilePath,
        )

        self.assertTrue(os.path.exists(saveFilePath))

    def testPlotStarInPixelOnDetector(self):

        detector = list(self.neighborStarMap.keys())[0]
        stars = self.starMap[detector]
        neighboringStar = self.neighborStarMap[detector]
        saveFilePath = os.path.join(self.dataDir, "starInPixel.png")
        plotStarInPixelOnDetector(
            stars, neighboringStar, xyDim=(4000, 4000), saveFilePath=saveFilePath
        )

        self.assertTrue(os.path.exists(saveFilePath))


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
