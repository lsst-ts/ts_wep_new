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
import tempfile

from lsst.ts.wep.Utility import getModulePath


class BaseBscTestCase(object):
    """Base class for the bright star catalog (BSC) tests"""

    def createBscTest(self):
        """Create the bright star catalog (BSC) used in the test."""

        # Create a temporary test directory
        testDir = os.path.join(getModulePath(), "tests")
        self._tempDir = tempfile.TemporaryDirectory(dir=testDir)

        # Copy the db3 database
        dbAdressSrc = os.path.join(testDir, "testData", "bsc.db3")
        self._dbAdress = os.path.join(self._tempDir.name, "bsc.db3")

        shutil.copy(dbAdressSrc, self._dbAdress)

    def removeBscTest(self):
        """Remove the bright star catalog (BSC) used in the test.

        Need to call the createBscTest() first.
        """

        self._tempDir.cleanup()

    def getPathOfBscTest(self):
        """Get the path of bright star catalog (BSC) used in the test.

        Need to call the createBscTest() first.
        """

        return self._dbAdress
