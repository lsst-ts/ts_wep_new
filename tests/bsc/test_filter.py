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

from lsst.ts.wep.Utility import FilterType
from lsst.ts.wep.bsc.Filter import Filter


class TestFilter(unittest.TestCase):
    """Test the Filter class."""

    def setUp(self):

        self.filter = Filter()

    def testGetFilter(self):

        self.assertEqual(self.filter.getFilter(), FilterType.U)

    def testSetFilter(self):

        filterType = FilterType.G
        self.filter.setFilter(filterType)
        self.assertEqual(self.filter.getFilter(), filterType)

    def testGetMagBoundaryOfFilterG(self):

        self.filter.setFilter(FilterType.G)

        lowMagnitude, highMagnitude = self.filter.getMagBoundary()
        self.assertEqual(lowMagnitude, Filter.G_LOW_MAG)
        self.assertEqual(highMagnitude, Filter.G_HIGH_MAG)

    def testGetMagBoundaryOfFilterU(self):

        self.filter.setFilter(FilterType.U)

        lowMagnitude, highMagnitude = self.filter.getMagBoundary()
        self.assertEqual(lowMagnitude, Filter.U_LOW_MAG)
        self.assertEqual(highMagnitude, Filter.U_HIGH_MAG)

    def testGetMagBoundaryOfFilterRef(self):

        self.filter.setFilter(FilterType.REF)

        lowMagnitude, highMagnitude = self.filter.getMagBoundary()
        self.assertEqual(lowMagnitude, Filter.G_LOW_MAG)
        self.assertEqual(highMagnitude, Filter.G_HIGH_MAG)

    def testGetMagBoundaryWithError(self):

        self.filter.setFilter("r")
        self.assertRaises(ValueError, self.filter.getMagBoundary)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
