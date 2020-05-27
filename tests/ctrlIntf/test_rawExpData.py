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

from lsst.ts.wep.ctrlIntf.RawExpData import RawExpData


class TestRawExpData(unittest.TestCase):
    """Test the RawExpData class."""

    def setUp(self):

        self.rawExpData = RawExpData()

    def testGetVisit(self):

        visit = self.rawExpData.getVisit()
        self.assertTrue(isinstance(visit, list))
        self.assertEqual(len(visit), 0)

    def testGetSnap(self):

        snap = self.rawExpData.getSnap()
        self.assertTrue(isinstance(snap, list))
        self.assertEqual(len(snap), 0)

    def testGetRawExp(self):

        rawExpDir = self.rawExpData.getRawExpDir()
        self.assertTrue(isinstance(rawExpDir, list))
        self.assertEqual(len(rawExpDir), 0)

    def testAppend(self):

        visit = 1
        snap = 0
        rawExpDir = "rawExpDir"
        self.rawExpData.append(visit, snap, rawExpDir)

        visitInObj = self.rawExpData.getVisit()
        snapInObj = self.rawExpData.getSnap()
        rawExpDirInObj = self.rawExpData.getRawExpDir()

        self.assertEqual(visitInObj[0], visit)
        self.assertEqual(len(visitInObj), 1)

        self.assertEqual(snapInObj[0], snap)
        self.assertEqual(len(snapInObj), 1)

        self.assertEqual(rawExpDirInObj[0], rawExpDir)
        self.assertEqual(len(rawExpDirInObj), 1)

    def testAppendWithNegativeVisit(self):

        self.assertRaises(ValueError, self.rawExpData.append, -1, 0, "temp")

    def testAppendWithNegativeSnap(self):

        self.assertRaises(ValueError, self.rawExpData.append, 0, -1, "temp")

    def testReset(self):

        self.rawExpData.append(1, 0, "temp")

        self.rawExpData.reset()
        visitInObj = self.rawExpData.getVisit()
        snapInObj = self.rawExpData.getSnap()
        rawExpDirInObj = self.rawExpData.getRawExpDir()

        self.assertEqual(len(visitInObj), 0)
        self.assertEqual(len(snapInObj), 0)
        self.assertEqual(len(rawExpDirInObj), 0)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
