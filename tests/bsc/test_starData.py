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

import numpy as np
import unittest

from lsst.ts.wep.bsc.StarData import StarData
from lsst.ts.wep.Utility import FilterType


class TestStarData(unittest.TestCase):
    """Test the StarData class."""

    def setUp(self):
        self.stars = StarData(
            [123, 456, 789],
            [0.1, 0.2, 0.3],
            [2.1, 2.2, 2.3],
            [2.0, 3.0, 4.0],
            [2.1, 2.1, 4.1],
            [2.2, 3.2, 4.2],
            [2.3, 3.3, 4.3],
            [2.4, 3.4, 4.4],
            [2.5, 3.5, 4.5],
        )

    def testGetId(self):

        starId = self.stars.getId()
        self.assertEqual(starId.dtype, int)
        self.assertEqual(starId.tolist(), [123, 456, 789])

    def testGetRA(self):

        self.assertEqual(self.stars.getRA().tolist(), [0.1, 0.2, 0.3])

    def testGetDecl(self):

        self.assertEqual(self.stars.getDecl().tolist(), [2.1, 2.2, 2.3])

    def testGetMag(self):

        self.assertEqual(self.stars.getMag(FilterType.U).tolist(), [2.0, 3.0, 4.0])
        self.assertEqual(self.stars.getMag(FilterType.G).tolist(), [2.1, 2.1, 4.1])
        self.assertEqual(self.stars.getMag(FilterType.R).tolist(), [2.2, 3.2, 4.2])
        self.assertEqual(self.stars.getMag(FilterType.I).tolist(), [2.3, 3.3, 4.3])
        self.assertEqual(self.stars.getMag(FilterType.Z).tolist(), [2.4, 3.4, 4.4])
        self.assertEqual(self.stars.getMag(FilterType.Y).tolist(), [2.5, 3.5, 4.5])

    def testSetMag(self):

        mag = [1, 3, 4, 5]
        self.stars.setMag(FilterType.U, mag)

        self.assertEqual(self.stars.getMag(FilterType.U).tolist(), mag)

    def testSetAndGetDetector(self):

        detector = "CCD"
        self.stars.setDetector(detector)
        self.assertEqual(self.stars.getDetector(), detector)

    def testSetAndGetRaInPixel(self):

        raInPixel = [1.0, 2.0]
        self.stars.setRaInPixel(raInPixel)

        self.assertEqual(self.stars.getRaInPixel().tolist(), raInPixel)

    def testSetRaInPixelWithFloatValue(self):

        raInPixel = 1.0
        self.stars.setRaInPixel(raInPixel)

        self.assertEqual(self.stars.getRaInPixel().tolist(), [raInPixel])

    def testSetRaInPixelWithNpArray(self):

        raInPixel = np.array([1.0, 2.0])
        self.stars.setRaInPixel(raInPixel)

        delta = np.sum(np.abs(self.stars.getRaInPixel() - raInPixel))
        self.assertEqual(delta, 0)

    def testSetAndGetDeclInPixel(self):

        declInPixel = [2.0, 3.0]
        self.stars.setDeclInPixel(declInPixel)

        self.assertEqual(self.stars.getDeclInPixel().tolist(), declInPixel)

    def testCheckCandidateStars(self):

        indexCandidateU = self.stars.checkCandidateStars(FilterType.U, 1.9, 2.1)
        indexCandidateG = self.stars.checkCandidateStars(FilterType.G, 0, 5)
        indexCandidateR = self.stars.checkCandidateStars(FilterType.R, 0, 1)
        indexCandidateI = self.stars.checkCandidateStars(FilterType.I, 2.1, 4.0)
        indexCandidateZ = self.stars.checkCandidateStars(FilterType.Z, 3.0, 5.0)
        indexCandidateY = self.stars.checkCandidateStars(FilterType.Y, 1.0, 2.0)

        self.assertEqual(indexCandidateU, [0])
        self.assertEqual(indexCandidateG, [0, 1, 2])
        self.assertEqual(indexCandidateR, [])
        self.assertEqual(indexCandidateI, [0, 1])
        self.assertEqual(indexCandidateZ, [1, 2])
        self.assertEqual(indexCandidateY, [])

    def testGetNeighboringStar(self):

        self._populateRaDeclInPixel()

        neighboringStarU = self.stars.getNeighboringStar(
            [0], 3, FilterType.U, maxNumOfNbrStar=99
        )
        neighboringStarG = self.stars.getNeighboringStar(
            [0, 1], 3, FilterType.G, maxNumOfNbrStar=99
        )
        neighboringStarR = self.stars.getNeighboringStar(
            [0], 1, FilterType.R, maxNumOfNbrStar=99
        )
        neighboringStarI = self.stars.getNeighboringStar(
            [], 3, FilterType.I, maxNumOfNbrStar=99
        )
        neighboringStarZ = self.stars.getNeighboringStar(
            [0, 1], 2, FilterType.Z, maxNumOfNbrStar=1
        )
        neighboringStarY = self.stars.getNeighboringStar(
            [1], 2, FilterType.Y, maxNumOfNbrStar=1
        )

        self.assertEqual(len(neighboringStarU.getId()[123]), 2)
        self.assertEqual(len(neighboringStarU.getRaDecl()), 3)

        self.assertEqual(len(neighboringStarG.getId()), 2)

        self.assertEqual(len(neighboringStarR.getId()[123]), 0)

        self.assertEqual(neighboringStarI.getId(), {})

        self.assertEqual(len(neighboringStarZ.getId()), 1)

        self.assertEqual(neighboringStarY.getId(), {})

    def _populateRaDeclInPixel(self):

        self.stars.setRaInPixel(self.stars.getRA() * 10)
        self.stars.setDeclInPixel(self.stars.getDecl() * 10)

    def testGetNeighboringStarWithNothing(self):

        self._populateRaDeclInPixel()
        neighboringStarU = self.stars.getNeighboringStar(
            [], 3.0, FilterType.U, maxNumOfNbrStar=99
        )

        self.assertEqual(len(neighboringStarU.getId()), 0)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
