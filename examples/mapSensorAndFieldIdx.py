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

from lsst.ts.wep.SourceProcessor import SourceProcessor


if __name__ == "__main__":

    # Configurate the source processor
    sourPro = SourceProcessor()

    # Get the list of sensor name and coordinate
    sensorFocaPlaneInDeg = sourPro.sensorFocaPlaneInDeg

    nameList = []
    xList = []
    yList = []
    for aKey, aItem in sensorFocaPlaneInDeg.items():
        nameList.append(aKey)
        xList.append(aItem[0])
        yList.append(aItem[1])

    xySensor = np.array([xList, yList]).T

    # Set the 35 field points
    nArm = 6
    armLen = [0.379, 0.841, 1.237, 1.535, 1.708]
    fieldWFSx = [1.176, -1.176, -1.176, 1.176]
    fieldWFSy = [1.176, 1.176, -1.176, -1.176]
    pointAngle = np.arange(nArm) * (2 * np.pi) / nArm
    fieldX = np.concatenate(
        [np.zeros(1), np.kron(armLen, np.cos(pointAngle)), fieldWFSx]
    )
    fieldY = np.concatenate(
        [np.zeros(1), np.kron(armLen, np.sin(pointAngle)), fieldWFSy]
    )

    xyField = np.array([fieldX, fieldY]).T

    # Calculate the distance matrix (sensor by field)
    disM = np.zeros((len(xList), len(fieldX)))
    for ii in range(len(fieldX)):
        vector = xySensor - xyField[ii, :]
        dis = np.linalg.norm(vector, axis=1)
        disM[:, ii] = dis

    # Find the minimun distance for each sensor and assign the field index
    idxList = np.zeros(len(nameList), dtype="int")
    for ii in range(len(idxList)):
        idxList[ii] = np.argmin(disM[ii, :])

    # Print the information
    for ii in range(len(idxList)):
        print(nameList[ii], idxList[ii])
