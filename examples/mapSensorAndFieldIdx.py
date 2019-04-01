import os
import numpy as np

from lsst.ts.wep.SourceProcessor import SourceProcessor
from lsst.ts.wep.Utility import getModulePath


if __name__ == "__main__":

    # Folder path of focal plane txt file
    folderPath2FocalPlane = os.path.join(getModulePath(), "tests", "testData")

    # Configurate the source processor
    sourPro = SourceProcessor(folderPath2FocalPlane=folderPath2FocalPlane)

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
    pointAngle = np.arange(nArm) * (2*np.pi)/nArm
    fieldX = np.concatenate([np.zeros(1), np.kron(armLen, np.cos(pointAngle)),
                            fieldWFSx])
    fieldY = np.concatenate([np.zeros(1), np.kron(armLen, np.sin(pointAngle)),
                            fieldWFSy])

    xyField = np.array([fieldX, fieldY]).T

    # Calculate the distance matrix (sensor by field)
    disM = np.zeros((len(xList), len(fieldX)))
    for ii in range(len(fieldX)):
        vector = xySensor-xyField[ii, :]
        dis = np.linalg.norm(vector, axis=1)
        disM[:, ii] = dis

    # Find the minimun distance for each sensor and assign the field index
    idxList = np.zeros(len(nameList), dtype="int")
    for ii in range(len(idxList)):
        idxList[ii] = np.argmin(disM[ii, :])

    # Print the information
    for ii in range(len(idxList)):
        print(nameList[ii], idxList[ii])
