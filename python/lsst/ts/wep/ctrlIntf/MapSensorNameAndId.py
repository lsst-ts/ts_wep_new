import os

from lsst.ts.wep.ParamReader import ParamReader
from lsst.ts.wep.Utility import getConfigDir


class MapSensorNameAndId(object):

    def __init__(self, sensorNameToIdFileName="sensorNameToId.yaml"):
        """Construct a MapSensorNameAndId object.

        Parameters
        ----------
        sensorNameToIdFileName : str, optional
            Configuration file name to map sensor name and Id. (the default is
            "sensorNameToId.yaml".)
        """

        sensorNameToIdFilePath = os.path.join(getConfigDir(),
                                              sensorNameToIdFileName)
        self._sensorNameToIdFile = ParamReader(filePath=sensorNameToIdFilePath)

    def mapSensorNameToId(self, sensorName):
        """Map the sensor name to sensor Id.

        Parameters
        ----------
        sensorName : list[str] or str
            List or string of abbreviated sensor names.

        Returns
        -------
        list[int]
            List of sensor Id.
        """

        sensorNameList = self._changeToListIfNeed(sensorName)

        sensorIdList = []
        for sensor in sensorNameList:
            sensorId = self._sensorNameToIdFile.getSetting(sensor)
            sensorIdList.append(sensorId)

        return sensorIdList

    def _changeToListIfNeed(self, inputArg):
        """Change the input argument to list type if needed.

        Parameters
        ----------
        inputArg : obj
            Input argument.

        Returns
        -------
        list
            Input argument as the list type.
        """

        if (not isinstance(inputArg, list)):
            inputArg = [inputArg]

        return inputArg

    def mapSensorIdToName(self, sensorId):
        """Map the sensor Id to sensor name.

        If no sensor name is found for a specific Id, there will be no returned
        value.

        Parameters
        ----------
        sensorId : list[int] or int
            List or integer of sensor Id.

        Returns
        -------
        list
            List of abbreviated sensor names.
        int
            Number of sensors.
        """

        sensorIdList = self._changeToListIfNeed(sensorId)

        sensorNameList = []
        content = self._sensorNameToIdFile.getContent()
        for sensor in sensorIdList:
            try:
                sensorName = self._getKeyFromValueInDict(content, sensor)
                sensorNameList.append(sensorName)
            except ValueError:
                pass

        return sensorNameList, len(sensorNameList)

    def _getKeyFromValueInDict(self, aDict, value):
        """Get the key from value in a dictionary object.

        Parameters
        ----------
        aDict : dict
            Dictionary object.
        value : str or int
            Value in the dictionary.

        Returns
        -------
        str
            Dictionary key.
        """

        return list(aDict.keys())[list(aDict.values()).index(value)]


if __name__ == "__main__":
    pass
