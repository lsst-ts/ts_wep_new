import os
import numpy as np
import yaml
import warnings


class ParamReader(object):

    def __init__(self, filePath=None):
        """Initialization of parameter reader of yaml format class.

        Parameters
        ----------
        filePath : str, optional
            File path. (the default is None.)
        """

        if (filePath is None):
            self.filePath = ""
        else:
            self.filePath = filePath

        self._content = self._readContent(self.filePath)

    def _readContent(self, filePath):
        """Read the content of file.

        Parameters
        ----------
        filePath : str
            File path.

        Returns
        -------
        dict
            Content of file.
        """

        if (os.path.exists(filePath)):
            with open(filePath, "r") as yamlFile:
                content = yaml.safe_load(yamlFile)
        else:
            content = dict()

        return content

    def getFilePath(self):
        """Get the parameter file path.

        Returns
        -------
        str
            Get the file path.
        """

        return self.filePath

    def setFilePath(self, filePath):
        """Set the file path and load the setting.

        Parameters
        ----------
        filePath : str
            File path.
        """

        self.__init__(filePath=filePath)

    def getContent(self):
        """Get the content.

        Returns
        -------
        list or dict
            Content.
        """

        return self._content

    def getSetting(self, param):
        """Get the setting value.

        Parameters
        ----------
        param : str
            Parameter name.

        Returns
        -------
        int, float, list, or dict
            Parameter value.

        Raises
        ------
        ValueError
            The parameter does not exist.
        """

        try:
            return self._content[param]
        except KeyError:
            raise ValueError("The '%s' does not exist." % param)

    @staticmethod
    def writeMatToFile(matrix, filePath):
        """Write the matrix data to file.

        Parameters
        ----------
        matrix : numpy.ndarray
            Matrix data.
        filePath : str
            Yaml file path.
        """

        ParamReader._writeDataToFile(matrix.tolist(), filePath)

    @staticmethod
    def _writeDataToFile(data, filePath):
        """Write the data to file.

        Parameters
        ----------
        data : Python object
            Data.
        filePath : str
            Yaml file path.

        Raises
        ------
        ValueError
            The file name should end with ".yaml".
        """

        if filePath.endswith(".yaml"):
            with open(filePath, "w") as yamlFile:
                yaml.safe_dump(data, stream=yamlFile, default_flow_style=None)
        else:
            raise ValueError("The file name should end with '.yaml'.")

    def getMatContent(self):
        """Get the matrix content.

        Returns
        -------
        numpy.ndarray
            Matrix content.
        """

        if (self._content == dict()):
            mat = np.array([])
        else:
            mat = np.array(self._content)

        return mat

    def updateSettingSeries(self, settingSeries):
        """Update the settings based on a serious of setting.

        Parameters
        ----------
        settingSeries : dict
            A serious of setting to update.
        """

        for param, value in settingSeries.items():
            self.updateSetting(param, value)

    def updateSetting(self, param, value):
        """Update the setting.

        Parameters
        ----------
        param : str
            Parameter name.
        value : int, float, list, or dict
            Updated value.
        """

        origVal = self.getSetting(param)
        if not isinstance(value, type(origVal)):
            warnings.warn("Update with the different type of value.",
                          category=UserWarning)

        self._content[param] = value

    def saveSetting(self, filePath=None):
        """Save the setting.

        Parameters
        ----------
        filePath : str, optional
            File path. If None, the loaded file path will be used. (the
            default is None.)
        """

        if (filePath is None):
            filePath = self.filePath
        else:
            self.filePath = filePath

        self._writeDataToFile(self._content, filePath)


if __name__ == "__main__":
    pass
