import numpy as np

from lsst.ts.wep.Middleware import Middleware


class MiddlewareWrapper(object):

    def __init__(self):

        self.middleWare = dict()
        self.topicList = dict()

    def setMiddleWare(self, topicList, moduleName="tcsWEP"):
        """
        
        Set the middle ware. 
        
        Arguments:
            topicList {[list]} -- List of topic.
        
        Keyword Arguments:
            moduleName {[str]} -- Module name. (default: {"tcsWEP"})
        """

        if moduleName not in list(self.topicList.keys()):
            self.topicList[moduleName] = np.array(topicList)
            middleWareList = []
            for ii in range(len(topicList)):
                middleWareList.append(Middleware(moduleName))
            self.middleWare[moduleName] = middleWareList
        else:
            print("The module: '%s' is in the list ready." % moduleName)

    def shutDownSal(self):
        """
        
        Shutdown the SAL.
        """

        for moduleName, middleWareList in self.middleWare.items():
            print("Shutdown SAL module: %s" % moduleName)
            for middleWare in middleWareList:
                middleWare.shutDownSal()

    def issueEvent(self, topic, newData, moduleName="tcsWEP"):
        """
        
        Issue the event.
        
        Arguments:
            topic {[str]} -- Topic name.
            newData {[dict]} -- New input data to issue.
        
        Keyword Arguments:
            moduleName {[str]} -- Module name. (default: {"tcsWEP"})
        """

        idx = self.__getIdxOfTopic(topic, moduleName)
        self.middleWare[moduleName][idx].issueEvent(topic, newData)

    def issueTelemetry(self, topic, newData, moduleName="tcsWEP"):
        """
        
        Issue the telemetry.
        
        Arguments:
            topic {[str]} -- Topic name.
            newData {[dict]} -- New input data to issue.
        
        Keyword Arguments:
            moduleName {[str]} -- Module name. (default: {"tcsWEP"})
        """

        idx = self.__getIdxOfTopic(topic, moduleName)
        self.middleWare[moduleName][idx].issueTelemetry(topic, newData)

    def __getIdxOfTopic(self, topic, moduleName):
        """
        
        Get the index of specific topic in specific module. 
        
        Arguments:
            topic {[str]} -- Topic name.
            moduleName {[str]} -- Module name.
        
        Returns:
            [int] -- Index of topic.
        """

        idx = np.where(self.topicList[moduleName] == topic)[0][0]

        return idx


if __name__ == "__main__":
    pass
