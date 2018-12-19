import unittest

from lsst.ts.wep.Middleware import Middleware
from lsst.ts.wep.MiddlewareWrapper import MiddlewareWrapper


class TestMiddlewareWrapper(unittest.TestCase):
    """Test the middleware wrapper class."""

    def setUp(self):
        
        self.topicList = ["WavefrontErrorCalculated", "WavefrontError"]
        self.moduleName = "tcsWEP"

        # Instantiate the WEP controller
        self.wepCntlr = MiddlewareWrapper()

        # Set another middleware client
        self.middlewareClient = Middleware(self.moduleName)

    def testSalFunction(self):

        # Set the middle ware
        self.wepCntlr.setMiddleWare(self.topicList,
                                    moduleName=self.moduleName)
        self.assertEqual(list(self.wepCntlr.topicList.keys())[0],
                         self.moduleName)
        self.assertEqual(len(self.wepCntlr.topicList[self.moduleName]), 2)

        # Issue the event
        sensorName = "R:2,2 S:1,1"
        timestamp = time.time()
        priority = 1
        eventData = {"sensorID": sensorName,
                     "timestamp": timestamp,
                     "priority": priority}
        self.wepCntlr.issueEvent("WavefrontErrorCalculated", eventData)

        # Issue the telemetry
        zkList = np.random.rand(19)
        telData = {"sensorID": sensorName,
                   "annularZerikePolynomials": zkList,
                   "timestamp": timestamp}
        self.wepCntlr.issueTelemetry("WavefrontError", telData)

        # Sleep 1 second
        time.sleep(1)

        # Accept the event
        self.middlewareClient.getEvent("WavefrontErrorCalculated")

        # Check the value
        self.assertEqual(self.middlewareClient.retData["sensorID"], sensorName)

        # Accept the telemetry
        self.middlewareClient.getEvent("WavefrontErrorCalculated")

        # Reset the topic
        self.middlewareClient.resetTopic()
        self.middlewareClient.getTelemetry("WavefrontError")

        # Check the value
        self.assertAlmostEqual(
            np.sum(self.middlewareClient.retData["annularZerikePolynomials"]), 
            np.sum(zkList))

    def tearDown(self):

        # Turn off the sal
        self.wepCntlr.shutDownSal()
        self.middlewareClient.shutDownSal()


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
