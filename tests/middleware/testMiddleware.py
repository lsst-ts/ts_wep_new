import time
import unittest

from lsst.ts.wep.middleware.Middleware import Middleware


class TestMiddleware(unittest.TestCase):
    """Test the middleware class."""

    def setUp(self):

        # Module name
        moduleName = "tcsWEP"

        # Declare the Middleware
        self.wepSalIssue = Middleware(moduleName)
        self.wepSalGet = Middleware(moduleName)

    def testTelemetry(self):

        # Reset the topic
        self.wepSalIssue.resetTopic()
        self.wepSalGet.resetTopic()

        # Set the time out
        timeOut = 15
        self.wepSalGet.setTimeOut(timeOut)
        
        # Test the set timeOut
        self.assertEqual(self.wepSalGet.timeOut, timeOut)

        # Set the telemetry topic
        topic = "timestamp"

        # Data information of "tcsWEP_timestamp"
        timestamp = 2.0
        newData = {"timestamp": timestamp}

        # Issue the telemetry
        self.wepSalIssue.issueTelemetry(topic, newData)

        # Sleep 1 sec
        time.sleep(1)

        # Get the telemetry
        self.wepSalGet.getTelemetry(topic)

        # Test to get the telemetry
        self.assertEqual(self.wepSalGet.retData["timestamp"], timestamp)

    def testEvent(self):

        # Reset the topic
        self.wepSalIssue.resetTopic()
        self.wepSalGet.resetTopic()

        # Set the event topic
        topic = "summaryState"

        # Event data
        summaryStateValue = 2
        priority = 1
        newData = {"summaryState": summaryStateValue,
                   "priority": priority}

        # Issue the event
        self.wepSalIssue.issueEvent(topic, newData)

        # Sleep 1 sec
        time.sleep(1)

        # Get the event
        self.wepSalGet.getEvent(topic)

        # Test to get the event information
        self.assertEqual(self.wepSalGet.retData["summaryState"],
                         summaryStateValue)

    def testCommand(self):

        # Reset the topic
        self.wepSalIssue.resetTopic()
        self.wepSalGet.resetTopic()

        # Set the command topic
        topic = "start"

        # Command data
        settingsToApply = "defaultSetting"
        newData = {"settingsToApply": settingsToApply}

        # Issue the event
        self.wepSalIssue.issueCommand(topic, newData)

        # Sleep 1 sec
        self.wepSalGet.getCommand(topic)

        # Test to get the command item list
        self.assertEqual(len(self.wepSalGet.retData), 5)

    def tearDown(self):

        # Turn off the sal
        self.wepSalIssue.shutDownSal()
        self.wepSalGet.shutDownSal()


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
