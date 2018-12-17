import time
import numpy as np

from lsst.ts.wep.WEPController import WEPController

if __name__ == "__main__":

    # This script is to show the communication of WEP controller by SAL
    
    # Initiate the WEP Controller
    wepCntlr = WEPController()

    # Set the middle ware
    topicList = ["normalTargetWfsList", "wavefrontError"]
    wepCntlr.setMiddleWare(topicList)

    # Issue the event
    sensorIdList = [1, 2, 3, -1]
    timestamp = time.time()
    priority = 1
    eventData = {"sensorIdList": sensorIdList,
                 "timestamp": timestamp,
                 "priority": priority}
    wepCntlr.issueEvent(topicList[0], eventData)
    time.sleep(1)

    # Issue the event 
    sensorId = 11
    zkList = np.random.rand(19)
    telData = {"sensorId": sensorId,
               "annularZernikePoly": zkList,
               "timestamp": timestamp, 
               "priority": priority}
    wepCntlr.issueEvent(topicList[1], telData)
    time.sleep(1)

    # Turn off the SAL
    wepCntlr.shutDownSal()