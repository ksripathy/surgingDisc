import os
import sys

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

#Add src folder to python path
sys.path.append(srcDir)

import numpy as np

from loadCell.src.forcesOld import forces
from loadCell.src.miscToolsOld import freqsSignal

noWindObj = forces(dataDir + "/clean/P45k2a5clean.tdms")
windObj = forces(dataDir + "/dynamic/P45k2a5.tdms")

freqsNoWind, ampNoWind = freqsSignal(noWindObj.forceSeries, noWindObj.samplingFreq)
sigFreqNoWind = freqsNoWind[np.argmax(np.abs(ampNoWind[1:]))+1]

freqsWind, ampWind = freqsSignal(windObj.forceSeries, windObj.samplingFreq)
sigFreqWind = freqsWind[np.argmax(np.abs(ampWind[1:]))+1]