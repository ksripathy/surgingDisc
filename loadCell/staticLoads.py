import os
import sys

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data/newNaming")

#Add src folder to python path
sys.path.append(srcDir)

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_1samp

from src.ioUtils import caseData
from src.forces import recordedForces

p70Case6 = caseData(dataDir+"/p70Case6.tdms")
p70Case7 = caseData(dataDir+"/p70Case7.tdms")

force = p70Case7.windObj.forceSeries
meanForce = np.mean(force)
res = ttest_1samp(force, meanForce)



