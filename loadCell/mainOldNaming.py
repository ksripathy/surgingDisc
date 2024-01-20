import os
import sys
#%matplotlib widget

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

#Add src folder to python path
sys.path.append(srcDir)

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

from src.forces import recordedForces
from src.forces import aeroForces

P45FstreamVelc = np.array([2.988726, 3.057711, 2.998582, 2.050731, 2.068581, 2.031872])
P45FstreamRho = np.array([1.188489, 1.186424, 1.185136, 1.186961, 1.186961, 1.201072])

P70FstreamVelc = np.array([3.045224, 3.030344, 3.035440, 2.025460, 2.049157, 2.007377])
P70FstreamRho = np.array([1.196899, 1.198764, 1.198306, 1.194673, 1.192900, 1.189482])

expSurgeFreq = np.array([2.387324, 4.774648, 4.774648, 4.774648, 3.183099, 4.297183])

noWindObjP45Case0 = recordedForces(dataDir + "/clean/P45k1a5clean.tdms", P45FstreamVelc[0], P45FstreamRho[0], 0.20, expSurgeFreq[0])
windObjP45Case0 = recordedForces(dataDir + "/dynamic/P45k1a5.tdms", P45FstreamVelc[0], P45FstreamRho[0], 0.20, expSurgeFreq[0])
aeroObjP45Case0 = aeroForces(windObjP45Case0, noWindObjP45Case0, 2)

fig1, ax1 = plt.subplots()
plt.plot(aeroObjP45Case0.cycleNdimTime,noWindObjP45Case0.cycleAvgSyncForce)
plt.plot(aeroObjP45Case0.cycleNdimTime,windObjP45Case0.cycleAvgForce)
plt.plot(aeroObjP45Case0.cycleNdimTime,noWindObjP45Case0.cycleAvgSyncForce)
