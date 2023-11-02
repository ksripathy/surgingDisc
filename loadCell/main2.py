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

from src.forces import recordedForces
from src.forces import aeroForces
from src.plot import cyclePlot

staticP45v2 = recordedForces(dataDir + "/static/P45v2.tdms", 2.052653, 1.194716, 0.2)
staticP45v3 = recordedForces(dataDir + "/static/P45v3.tdms", 2.999658, 1.181605, 0.2)
staticP45v5 = recordedForces(dataDir + "/static/P45v5.tdms", 5.015109, 1.182036, 0.2)
staticP45v8 = recordedForces(dataDir + "/static/P45v8.tdms", 7.966258, 1.181688, 0.2)
staticP45v12 = recordedForces(dataDir + "/static/P45v12.tdms", 12.066110, 1.180721, 0.2)
staticP45v16 = recordedForces(dataDir + "/static/P45v16.tdms", 15.992357, 1.177549, 0.2)
staticP70v2 = recordedForces(dataDir + "/static/P70v2.tdms", 2.045414, 1.195081, 0.2)
staticP70v3 = recordedForces(dataDir + "/static/P70v3.tdms", 3.096507, 1.194618, 0.2)

avgP45v2 = np.mean(staticP45v2.forceSeries)
avgP45v3 = np.mean(staticP45v3.forceSeries)
avgP70v2 = np.mean(staticP70v2.forceSeries)
avgP70v3 = np.mean(staticP70v3.forceSeries)

#expSurgeFreq = 2.387324

#noWindObj = recordedForces(dataDir + "/clean/P70k1a5clean.tdms", expSurgeFreq)
#windObj = recordedForces(dataDir + "/dynamic/P70k1a5.tdms", expSurgeFreq)
#aeroObj = aeroForces(windObj, noWindObj, 2)

noWindObjP45 = np.empty(6, dtype=np.ndarray)
windObjP45 = np.empty(6, dtype = np.ndarray)
aeroObjP45 = np.empty(6, dtype = np.ndarray)
plotObjP45 = np.empty(6, dtype = np.ndarray)
staticObjP45 = np.array([staticP45v2, staticP45v2, staticP45v2, staticP45v2, staticP45v2, staticP45v2])
P45FstreamVelc = np.array([2.064786, 2.057820, 2.046670, 2.060026, 2.068581, 2.059627])
P45FstreamRho = np.array([1.187888, 1.187394, 1.186354, 1.186983, 1.186961, 1.186815])

noWindObjP70 = np.empty(6, dtype=np.ndarray)
windObjP70 = np.empty(6, dtype=np.ndarray)
aeroObjP70 = np.empty(6, dtype = np.ndarray)
plotObjP70 = np.empty(6, dtype = np.ndarray)
staticObjP70 = np.array([staticP70v3, staticP70v3, staticP70v3, staticP70v2, staticP70v2, staticP70v2])
P70FstreamVelc = np.array([3.045224, 3.030344, 3.035440, 2.025460, 2.049157, 2.007377])
P70FstreamRho = np.array([1.196899, 1.198764, 1.198306, 1.194673, 1.192900, 1.189482])

expSurgeFreqP70 = np.array([2.387324, 4.774648, 4.774648, 4.774648, 3.183099, 4.297183])
expSurgeFreqP45 = np.array([1.273240, 1.273240, 1.273240, 3.183099, 4.774648, 6.366198])

#Data from clean2 and dynamic2 directories is from additional loadCell tests without PIV data.
#Comparison is made between non standard surge motion scenarios
noWindObjP45[0] = recordedForces(dataDir + "/clean2/P45k1a2p5clean.tdms", P45FstreamVelc[0], P45FstreamRho[0], 0.20, expSurgeFreqP45[0])
noWindObjP45[1] = recordedForces(dataDir + "/clean2/P45k1a5clean.tdms", P45FstreamVelc[1], P45FstreamRho[1], 0.20, expSurgeFreqP45[1])
noWindObjP45[2] = recordedForces(dataDir + "/clean2/P45k1a7p5clean.tdms", P45FstreamVelc[2], P45FstreamRho[2], 0.20, expSurgeFreqP45[2])
noWindObjP45[3] = recordedForces(dataDir + "/clean2/P45k2a5clean.tdms", P45FstreamVelc[3], P45FstreamRho[3], 0.20, expSurgeFreqP45[3])
noWindObjP45[4] = recordedForces(dataDir + "/clean2/P45k3a5clean.tdms", P45FstreamVelc[4], P45FstreamRho[4], 0.20, expSurgeFreqP45[4])
noWindObjP45[5] = recordedForces(dataDir + "/clean2/P45k4a5clean.tdms", P45FstreamVelc[5], P45FstreamRho[5], 0.20, expSurgeFreqP45[5])

windObjP45[0] = recordedForces(dataDir + "/dynamic2/P45k1a2p5.tdms", P45FstreamVelc[0], P45FstreamRho[0], 0.20, expSurgeFreqP45[0])
windObjP45[1] = recordedForces(dataDir + "/dynamic2/P45k1a5.tdms", P45FstreamVelc[1], P45FstreamRho[1], 0.20, expSurgeFreqP45[1])
windObjP45[2] = recordedForces(dataDir + "/dynamic2/P45k1a7p5.tdms", P45FstreamVelc[2], P45FstreamRho[2], 0.20, expSurgeFreqP45[2])
windObjP45[3] = recordedForces(dataDir + "/dynamic2/P45k2a5.tdms", P45FstreamVelc[3], P45FstreamRho[3], 0.20, expSurgeFreqP45[3])
windObjP45[4] = recordedForces(dataDir + "/dynamic2/P45k3a5.tdms", P45FstreamVelc[4], P45FstreamRho[4], 0.20, expSurgeFreqP45[4])
windObjP45[5] = recordedForces(dataDir + "/dynamic2/P45k4a5.tdms", P45FstreamVelc[5], P45FstreamRho[5], 0.20, expSurgeFreqP45[5])

noWindObjP70[0] = recordedForces(dataDir + "/clean/P70k1a5clean.tdms", P70FstreamVelc[0], P70FstreamRho[0], 0.20, expSurgeFreqP70[0])
noWindObjP70[1] = recordedForces(dataDir + "/clean/P70k2a2p5clean.tdms", P70FstreamVelc[1], P70FstreamRho[1], 0.20, expSurgeFreqP70[1])
noWindObjP70[2] = recordedForces(dataDir + "/clean/P70k2a5clean.tdms", P70FstreamVelc[2], P70FstreamRho[2], 0.20, expSurgeFreqP70[2])
noWindObjP70[3] = recordedForces(dataDir + "/clean/P70k2a7p5clean.tdms", P70FstreamVelc[3], P70FstreamRho[3], 0.20, expSurgeFreqP70[3])
noWindObjP70[4] = recordedForces(dataDir + "/clean/P70k3a5clean.tdms", P70FstreamVelc[4], P70FstreamRho[4], 0.20, expSurgeFreqP70[4])
noWindObjP70[5] = recordedForces(dataDir + "/clean/P70vmax1clean.tdms", P70FstreamVelc[5], P70FstreamRho[5], 0.20, expSurgeFreqP70[5])

windObjP70[0] = recordedForces(dataDir + "/dynamic/P70k1a5.tdms", P70FstreamVelc[0], P70FstreamRho[0], 0.20, expSurgeFreqP70[0])
windObjP70[1] = recordedForces(dataDir + "/dynamic/P70k2a2p5.tdms", P70FstreamVelc[1], P70FstreamRho[1], 0.20, expSurgeFreqP70[1])
windObjP70[2] = recordedForces(dataDir + "/dynamic/P70k2a5.tdms", P70FstreamVelc[2], P70FstreamRho[2], 0.20,  expSurgeFreqP70[2])
windObjP70[3] = recordedForces(dataDir + "/dynamic/P70k2a7p5.tdms", P70FstreamVelc[3], P70FstreamRho[3], 0.20, expSurgeFreqP70[3])
windObjP70[4] = recordedForces(dataDir + "/dynamic/P70k3a5.tdms", P70FstreamVelc[4], P70FstreamRho[4], 0.20, expSurgeFreqP70[4])
windObjP70[5] = recordedForces(dataDir + "/dynamic/P70vmax1.tdms", P70FstreamVelc[5], P70FstreamRho[5], 0.20, expSurgeFreqP70[5])

#Generate aero and plotObj
for i in range(6):
    
    aeroObjP45[i] = aeroForces(windObjP45[i], noWindObjP45[i], 1)
    aeroObjP70[i] = aeroForces(windObjP70[i], noWindObjP70[i], 1)
    
    #print(windObjP45[i].uid + f"Mean inertial loading: {np.mean(noWindObjP45[i].syncForceTrim)}")
    #print(windObjP70[i].uid + f"Mean inertial loading: {np.mean(noWindObjP70[i].syncForceTrim)}")
    
    #print(windObjP45[i].uid + f"Mean filtered inertial loading: {np.mean(noWindObjP45[i].syncFiltForceTrim)}")
    #print(windObjP70[i].uid + f"Mean filtered inertial loading: {np.mean(noWindObjP70[i].syncFiltForceTrim)}")
    
    print(windObjP45[i].uid + f"Mean aero loading: {np.mean(aeroObjP45[i].forceSeries)}")
    #print(windObjP70[i].uid + f"Mean aero loading: {np.mean(aeroObjP70[i].filtForce)}")

    #plotObjP45[i] = cyclePlot(windObjP45[i], noWindObjP45[i], aeroObjP45[i], staticObjP45[i])
    #plotObjP70[i] = cyclePlot(windObjP70[i], noWindObjP70[i], aeroObjP70[i], staticObjP70[i])
    
    #animP45 = plotObjP45[i].genAnim(plotDir+"/P45anims")
    #animP70 = plotObjP70[i].genAnim(plotDir+"/P70anims")
    
    #plotP45 = plotObjP45[i].genAnim(plotDir+"/P45plots")
    #plotP70 = plotObjP70[i].genAnim(plotDir+"/P70plots")

plotColors = ["tab:red", "tab:green", "tab:blue", "tab:cyan", "tab:purple", "tab:brown"]

P45Fig, P45Axs = plt.subplots()
P70Fig, P70Axs = plt.subplots()

for i in range(6):
    
    P45Axs.plot(aeroObjP45[i].cycleNdimTime, aeroObjP45[i].cycleAvgFiltForce, color = plotColors[i], label=windObjP45[i].uid)
    P70Axs.plot(aeroObjP70[i].cycleNdimTime, aeroObjP70[i].cycleAvgFiltForce, color = plotColors[i], label=windObjP70[i].uid)
    #P45Axs.plot(aeroObjP45[i].cycleNdimTime, aeroObjP45[i].cycleAvgFiltForceV3, color = plotColors[i], linestyle = "--")
    #P70Axs.plot(aeroObjP70[i].cycleNdimTime, aeroObjP70[i].cycleAvgFiltForceV3, color = plotColors[i], linestyle = "--")
    
P45Axs.hlines(np.mean(staticP45v2.forceSeries), xmin = 0, xmax = 1, color='k', label="staticAerov2")
#P45Axs.hlines(np.mean(staticP45v3.forceSeries), xmin = 0, xmax = 1, color='k', linestyle="--", label="staticAerov3")
P70Axs.hlines(np.mean(staticP70v2.forceSeries), xmin = 0, xmax = 1, color='k', label="staticAerov2")
P70Axs.hlines(np.mean(staticP70v3.forceSeries), xmin = 0, xmax = 1, color='k', linestyle="--", label="staticAerov3")

P45Axs.set_xlabel("Time [-]")
P45Axs.set_ylabel(F"$C_T$ [-]")
P70Axs.set_xlabel("Time [-]")
P70Axs.set_ylabel(F"$C_T$ [-]")

P45Axs.grid()
P45Axs.legend()

P70Axs.grid()
P70Axs.legend()

plt.show()



#plotObj = cyclePlot(windObj, noWindObj, aeroObj, staticP70v3)
#nims = plotObj.genAnim(plotDir)
#plotObj.genPlot(25,plotDir)


