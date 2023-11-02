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
from src.plot import cyclePlot
from src.miscTools import lagShift

staticP45v2 = recordedForces(dataDir + "/static/P45v2.tdms", 2.052653, 1.194716, 0.2)
staticP45v3 = recordedForces(dataDir + "/static/P45v3.tdms", 3.092307, 1.194394, 0.2)
staticP45v12 = recordedForces(dataDir + "/static/P45v12.tdms", 12.066110, 1.180721, 0.2)
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
staticObjP45 = np.array([staticP45v12, staticP45v12, staticP45v12, staticP45v12, staticP45v12, staticP45v12])
P45FstreamVelc = np.array([2.988726, 3.057711, 2.998582, 2.050731, 2.068581, 2.031872])
P45FstreamRho = np.array([1.188489, 1.186424, 1.185136, 1.186961, 1.186961, 1.201072])

noWindObjP70 = np.empty(6, dtype=np.ndarray)
windObjP70 = np.empty(6, dtype=np.ndarray)
aeroObjP70 = np.empty(6, dtype = np.ndarray)
plotObjP70 = np.empty(6, dtype = np.ndarray)
staticObjP70 = np.array([staticP70v3, staticP70v3, staticP70v3, staticP70v2, staticP70v2, staticP70v2])
P70FstreamVelc = np.array([3.045224, 3.030344, 3.035440, 2.025460, 2.049157, 2.007377])
P70FstreamRho = np.array([1.196899, 1.198764, 1.198306, 1.194673, 1.192900, 1.189482])

expSurgeFreq = np.array([2.387324, 4.774648, 4.774648, 4.774648, 3.183099, 4.297183])

#Data from clean2 and dynamic2 directories is from additional loadCell tests without PIV data
noWindObjP45[0] = recordedForces(dataDir + "/clean/P45k1a5clean.tdms", P45FstreamVelc[0], P45FstreamRho[0], 0.20, expSurgeFreq[0])
noWindObjP45[1] = recordedForces(dataDir + "/clean/P45k2a2p5clean.tdms", P45FstreamVelc[1], P45FstreamRho[1], 0.20, expSurgeFreq[1])
noWindObjP45[2] = recordedForces(dataDir + "/clean/P45k2a5clean.tdms", P45FstreamVelc[2], P45FstreamRho[2], 0.20, expSurgeFreq[2])
noWindObjP45[3] = recordedForces(dataDir + "/clean2/P45k2a7p5clean.tdms", P45FstreamVelc[3], P45FstreamRho[3], 0.20, expSurgeFreq[3])
noWindObjP45[4] = recordedForces(dataDir + "/clean2/P45k3a5clean.tdms", P45FstreamVelc[4], P45FstreamRho[4], 0.20, expSurgeFreq[4])
noWindObjP45[5] = recordedForces(dataDir + "/clean/P45vmax1clean.tdms", P45FstreamVelc[5], P45FstreamRho[5], 0.20, expSurgeFreq[5])

windObjP45[0] = recordedForces(dataDir + "/dynamic/P45k1a5.tdms", P45FstreamVelc[0], P45FstreamRho[0], 0.20, expSurgeFreq[0])
windObjP45[1] = recordedForces(dataDir + "/dynamic/P45k2a2p5.tdms", P45FstreamVelc[1], P45FstreamRho[1], 0.20, expSurgeFreq[1])
windObjP45[2] = recordedForces(dataDir + "/dynamic/P45k2a5.tdms", P45FstreamVelc[2], P45FstreamRho[2], 0.20, expSurgeFreq[2])
windObjP45[3] = recordedForces(dataDir + "/dynamic2/P45k2a7p5.tdms", P45FstreamVelc[3], P45FstreamRho[3], 0.20, expSurgeFreq[3])
windObjP45[4] = recordedForces(dataDir + "/dynamic2/P45k3a5.tdms", P45FstreamVelc[4], P45FstreamRho[4], 0.20, expSurgeFreq[4])
windObjP45[5] = recordedForces(dataDir + "/dynamic/P45vmax1.tdms", P45FstreamVelc[5], P45FstreamRho[5], 0.20, expSurgeFreq[5])

noWindObjP70[0] = recordedForces(dataDir + "/clean/P70k1a5clean.tdms", P70FstreamVelc[0], P70FstreamRho[0], 0.20, expSurgeFreq[0])
noWindObjP70[1] = recordedForces(dataDir + "/clean/P70k2a2p5clean.tdms", P70FstreamVelc[1], P70FstreamRho[1], 0.20, expSurgeFreq[1])
noWindObjP70[2] = recordedForces(dataDir + "/clean/P70k2a5clean.tdms", P70FstreamVelc[2], P70FstreamRho[2], 0.20, expSurgeFreq[2])
noWindObjP70[3] = recordedForces(dataDir + "/clean/P70k2a7p5clean.tdms", P70FstreamVelc[3], P70FstreamRho[3], 0.20, expSurgeFreq[3])
noWindObjP70[4] = recordedForces(dataDir + "/clean/P70k3a5clean.tdms", P70FstreamVelc[4], P70FstreamRho[4], 0.20, expSurgeFreq[4])
noWindObjP70[5] = recordedForces(dataDir + "/clean/P70vmax1clean.tdms", P70FstreamVelc[5], P70FstreamRho[5], 0.20, expSurgeFreq[5])

windObjP70[0] = recordedForces(dataDir + "/dynamic/P70k1a5.tdms", P70FstreamVelc[0], P70FstreamRho[0], 0.20, expSurgeFreq[0])
windObjP70[1] = recordedForces(dataDir + "/dynamic/P70k2a2p5.tdms", P70FstreamVelc[1], P70FstreamRho[1], 0.20, expSurgeFreq[1])
windObjP70[2] = recordedForces(dataDir + "/dynamic/P70k2a5.tdms", P70FstreamVelc[2], P70FstreamRho[2], 0.20,  expSurgeFreq[2])
windObjP70[3] = recordedForces(dataDir + "/dynamic/P70k2a7p5.tdms", P70FstreamVelc[3], P70FstreamRho[3], 0.20, expSurgeFreq[3])
windObjP70[4] = recordedForces(dataDir + "/dynamic/P70k3a5.tdms", P70FstreamVelc[4], P70FstreamRho[4], 0.20, expSurgeFreq[4])
windObjP70[5] = recordedForces(dataDir + "/dynamic/P70vmax1.tdms", P70FstreamVelc[5], P70FstreamRho[5], 0.20, expSurgeFreq[5])

#Generate aero and plotObj
for i in range(6):
    
    aeroObjP45[i] = aeroForces(windObjP45[i], noWindObjP45[i], 2)
    aeroObjP70[i] = aeroForces(windObjP70[i], noWindObjP70[i], 2)
    
    #print(windObjP45[i].uid + f"Mean inertial loading: {np.mean(noWindObjP45[i].syncForceTrim)}")
    #print(windObjP70[i].uid + f"Mean inertial loading: {np.mean(noWindObjP70[i].syncForceTrim)}")
    
    #print(windObjP45[i].uid + f"Mean filtered inertial loading: {np.mean(noWindObjP45[i].syncFiltForceTrim)}")
    #print(windObjP70[i].uid + f"Mean filtered inertial loading: {np.mean(noWindObjP70[i].syncFiltForceTrim)}")
    
    print(windObjP45[i].uid + f"Mean aero loading: {np.mean(aeroObjP45[i].forceSeries)}, {np.std(aeroObjP45[i].forceSeries)}, {0.5*(max(aeroObjP45[i].cycleAvgFiltForce) - min(aeroObjP45[i].cycleAvgFiltForce))}")
    print(windObjP70[i].uid + f"Mean aero loading: {np.mean(aeroObjP70[i].forceSeries)}, {np.std(aeroObjP70[i].forceSeries)}, {0.5*(max(aeroObjP70[i].cycleAvgFiltForce) - min(aeroObjP45[i].cycleAvgFiltForce))}")

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
    
P45Axs.hlines(np.mean(staticP45v12.forceSeries), xmin = 0, xmax = 1, color='k', label="staticAerov2")
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

#Final plots
P45ConstFreqFig, P45ConstFreqAxs = plt.subplots()
P45ConstFreqAxs.set_xlabel("Time [-]")
P45ConstFreqAxs.set_ylabel(r"$C_T$ [-]")
P45ConstFreqAxs.set_xlim(0,1)

P45ConstAmpFig, P45ConstAmpAxs = plt.subplots()
P45ConstAmpAxs.set_xlabel(r"Time [-]")
P45ConstAmpAxs.set_ylabel(r"$C_T$ [-]")
P45ConstAmpAxs.set_xlim(0,1)

P70ConstFreqFig, P70ConstFreqAxs = plt.subplots()
P70ConstFreqAxs.set_xlabel(r"Time [-]")
P70ConstFreqAxs.set_ylabel(r"$C_T$ [-]")
P70ConstFreqAxs.set_xlim(0,1)

P70ConstAmpFig, P70ConstAmpAxs = plt.subplots()
P70ConstAmpAxs.set_xlabel(r"Time [-]")
P70ConstAmpAxs.set_ylabel(r"$C_T$ [-]")
P70ConstAmpAxs.set_xlim(0,1)

P45ConstFreqData = [aeroObjP45[1], aeroObjP45[2], aeroObjP45[3]]
P70ConstFreqData = [aeroObjP70[1], aeroObjP70[2], aeroObjP70[3]]

P45ConstAmpData = [aeroObjP45[0], aeroObjP45[2], aeroObjP45[4]]
P70ConstAmpData = [aeroObjP70[0], aeroObjP70[2], aeroObjP70[4]]

constFreqLabel = [r"$\frac{x_{a}}{D} = 0.125$",r"$\frac{x_{a}}{D} = 0.25$", r"$\frac{x_{a}}{D} = 0.375$"]
constAmpLabel = [r"$k = 1$", r"$k = 2$", r"$k = 3$"]

plotLineStyles = ["-", "-", "--"]

P45ConstAmpAxs.hlines(np.mean(staticP45v3.forceSeries), xmin = 0, xmax = 1, color='k', label=r"$k = 0, V_\infty = 3 \ ms^{-1}$")
P45ConstAmpAxs.hlines(np.mean(staticP45v2.forceSeries), xmin = 0, xmax = 1, color='k', linestyle="--", label=r"$k = 0, V_\infty = 2 \ ms^{-1}$")
P70ConstAmpAxs.hlines(np.mean(staticP70v3.forceSeries), xmin = 0, xmax = 1, color='k', label=r"$k = 0, V_\infty = 3 \ ms^{-1}$")
P70ConstAmpAxs.hlines(np.mean(staticP70v2.forceSeries), xmin = 0, xmax = 1, color='k', linestyle="--", label=r"$k = 0, V_\infty = 2 \ ms^{-1}$")

P45ConstFreqAxs.hlines(np.mean(staticP45v3.forceSeries), xmin = 0, xmax = 1, color='k', label = r"$\frac{x_{a}}{D} = 0, V_\infty = 3 \ ms^{-1}$")
P45ConstFreqAxs.hlines(np.mean(staticP45v2.forceSeries), xmin = 0, xmax = 1, color='k', linestyle="--", label=r"$\frac{x_{a}}{D} = 0, V_\infty = 2 \ ms^{-1}$")
P70ConstFreqAxs.hlines(np.mean(staticP70v3.forceSeries), xmin = 0, xmax = 1, color='k', label = r"$\frac{x_{a}}{D} = 0, V_\infty = 3 \ ms^{-1}$")
P70ConstFreqAxs.hlines(np.mean(staticP70v2.forceSeries), xmin = 0, xmax = 1, color='k', linestyle="--", label=r"$\frac{x_{a}}{D} = 0, V_\infty = 2 \ ms^{-1}$")

for i in range(3):
            
    P45ConstFreqAxs.plot(P45ConstFreqData[i].cycleNdimTime, P45ConstFreqData[i].cycleAvgFiltForce, color=plotColors[i], linestyle=plotLineStyles[i], label=constFreqLabel[i])
    P45ConstAmpAxs.plot(P45ConstAmpData[i].cycleNdimTime, P45ConstAmpData[i].cycleAvgFiltForce, color=plotColors[i], linestyle=plotLineStyles[i], label=constAmpLabel[i])
    
    P70ConstFreqAxs.plot(P70ConstFreqData[i].cycleNdimTime, P70ConstFreqData[i].cycleAvgFiltForce, color=plotColors[i], linestyle=plotLineStyles[i], label=constFreqLabel[i])
    P70ConstAmpAxs.plot(P70ConstAmpData[i].cycleNdimTime, P70ConstAmpData[i].cycleAvgFiltForce, color=plotColors[i], linestyle=plotLineStyles[i], label=constAmpLabel[i])

P45ConstFreqAxs.grid()
P45ConstAmpAxs.grid()
P70ConstFreqAxs.grid()
P70ConstAmpAxs.grid()

P45ConstFreqAxs.legend()
P45ConstAmpAxs.legend()
P70ConstFreqAxs.legend()
P70ConstAmpAxs.legend()

#P45ConstFreqFig.savefig(plotDir + "/P45plots/P45constFreq.png", dpi=300)
#P45ConstAmpFig.savefig(plotDir + "/P45plots/P45constAmp.png", dpi=300)

plt.show()



#plotObj = cyclePlot(windObj, noWindObj, aeroObj, staticP70v3)
#nims = plotObj.genAnim(plotDir)
#plotObj.genPlot(25,plotDir)


