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

from src.ioUtils import caseData
from src.forces import recordedForces

p45Case0 = caseData(dataDir+"/p45Case0.tdms")
p45Case1 = caseData(dataDir+"/p45Case1.tdms")
p45Case2 = caseData(dataDir+"/p45Case2.tdms")
p45Case3 = caseData(dataDir+"/p45Case3.tdms")
p45Case4 = caseData(dataDir+"/p45Case4.tdms")
p45Case5 = caseData(dataDir+"/p45Case5.tdms")
p45Case6 = caseData(dataDir+"/p45Case6.tdms")
p45Case7 = caseData(dataDir+"/p45Case7.tdms")

p70Case0 = caseData(dataDir+"/p70Case0.tdms")
p70Case1 = caseData(dataDir+"/p70Case1.tdms")
p70Case2 = caseData(dataDir+"/p70Case2.tdms")
p70Case3 = caseData(dataDir+"/p70Case3.tdms")
p70Case4 = caseData(dataDir+"/p70Case4.tdms")
p70Case5 = caseData(dataDir+"/p70Case5.tdms")
p70Case6 = caseData(dataDir+"/p70Case6.tdms")
p70Case7 = caseData(dataDir+"/p70Case7.tdms")

plotData2 = "aeroCt"

fig, ax = plt.subplots(ncols=2)
fig.set_size_inches(26/2.54, 10/2.54)
ax[0].plot(p45Case0.aeroObj.cycleNdimTime, getattr(p45Case0, plotData2), label="case0")
ax[0].plot(p45Case1.aeroObj.cycleNdimTime, getattr(p45Case1, plotData2), label="case1")
ax[0].plot(p45Case2.aeroObj.cycleNdimTime, getattr(p45Case2, plotData2), label="case2")
ax[0].plot(p45Case3.aeroObj.cycleNdimTime, getattr(p45Case3, plotData2), label="case3")
ax[0].plot(p45Case4.aeroObj.cycleNdimTime, getattr(p45Case4, plotData2), label="case4")
ax[0].plot(p45Case5.aeroObj.cycleNdimTime, getattr(p45Case5, plotData2), label="case5")
ax[0].axhline(p45Case6.aeroCt, label="case6", color = "k")
ax[0].axhline(p45Case7.aeroCt, label="case7", color = "k", linestyle = "--")
ax[0].grid()

plt.figlegend(loc="upper center", ncols=4, bbox_to_anchor=(0.5,1.05))

ax[1].plot(p70Case0.aeroObj.cycleNdimTime, getattr(p70Case0, plotData2), label="case0")
ax[1].plot(p70Case1.aeroObj.cycleNdimTime, getattr(p70Case1, plotData2), label="case1")
ax[1].plot(p70Case2.aeroObj.cycleNdimTime, getattr(p70Case2, plotData2), label="case2")
ax[1].plot(p70Case3.aeroObj.cycleNdimTime, getattr(p70Case3, plotData2), label="case3")
ax[1].plot(p70Case4.aeroObj.cycleNdimTime, getattr(p70Case4, plotData2), label="case4")
ax[1].plot(p70Case5.aeroObj.cycleNdimTime, getattr(p70Case5, plotData2), label="case5")
ax[1].axhline(p70Case6.aeroCt, label="case6", color = "k")
ax[1].axhline(p70Case7.aeroCt, label="case7", linestyle = "--")
ax[1].grid()

plotData = "cycleAvgFiltForce"
loadObj = "aeroObj"

fig1, ax1 = plt.subplots()
ax1.plot(p45Case0.aeroObj.cycleNdimTime,getattr(getattr(p45Case0,loadObj),plotData),label="case0")
ax1.plot(p45Case1.aeroObj.cycleNdimTime,getattr(getattr(p45Case1,loadObj),plotData),label="case1")
#ax1.plot(p45Case2.aeroObj.cycleNdimTime,getattr(getattr(p45Case2,loadObj),plotData),label="case2")
ax1.plot(p45Case3.aeroObj.cycleNdimTime,getattr(getattr(p45Case3,loadObj),plotData),label="case3")
#ax1.plot(p45Case4.aeroObj.cycleNdimTime,getattr(getattr(p45Case4,loadObj),plotData),label="case4")
#ax1.plot(p45Case5.aeroObj.cycleNdimTime,getattr(getattr(p45Case5,loadObj),plotData),label="case5")
#ax1.axhline(np.mean(p45Case6.windObj.forceSeries))
#ax1.axhline(np.mean(p45Case7.windObj.forceSeries))
plt.legend()

'''fig2, ax2 = plt.subplots()
ax2.plot(p45Case5.aeroObj.cycleNdimTime,p45Case5.windObj.cycleAvgFiltForce)
ax2.plot(p45Case5.aeroObj.cycleNdimTime,p45Case5.noWindObj.cycleAvgSyncForce)
ax2.plot(p45Case5.aeroObj.cycleNdimTime,p45Case5.aeroObj.cycleAvgFiltForce)

from src.plot import cyclePlot

plotObj = cyclePlot(p45Case4.windObj, p45Case4.noWindObj, p45Case4.aeroObj, p45Case7.windObj)
anims = plotObj.genAnim(plotDir)
plotObj.genPlot(25,plotDir)'''

#testForceObj = recordedForces(dataDir+"/p45Case0.tdms",1,1,0.2,1)

obj = p45Case5
fig3, ax3 = plt.subplots()
#ax3.plot(obj.aeroObj.cycleNdimTime, obj.windObj.cycleAvgFiltForce)
#ax3.plot(obj.aeroObj.cycleNdimTime, obj.noWindObj.cycleAvgSyncFiltForce)
ax3.plot(obj.aeroObj.cycleNdimTime, obj.windObj.cycleAvgFiltForce - obj.noWindObj.cycleAvgSyncFiltForce)
#ax3.plot(obj.aeroObj.cycleNdimTime, obj.aeroObj.cycleAvgFiltForce/obj.ctNdimTerm)

print("delta arg max", np.argmax(obj.windObj.cycleAvgFiltForce) - np.argmax(obj.noWindObj.cycleAvgSyncFiltForce))
print("delta arg min", np.argmin(obj.windObj.cycleAvgFiltForce) - np.argmin(obj.noWindObj.cycleAvgSyncFiltForce))