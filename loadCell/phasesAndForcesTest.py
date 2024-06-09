import os
import sys

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

import numpy as np
import glob
import matplotlib.pyplot as plt

from src.miscTools import freqsSignal
from src.phaseAndForces import recordedData


# discPor = "70" #in percentage
# kword = "Case[0-9][0-9].txt"
# loadCell = "2N"

def filePaths(discPor, kword, loadCell):
    
    totalLoadsPath = dataDir + f"/testMatrixTotal/curved/{loadCell}/p{discPor}Load{kword}.txt"
    totalLogPath = dataDir + f"/testMatrixTotal/curved/{loadCell}/p{discPor}Log{kword}.txt"

    inertialLoadsPath = dataDir + f"/testMatrixInertial/{loadCell}/p{discPor}Load{kword}.txt"
    inertialLogPath = dataDir + f"/testMatrixInertial/{loadCell}/p{discPor}Log{kword}.txt"
    
    return totalLoadsPath,totalLogPath,inertialLoadsPath,inertialLogPath

discPor = "45"
loadCell="2N"
cases=["Case01"]

fig1, ax1 = plt.subplots()
# fig3, ax3 = plt.subplots()

objs = []

for case in cases:
    
    print(case)
    
    obj = recordedData(*filePaths(discPor,case,loadCell))
    objs.append(obj)
    
    cycleNo = 0
    plotSlice = np.arange(cycleNo*obj.surgeCycleSize,(cycleNo+200)*obj.surgeCycleSize).astype(int)
    
    #ax1.plot(obj.inertialLoadPhases[plotSlice],obj.totalCt[plotSlice],label=case)
    ax1.plot(obj.inertialLoadPhases[plotSlice],obj.filtAeroCt[plotSlice],label=case)
    # ax3.plot(obj.userDefinedCyclePhases, obj.cycleAvgFiltAeroCt,label=case)
    
if discPor == "45":
    
    staticCt = 0.87
    
elif discPor == "70":
    
    staticCt = 0.53
    
ax1.axhline(staticCt,linestyle="--",color="k",label="static")
ax1.legend()
ax1.grid()
ax1.set_xlabel(r"$\phi [\degree]$")
ax1.set_ylabel(r"$C_T [-]$")

# ax3.axhline(staticCt,linestyle="--",color="k",label="static")
# ax3.legend()
# ax3.grid()
# ax3.set_xlabel(r"$\phi [\degree]$")
# ax3.set_ylabel(r"$C_T [-]$")


objAttrb = getattr(objs[0],"totalCt")
freq, amp = freqsSignal(objAttrb - np.mean(objAttrb),2e3)
fig2, ax2 = plt.subplots()
ax2.stem(freq[:5000],amp[:5000])
    



# obj = recordedData(*filePaths(70,"Case00","20N"))

# totalCycleSizes = []
# inertialCycleSizes = []

# for i in range(len(obj.totalCycleCt)):
    
#     totalCycleSizes.append(len(obj.totalCycleLoadPhases[i]))
#     inertialCycleSizes.append(len(obj.inertialCycleLoadPhases[i]))
    
# deltaPhases = obj.totalLoadPhases - obj.inertialLoadPhases

# for i,deltaPhase in enumerate(deltaPhases):
    
#     if deltaPhase >= 350:
        
#         deltaPhases[i] = 360 - deltaPhase
        
#     elif deltaPhase <= -350:
        
#         deltaPhases[i] = 360 + deltaPhase
         
# rmsDelta = np.sqrt(np.mean(deltaPhases**2))
    
# print("Delta total cycles zize:", max(totalCycleSizes) - min(totalCycleSizes))
# print("Delta inertial cycles zize:", max(inertialCycleSizes) - min(inertialCycleSizes))
# print(totalCycleSizes == inertialCycleSizes)
# print("length total loads size:", len(obj.totalCt))
# print("length inertial loads size:", len(obj.inertialCt))
# print("Sum total cycles size:", sum(totalCycleSizes))
# print("Sum inertial cycles size:", sum(inertialCycleSizes))
# print("Surge cycle size:", obj.surgeCycleSize)
# print("rmsDelta of phases:", rmsDelta)
# print("Max deltaPhases", max(deltaPhases))
# print("Min deltaPhases", min(deltaPhases))

# cycleNo = 100

# plotSlice = np.arange(cycleNo*obj.surgeCycleSize,(cycleNo+1)*obj.surgeCycleSize).astype(int)

# fig1, ax1 = plt.subplots()
# ax1.plot(obj.totalCt[plotSlice])
# ax1.plot(obj.inertialCt[plotSlice])
# ax1.plot(obj.aeroCt[plotSlice])

# fig2, ax2 = plt.subplots()
# ax2.plot(obj.cycleAvgTotal)
# ax2.plot(obj.cycleAvgInertial)
# ax2.plot(obj.cycleAvgAero)

# fig3, ax3 = plt.subplots()
# ax3.plot(obj.filtTotalCt[plotSlice])
# ax3.plot(obj.filtInertialCt[plotSlice])
# ax3.plot(obj.filtAeroCt[plotSlice])
