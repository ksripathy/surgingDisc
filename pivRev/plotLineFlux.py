#%%
import os
import sys
import pickle
import math

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from scipy import stats

from src.ioUtilsPickle import pickleData

discPor = "70"
caseID = "06"
phaseID = "00"
    
caseObj = pickleData(dataDir+"/pickledFilledNdimData",discPor,caseID)
phaseObj = getattr(caseObj,f"phase{phaseID}")

plotX = caseObj.gridX - phaseObj.discXc
plotY = caseObj.gridY - phaseObj.discYc

xLocArr = np.arange(0.700,1.710,0.010)
yLoc1Arr = -0.8 * np.ones(len(xLocArr))
yLoc2Arr = -0.25 * np.ones(len(xLocArr))
yLoc3Arr = 0.25 * np.ones(len(xLocArr))
yLoc4Arr = 0.8 * np.ones(len(xLocArr))

omegaFlux12Arr = np.zeros(len(xLocArr))
omegaFlux34Arr = np.zeros(len(xLocArr))

fig, ax = plt.subplots()
vortPlot = ax.contourf(plotX, plotY, phaseObj.gridWz, np.linspace(-10,10,101), extend="both", cmap="seismic")
ax.set_aspect("equal")
ax.set_xlabel("x/D [-]")
ax.set_ylabel("y/D [-]")
vortCbar = fig.colorbar(vortPlot, orientation="horizontal", aspect=40)
vortCbar.set_label(r"$\omega_z D/V_\infty [-]$")

for i in range(len(xLocArr)):
    
    pt1 = [xLocArr[i],yLoc1Arr[i]]
    pt2 = [xLocArr[i],yLoc2Arr[i]]
    pt3 = [xLocArr[i],yLoc3Arr[i]]
    pt4 = [xLocArr[i],yLoc4Arr[i]]
    
    line12Data = caseObj.getWzFlux1D(pt1, pt2)
    line34Data = caseObj.getWzFlux1D(pt3, pt4)
    
    line12Pts = line12Data[0][0]
    omegaFlux12Arr[i] = line12Data[1][0][int(phaseID)] 
    
    line34Pts = line34Data[0][0]
    omegaFlux34Arr[i] = line34Data[1][0][int(phaseID)] 
    
    if i % 10 == 0:
        
        line12 = np.row_stack((line12Pts[0],line12Pts[1]))
        ax.plot(line12[:,0] - phaseObj.discXc,line12[:,1] - phaseObj.discYc,color="k")
        
        line34 = np.row_stack((line34Pts[0],line34Pts[1]))
        ax.plot(line34[:,0] - phaseObj.discXc,line34[:,1] - phaseObj.discYc,color="k")
        
CT12Arr = omegaFlux12Arr * 2
CT34Arr = omegaFlux34Arr * 2

print("CT Lower-Side Mean:", np.mean(CT12Arr))
print("CT Lower-Side std:", np.std(CT12Arr))

print("CT Upper-Side Mean:", np.mean(CT34Arr))
print("CT Upper-Side std:", np.std(CT34Arr))

loadCell45 = [0,0,0,0,0,0,0.85,0.84]
loadCell70 = [0,0,0,0,0,0,0.46,0.47]

fig2, ax2 = plt.subplots()
ax2.plot(xLocArr, CT12Arr,label="lower-half")
if discPor == "45":
    ax2.axhline(loadCell45[int(caseID)],color="k",linestyle="--",label="load-cell")
elif discPor == "70":
    ax2.axhline(loadCell70[int(caseID)],color="k",linestyle="--",label="load-cell")
ax2.legend()
ax2.set_xlabel("xLoc/D")
ax2.set_ylabel("CT")
    
fig3, ax3 = plt.subplots()
ax3.plot(xLocArr, np.abs(CT34Arr), color="tab:orange",label="upper-half")
if discPor == "45":
    ax3.axhline(loadCell45[int(caseID)],color="k",linestyle="--",label="load-cell")
elif discPor == "70":
    ax3.axhline(loadCell70[int(caseID)],color="k",linestyle="--",label="load-cell")
ax3.legend()
ax3.set_xlabel("xLoc/D")
ax3.set_ylabel("CT")
    
# %%
fig3, ax3 = plt.subplots()
velPlot = ax3.contourf(plotX, plotY, phaseObj.gridVx, cmap="seismic", levels=np.linspace(0.6,1,9),extend="both")
fig3.colorbar(velPlot, orientation="horizontal", aspect=40)


# %%
xArr = np.linspace(-1,1,51)
VxArr = caseObj.getDiscVelExtpV2()[0][2](xArr)
plt.plot(xArr, VxArr)

# %%
