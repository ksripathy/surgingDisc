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

discPor = "45"
caseID = "07"
phaseID = "00"

with open(dataDir + f"/pickledPlotNdimDataSegs/p{discPor}Case{caseID}.pickle","rb") as handle:
    pickleObj = pickle.load(handle)
    
phaseObj = getattr(pickleObj,f"phase{phaseID}")

plotX = pickleObj.gridX - phaseObj.discXc
plotY = pickleObj.gridY - phaseObj.discYc

fig3, ax3 = plt.subplots()
plot = ax3.contourf(plotX, plotY, phaseObj.gridWz, np.linspace(-10,10,101), extend="both")
ax3.set_aspect("equal")
fig3.colorbar(plot, ax=ax3)
    
CVArr = pickleObj.lowerCVArr
CVArrRef = CVArr[90:,:].copy()

#%%Functions

def getHist(CVArr, attrb="omegaFlux", phaseID="00",binQty=100):
    
    data = pickleObj.getCntrlVolArrData(pickleObj, CVArr, attrb=attrb,plotPhaseID=phaseID)
    hist = stats.relfreq(data.flatten(),numbins=binQty)
    freq = hist.frequency
    bins = hist.lowerlimit + np.linspace(0,hist.binsize*freq.size,freq.size)
    
    return data, bins, freq

def lowPassCV(CVArr, value, attrb="omegaFlux", phaseID="00"):
    
    data = pickleObj.getCntrlVolArrData(pickleObj, CVArr, attrb=attrb,plotPhaseID=phaseID)
    indxs = np.argwhere(data < value)
    
    return CVArr[indxs[:,0],indxs[:,1]]

def highPassCV(CVArr, value, attrb="omegaFlux", phaseID="00"):
    
    data = pickleObj.getCntrlVolArrData(pickleObj, CVArr, attrb=attrb,plotPhaseID=phaseID)
    indxs = np.argwhere(data > value)
    
    return CVArr[indxs[:,0],indxs[:,1]]

def bandPassCV(CVArr, value, attrb="omegaFlux", phaseID="00", cutOffCrit=5):
    
    data = pickleObj.getCntrlVolArrData(pickleObj, CVArr, attrb=attrb,plotPhaseID=phaseID)
    binWidth = (np.max(data) - np.min(data)) / 100
    print("Band pass CV bin width:", binWidth)
    cutOff1 = value - (binWidth*cutOffCrit)
    cutOff2 = value + (binWidth*cutOffCrit)
    indxs = np.argwhere(np.logical_and(data >= cutOff1,data <= cutOff2))
    
    return CVArr[indxs[:,0],indxs[:,1]]

def plotCV(axObj,CVArr,color="black"):
    
    for CV in CVArr:
        
        xLim = CV.xLimits
        yLim = CV.yLimits
        poly_coords = [(xLim[0],yLim[0]),(xLim[1],yLim[0]),(xLim[1],yLim[1]),(xLim[0],yLim[1])]
        axObj.add_patch(Polygon(poly_coords,edgecolor=color,fill=None))
        axObj.text(0.5*(xLim[0]+xLim[1]),0.5*(yLim[0]+yLim[0])-0.09,"1")
        axObj.text(0.5*(xLim[1]+xLim[1])+0.01,0.5*(yLim[0]+yLim[1]),"2")
        axObj.text(-0.035,0.5*(yLim[1]+yLim[1])+0.01,"3")
        axObj.text(0.5*(xLim[0]+xLim[0])-0.08,0.5*(yLim[0]+yLim[1]),"4")


# %% OmegaFlux histogram
if int(caseID) > 5:
    for CV in CVArr.flatten():    
        CV.curlForce = - CV.omegaFluxSegs[:,1]
        
    for CV in CVArrRef.flatten():    
        CV.curlForce = - CV.omegaFluxSegs[:,1]
    
    data, bins, freq = getHist(CVArr,attrb="curlForce",phaseID=phaseID)
    dataRef, binsRef, freqRef = getHist(CVArrRef,attrb="curlForce",phaseID=phaseID, binQty=20)
    
else:
    
    for CV in CVArr.flatten():    
        CV.curlForce = np.sum(CV.DGamma_DTSegs,axis=1) - CV.omegaFluxSegs[:,1]
        
    for CV in CVArrRef.flatten():    
        CV.curlForce = np.sum(CV.DGamma_DTSegs,axis=1) - CV.omegaFluxSegs[:,1]
        
    data, bins, freq = getHist(CVArr,attrb="curlForce",phaseID=phaseID)
    dataRef, binsRef, freqRef = getHist(CVArrRef,attrb="curlForce",phaseID=phaseID, binQty=20)

fig1, ax1 = plt.subplots()
ax1.plot(bins, freq, label = "all CV")
ax1.plot(binsRef, freqRef, label = "CV 1.2 to 1.5")
ax1.legend()
ax1.set_xlabel(r"$\int_V \nabla \times \frac{f}{\rho} dV$")
ax1.set_ylabel("Relative frequency")
ax1.grid()

curlForce = bins[np.argmax(freq)]*phaseObj.VInf**2 #Dimensionalizing
forceDensity = curlForce * phaseObj.rhoInf / 3e-3
force = forceDensity * np.pi * 0.1**2 * 3e-3
CT = force / (0.5 * phaseObj.rhoInf * phaseObj.VInf**2 * np.pi * 0.1**2)

print("CT :", CT)


curlForceRef = binsRef[np.argmax(freqRef)]*phaseObj.VInf**2 #Dimensionalizing
forceDensityRef = curlForceRef * phaseObj.rhoInf / 3e-3
forceRef = forceDensityRef * np.pi * 0.1**2 * 3e-3
CTRef = forceRef / (0.5 * phaseObj.rhoInf * phaseObj.VInf**2 * np.pi * 0.1**2)

print("CT Ref:", CTRef)

#%%Control volume plots
CVArr1 = lowPassCV(CVArr, value=0.15, attrb="curlForce", phaseID=phaseID)
CVArr2 = highPassCV(CVArr, value=0.30, attrb="curlForce", phaseID=phaseID)
CVArr3 = bandPassCV(CVArr, value=-0.48, attrb="curlForce", phaseID=phaseID, cutOffCrit=1)

fig2, ax2 = plt.subplots()
vortPlot = ax2.contourf(plotX[:,75:-75], plotY[:,75:-75], phaseObj.gridWz[:,75:-75], np.linspace(-10,10,101), extend="both", cmap="seismic")
ax2.set_aspect("equal")
ax2.set_xlabel("x/D [-]")
ax2.set_ylabel("y/D [-]")
vortCbar = fig2.colorbar(vortPlot, orientation="horizontal", aspect=40)
vortCbar.set_label(r"$\omega_z D/V_\infty [-]$")
# plotCV(ax2,CVArr1,"blue")
# plotCV(ax2,CVArr2,"red")
plotCV(ax2,np.array([CVArr[40,20]]),"black")
# fig2.savefig(plotDir + "/windTechField.png",dpi=600,pad_inches=0,bbox_inches="tight")

# %%Zero vorticity flux analysis
caseObj = pickleData(dataDir+"/pickledFilledNdimData",discPor,caseID)
lowerCV2, upperCV2 = caseObj.circWzParamStudy(xMin=0.3,xMax=0.6,deltaX=0.01,deltaY=0.01)
lowerCV3, upperCV3 = caseObj.circWzParamStudy(xMin=0.6,xMax=0.9,deltaX=0.01,deltaY=0.01)
lowerCV4, upperCV4 = caseObj.circWzParamStudy(xMin=0.9,xMax=1.2,deltaX=0.01,deltaY=0.01)
lowerCV5, upperCV5 = caseObj.circWzParamStudy(xMin=1.2,xMax=1.5,deltaX=0.01,deltaY=0.01)


# %%
CV2Arr = upperCV2
CV3Arr = upperCV3
CV4Arr = upperCV4
CV5Arr = upperCV5

for CV2 in CV2Arr.flatten():
    CV2.curlForce = -np.sum(CV2.omegaFluxSegs,axis=1)
for CV3 in CV3Arr.flatten():
    CV3.curlForce = -np.sum(CV3.omegaFluxSegs,axis=1)
for CV4 in CV4Arr.flatten():
    CV4.curlForce = -np.sum(CV4.omegaFluxSegs,axis=1)
for CV5 in CV5Arr.flatten():
    CV5.curlForce = -np.sum(CV5.omegaFluxSegs,axis=1)
    
bins2, freq2 = getHist(CV2Arr,attrb="curlForce",binQty=10)
bins3, freq3 = getHist(CV3Arr,attrb="curlForce",binQty=10)
bins4, freq4 = getHist(CV4Arr,attrb="curlForce",binQty=10)
bins5, freq5 = getHist(CV5Arr,attrb="curlForce",binQty=10)

fig4, ax4 = plt.subplots()
ax4.plot(bins2,freq2,"b",label="CV1L")
ax4.plot(bins3,freq3,"g",label="CV2L")
ax4.plot(bins4,freq4,"orange",label="CV3L")
ax4.plot(bins5,freq5,"r",label="CV4L")
ax4.grid()
ax4.legend()

fig6,ax6 = plt.subplots()
ax6.contourf(plotX, plotY, phaseObj.gridWz, np.linspace(-10,10,101), extend="both")
ax6.set_aspect("equal")
plotCV(ax6,np.array([CV2Arr[0,-1]]),"b")
plotCV(ax6,np.array([CV3Arr[0,-1]]),"g")
plotCV(ax6,np.array([CV4Arr[0,-1]]),"orange")
plotCV(ax6,np.array([CV5Arr[0,-1]]),"r")


# %%
CV2Arr3 = bandPassCV(CV2Arr, value=0, attrb="curlForce", phaseID=phaseID, cutOffCrit=10)

fig5, ax5 = plt.subplots()
ax5.contourf(plotX, plotY, phaseObj.gridWz, np.linspace(-10,10,101), extend="both")
ax5.set_aspect("equal")
# plotCV(ax2,CVArr1,"blue")
# plotCV(ax2,CVArr2,"red")
plotCV(ax5,CV2Arr3,"black")
# %%
