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
from scipy.interpolate import CubicSpline
from scipy import signal as sig
from scipy.ndimage import median_filter
from matplotlib import colors

from src.ioUtilsPickle import pickleData

from src.miscTools import genPolygonLineBndry
from src.miscTools import genVector
from src.miscTools import nrmlVector2D
from src.miscTools import lowPassFilt

#Class for plotting axis ticks as multiple of pi
#https://stackoverflow.com/questions/40642061/how-to-set-axis-ticks-in-multiples-of-pi-python-matplotlib
class Multiple:
    def __init__(self, denominator=2, number=np.pi, latex='\pi'):
        self.denominator = denominator
        self.number = number
        self.latex = latex
        
    def locator(self):
        return plt.MultipleLocator(self.number / self.denominator)
    
    def formatter(self):
        return plt.FuncFormatter(multiple_formatter(self.denominator, self.number, self.latex))
    
def multiple_formatter(denominator=2, number=np.pi, latex='\pi'):
    
    def gcd(a, b):
        while b:
            a, b = b, a%b
        return a
    
    def _multiple_formatter(x, pos):
        den = denominator
        num = int(np.rint(den*x/number))
        com = gcd(num,den)
        (num,den) = (int(num/com),int(den/com))
        if den==1:
            if num==0:
                return r'$0$'
            if num==1:
                return r'$%s$'%latex
            elif num==-1:
                return r'$-%s$'%latex
            else:
                return r'$%s%s$'%(num,latex)
        else:
            if num==1:
                return r'$\frac{%s}{%s}$'%(latex,den)
            elif num==-1:
                return r'$\frac{-%s}{%s}$'%(latex,den)
            else:
                return r'$\frac{%s%s}{%s}$'%(num,latex,den)
            
    return _multiple_formatter

discPor = "70"
staticID = "06"
staticObj = pickleData(dataDir+"/pickledFilledNdimData",discPor,staticID)
staticSegs, staticSegCirc, staticSegWzFlux = staticObj.getCircWzFlux2D(xBounds=[-0.25,1.75],yBounds=[0.2,0.8])
staticDiscVel = staticObj.getDiscVelRev()[10,0]
print("Static Circulation:", np.sum(staticSegCirc,axis=1))
print("Static Wz Flux:", np.sum(staticSegWzFlux,axis=1))

staticFig, staticAx = plt.subplots()
maskedWz = np.ma.masked_array(staticObj.phase00.gridWz,staticObj.phase00.discMask)
maskedVx = np.ma.masked_array(staticObj.phase00.gridVx,staticObj.phase00.discMask)
discMask = np.ma.masked_array(staticObj.phase00.discMask,np.logical_not(staticObj.phase00.discMask))

plotX = staticObj.gridX - staticObj.phase00.discXc
plotY = staticObj.gridY - staticObj.phase00.discYc

allPts = np.column_stack((staticObj.gridX.flatten(),staticObj.gridY.flatten()))
vortZ = staticObj.phase00.gridWzIntpr(allPts)

fig10, ax10 = plt.subplots()
ax10.contourf(plotX, plotY,vortZ.reshape(plotX.shape), levels=np.linspace(-100,100,101), extend="both")

staticAx.contourf(plotX, plotY, discMask, cmap="gray")
staticPlot = staticAx.contourf(plotX, plotY, maskedWz, levels=np.linspace(-100,100,101), extend="both")
staticAx.plot(staticSegs[0][0][:,0] - staticObj.phase00.discXc,staticSegs[0][0][:,1] - staticObj.phase00.discYc,"k")
staticAx.plot(staticSegs[0][1][:,0] - staticObj.phase00.discXc,staticSegs[0][1][:,1] - staticObj.phase00.discYc,"k")
staticAx.plot(staticSegs[0][2][:,0] - staticObj.phase00.discXc,staticSegs[0][2][:,1] - staticObj.phase00.discYc,"k")
staticAx.plot(staticSegs[0][3][:,0] - staticObj.phase00.discXc,staticSegs[0][3][:,1] - staticObj.phase00.discYc,"k")


#%%

fig, ax = plt.subplots()
fig2, ax2 = plt.subplots()
plotPhases = np.linspace(0,2*np.pi,11)

# for dynamicID in ["00","01","02","03","04","05"]:
for dynamicID in ["00"]:
    
    dynamicObj = pickleData(dataDir+"/pickledFilledNdimData",discPor,dynamicID)
    dynamicSegs, dynamicSegCirc, dynamicSegWzFlux = dynamicObj.getCircWzFlux2D(xBounds=[-0.25,1])
    dynamicCirc = np.sum(dynamicSegCirc,axis=1)
    dynamicCirc = np.append(dynamicCirc, [dynamicCirc[2],dynamicCirc[3]])
    DGamma_DT1 = (np.roll(dynamicCirc,-1) - (np.roll(dynamicCirc,1)))/dynamicObj.DT
    DGamma_DT2 = np.gradient(dynamicCirc, 0.5*dynamicObj.DT,edge_order=2)
    dynamicWzFlux = np.sum(dynamicSegWzFlux,axis=1)
    dynamicWzFlux = np.append(dynamicWzFlux, [dynamicWzFlux[2]])
    plotCirc = dynamicCirc[2:-1]
    plotDGamma_DT1 = DGamma_DT1[2:-1]
    plotDGamma_DT2 = DGamma_DT2[2:-1]
    
    ax.plot(plotPhases,plotDGamma_DT1)
    ax.plot(plotPhases,plotDGamma_DT2)
        
    
    
ax.set_xticks(plotPhases)
ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 5))
ax.xaxis.set_major_formatter(Multiple(5).formatter())
ax.set_xlabel(r"$\phi\ [rad]$")
ax.set_ylabel(r"$\Gamma^*\ [-]$")

ax.axvline(0.5*np.pi,color="g",linestyle="--")
ax.axvline(1.5*np.pi,color="r",linestyle="--")    
# ax.axhline(np.sum(staticSegCirc,axis=1),color="k",linestyle="--")
ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.15), ncol=6, fontsize="x-small")

ax2.set_xticks(plotPhases)
ax2.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 5))
ax2.xaxis.set_major_formatter(Multiple(5).formatter())
ax2.set_xlabel(r"$\phi\ [rad]$")
ax2.set_ylabel(r"$V_x^*\ [-]$")

ax2.axvline(0.5*np.pi,color="g",linestyle="--")
ax2.axvline(1.5*np.pi,color="r",linestyle="--")    
ax2.axhline(staticDiscVel,color="k",linestyle="--")
ax2.legend(loc="upper center",bbox_to_anchor=(0.5,-0.15), ncol=6, fontsize="x-small")

#%% plot field
discPor = "70"
caseID = "06"
phaseID = "00"

plotObj = pickleData(dataDir+"/pickledFilledData",discPor,caseID)
phaseObj = getattr(plotObj,f"phase{phaseID}")

discXc = phaseObj.discXc
discYc = phaseObj.discYc

plotX = plotObj.gridX - discXc
plotY = plotObj.gridY - discYc

vortFig, vortAx = plt.subplots()
maskedWz = np.ma.masked_array(phaseObj.gridWz,phaseObj.discMask)
maskedVx = np.ma.masked_array(phaseObj.gridVx,phaseObj.discMask)
discMask = np.ma.masked_array(phaseObj.discMask,np.logical_not(phaseObj.discMask))

vortAx.contourf(plotX, plotY, discMask, cmap="gray")
vortPlot = vortAx.contourf(plotX, plotY, maskedWz, levels=np.linspace(-100,100,101), extend="both")
vortAx.plot([-0.25,0.25],[-0.8,-0.8],"k")
vortAx.plot([0.25,0.25],[-0.8,0],"k")
vortAx.plot([0.25,-0.25],[0,0],"k")
vortAx.plot([-0.25,-0.25],[0,-0.8],"k")
vortAx.set_aspect("equal")

vortAx.set_xlabel(r"$x^{*}\ [-]$")
vortAx.set_ylabel(r"$y^{*}\ [-]$")
vortCbar = vortFig.colorbar(vortPlot, orientation="horizontal",aspect=40)
vortCbar.set_label(r"$\omega_z^{*}\ [-]$")

velFig, velAx = plt.subplots()
divnorm = colors.TwoSlopeNorm(vmin=np.min(maskedVx), vcenter=1.0, vmax=np.max(maskedVx))
velAx.contourf(plotX, plotY, discMask, cmap="gray")
velPlot = velAx.contourf(plotX, plotY, maskedVx, levels=np.linspace(0,np.max(maskedVx),21),norm=divnorm)
velAx.plot([-0.25,1],[-0.8,-0.8],"k")
velAx.plot([1,1],[-0.8,0],"k")
velAx.plot([1,-0.25],[0,0],"k")
velAx.plot([-0.25,-0.25],[0,-0.8],"k")

velAx.set_xlabel(r"$x^{*}\ [-]$")
velAx.set_ylabel(r"$y^{*}\ [-]$")
velCbar = vortFig.colorbar(velPlot, orientation="horizontal",aspect=40)
velCbar.set_label(r"$V_x^{*}\ [-]$")


    
# %%
