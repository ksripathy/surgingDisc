#%%
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

from src.ioUtilsPickle import pickleData
from src.miscTools import genPolygonLineBndry

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

# pickleObj = pickleData(dataDir+"/pickledRawDataIntpr","45","00")
# staticPickleObj = pickleData(dataDir+"/pickledRawDataIntpr","45","06")

staticPickleObj = pickleData(dataDir+"/pickledFilledData","45","06")
boundLines = staticPickleObj.getCirc2D(xBounds=[-0.25,1],yBounds=[-0.8,0])[0]
# boundLines = staticPickleObj.getCirc1D(pt1=[0,-0.45],pt2=[0,-0.55])[0]
boundLines = np.append(boundLines,boundLines[0,:].reshape(1,2),axis=0)
print("Circulation:",np.sum(staticPickleObj.circ))

fig5, ax5 = plt.subplots()
ax5.contourf(staticPickleObj.gridX - staticPickleObj.phase00.discXc,staticPickleObj.gridY - staticPickleObj.phase00.discYc,getattr(staticPickleObj,f"phase00").discMask,cmap="gray")
vortPlot = ax5.contourf(staticPickleObj.gridX - staticPickleObj.phase00.discXc,staticPickleObj.gridY - staticPickleObj.phase00.discYc,np.ma.masked_array(getattr(staticPickleObj,f"phase00").gridWz,getattr(staticPickleObj,f"phase00").discMask),levels=np.linspace(-100,100,101),extend="both")
ax5.plot(boundLines[:,0] - staticPickleObj.phase00.discXc,boundLines[:,1] - staticPickleObj.phase00.discYc, "k")
ax5.set_aspect("equal")
ax5.set_xlabel(r"$x^{*}\ [-]$")
ax5.set_ylabel(r"$y^{*}\ [-]$")
vortCbar = fig5.colorbar(vortPlot, orientation="horizontal",aspect=40)
vortCbar.set_label(r"$\omega_z^{*}\ [-]$")

#%%
circArr = np.zeros((6,12))
plotPhases = np.linspace(0,2*np.pi,11)

fig4, ax4 = plt.subplots()

ptList = np.empty(6,dtype=np.ndarray)

for i, caseID in enumerate(["00","01","02","03","04","05"]):
# for i, caseID in enumerate(["05"]):
    pickleObj = pickleData(dataDir+"/pickledFilledData","70",caseID)
    staticPickleObj = pickleData(dataDir+"/pickledFilledData","70","06")

    staticPickleObj.getCirc2D(xBounds=[-0.25,1],yBounds=[-0.8,0])

    ptList[i] = pickleObj.getCirc2D(xBounds=[-0.25,1],yBounds=[-0.8,0])
    circ = np.sum(pickleObj.circ,axis=1)
    circArr[i,:] = circ
    
    plotCirc = circ[2:]
    plotCirc = np.append(plotCirc,plotCirc[0])
    
    ax4.plot(plotPhases,plotCirc,label=f"case{caseID}",marker="x")
    
ax4.set_xticks(plotPhases)
ax4.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 5))
ax4.xaxis.set_major_formatter(Multiple(5).formatter())
ax4.set_xlabel(r"$\phi\ [rad]$")
ax4.set_ylabel(r"$\Gamma^*\ [-]$")

ax4.axvline(0.5*np.pi,color="g",linestyle="--")
ax4.axvline(1.5*np.pi,color="r",linestyle="--")    
ax4.axhline(np.sum(staticPickleObj.circ),color="k",linestyle="--")
ax4.legend(loc="upper center",bbox_to_anchor=(0.5,-0.15), ncol=6, fontsize="x-small")

#%%
fig3, ax3 = plt.subplots()
ax3.plot(plotPhases,plotCirc)
ax3.axhline(np.sum(staticPickleObj.circ),color="k",linestyle="--")
ax3.axvline(0.25,color="g",linestyle="--")
ax3.axvline(0.75,color="r",linestyle="--")

caseID = "04"
phaseID = "05"

boundLines = ptList[int(caseID)][int(phaseID)]
boundLines = np.append(boundLines,boundLines[0,:].reshape(1,2),axis=0)

pickleObj = pickleData(dataDir+"/pickledFilledData","70",caseID)
discXc = getattr(pickleObj,f"phase{phaseID}").discXc
discYc = getattr(pickleObj,f"phase{phaseID}").discYc

fig2, ax2 = plt.subplots()
# ax2.contourf(pickleObj.gridX,pickleObj.gridY,np.ma.masked_array(getattr(pickleObj,f"phase{phaseID}").gridWz,getattr(pickleObj,f"phase{phaseID}").gridMask),levels=np.linspace(-100,100,101),extend="both")
ax2.contourf(pickleObj.gridX - discXc,pickleObj.gridY - discYc,getattr(pickleObj,f"phase{phaseID}").discMask,cmap="gray")
# ax2.contourf(pickleObj.gridX,pickleObj.gridY,np.ma.masked_array(getattr(pickleObj,f"phase{phaseID}").gridVx,getattr(pickleObj,f"phase{phaseID}").discMask))
vortPlot2 = ax2.contourf(pickleObj.gridX - discXc,pickleObj.gridY - discYc,np.ma.masked_array(getattr(pickleObj,f"phase{phaseID}").gridWz,getattr(pickleObj,f"phase{phaseID}").discMask),levels=np.linspace(-100,100,101),extend="both")
ax2.plot(boundLines[:,0] - discXc,boundLines[:,1] - discYc, "k")
ax2.set_aspect("equal")

ax2.set_xlabel(r"$x^{*}\ [-]$")
ax2.set_ylabel(r"$y^{*}\ [-]$")
vortCbar2 = fig2.colorbar(vortPlot2, orientation="horizontal",aspect=40)
vortCbar2.set_label(r"$\omega_z^{*}\ [-]$")

 #%%

fig6, ax6 = plt.subplots()
velTang = pickleObj.phase00.gridVxIntpr(boundLines)


#%%

points = np.array([[-1,-0.85],
          [1,-0.85],
          [1,0.85],
          [-1,0.85]])

bndryLineSegs, cntrlPtSegs = genPolygonLineBndry(points,pickleObj.gridX[0,1] - pickleObj.gridX[0,0])

phaseID = "00"

circSegs, circPts = pickleObj.getCirc2D(phaseID, bndryLineSegs, cntrlPtSegs)

fig1, ax1 = plt.subplots()
ax1.contourf(pickleObj.gridX,pickleObj.gridY,np.ma.masked_array(getattr(pickleObj,f"phase{phaseID}").gridWz,getattr(pickleObj,f"phase{phaseID}").gridMask),levels=np.linspace(-100,100,101),extend="both")
ax1.plot(bndryLineSegs[0][:,0],bndryLineSegs[0][:,1], linewidth=3)
ax1.plot(bndryLineSegs[1][:,0],bndryLineSegs[1][:,1], linewidth=3)
ax1.plot(bndryLineSegs[2][:,0],bndryLineSegs[2][:,1], linewidth=3)
ax1.plot(bndryLineSegs[3][:,0],bndryLineSegs[3][:,1], linewidth=3)
ax1.set_aspect("equal")


#%%
# pickleObj.genInterpolator()
pickleObj.getDiscVelRaw()
temp = pickleObj.VxuRaw
polyFitObj = pickleObj.getDiscVelExtp()

#%%
xArr = pickleObj.xArrVxcExtp[:,2].flatten()
y = pickleObj.phase03.discYc

line = np.column_stack((xArr,np.tile(y,len(xArr))))
vxArr = pickleObj.phase02.gridVxIntpr(line)
polyFit3 = np.polyfit(xArr, vxArr, deg=3)
vxFit = np.poly1d(polyFit3)(xArr)

fig, ax = plt.subplots()
ax.plot(xArr, vxArr)
ax.plot(xArr, polyFitObj[10,2](xArr))
ax.plot(xArr, vxFit)

# pickleContainer = pickleObj.packData()
# temp = pickleContainer.discVxcExtp
# plt.contourf(pickleContainer.gridX,pickleContainer.gridY,pickleContainer.phase00.gridWz)



# %%
