import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
from scipy.io import loadmat
import sys

sys.path.append("..")

from src.miscTools import lagShift
from src.miscTools import symAvgArray
from piv.src.ioUtilsNoStitch import loadData

class dynInflow:
    
    def __init__(self, caseNo):
        
        caseDict = np.array([[1.0, 0.250, 2.387324, 5.0e-2, 1.196],
                             [2.0, 0.125, 4.774648, 2.5e-2, 1.199], 
                             [2.0, 0.250, 4.774648, 5.0e-2, 1.198], 
                             [2.0, 0.375, 3.183099, 7.5e-2, 1.192], 
                             [3.0, 0.250, 4.774648, 5.0e-2, 1.194],
                             [2.7, 0.375, 4.297183, 7.5e-2, 1.189],
                             [0.0, 0.000, 0.000000, 0.0e-2, 1.195],
                             [0.0, 0.000, 0.000000, 0.0e-2, 1.195]])
        
        discArea = 0.25 * np.pi * (20e-2)**2
        p70FstreamRho = np.array([1.196899, 1.198764, 1.198306, 1.194673, 1.192900, 1.189482, 1.194618, 1.195081])
        
        staticCtLoadCell = [0.60, 0.60, 0.60, 0.56]
        staticCtPiv = [0.56, 0.56, 0.56, 0.56, 0.56, 0.56]
        self.dynCtPiv = np.abs(loadmat("/home/ksripathy/phd/surgingDisc/piv/data/loads/"+f"p70Case{caseNo}")["dynCT"][0])
        self.dynCt = np.zeros(len(self.dynCtPiv)) #Constructed using piv velocity data
        
        if caseNo < 3:
                
            staticCase = loadData("/home/ksripathy/phd/surgingDisc/piv/plots/results_p70Case6.json")
            rhoStatic = caseDict[6,4]
        
        else:
        
            staticCase = loadData("/home/ksripathy/phd/surgingDisc/piv/plots/results_p70Case7.json")
            rhoStatic = caseDict[7,4]
            
        dynCase = loadData(f"/home/ksripathy/phd/surgingDisc/piv/plots/results_p70Case{caseNo}.json")
        
        vinfStatic = np.abs(staticCase.phaseVinf)
        qinfStatic = 0.5 * rhoStatic * vinfStatic**2
        self.discVelStatic = 1 - np.array(staticCase.cycleAxialInd["intp200"][0])#Ndim
        
        vinfDyn = np.abs(np.array(dynCase.phaseVinf))
        rhoDyn = caseDict[caseNo,4]
        qinfDyn = 0.5 * rhoDyn * vinfDyn**2
        qinfDynTheory = 0.5 * rhoDyn * np.mean(vinfDyn)**2

        self.phasesTheory = 2*np.pi * np.linspace(0.05,1.15,241)
        freq = caseDict[caseNo,2]
        amp = caseDict[caseNo,3]
        omega = 2 * np.pi * freq
        self.surgeVelTheory = amp * omega * np.cos(self.phasesTheory)/np.mean(vinfDyn)
        cenLineInd = staticCase.cycleAxialInd["intp200"][0][40] #Centreline disc induction
        self.discVelDynTheory = (1 - cenLineInd) + (cenLineInd * self.surgeVelTheory)#ndim
        self.discVelDynTheoryOld = (1 - cenLineInd) + self.surgeVelTheory#ndim
        self.dynCtTheory = staticCtPiv[caseNo] + 4 * cenLineInd * self.surgeVelTheory * ((1 + cenLineInd * self.surgeVelTheory) - 2 * self.discVelDynTheory)
        
        self.phases = 2*np.pi * np.linspace(0.05, 1.15, 12)
        self.surgeVel = np.zeros(len(self.phases))
        self.discVelDyn = np.empty(shape=len(self.discVelStatic), dtype=np.ndarray) #Ndim
        
        for radius in range(len(self.discVelStatic)):
            
            self.discVelDyn[radius] = np.zeros(12)
            localDiscInd = staticCase.cycleAxialInd["intp200"][0][radius]
            
            for phase in range(12):
                
                self.discVelDyn[radius][phase] = 1 - dynCase.cycleAxialInd["intp200"][phase][radius]
                
                if radius == 40:
                    self.surgeVel[phase] = amp * omega * np.cos(self.phases[phase])/vinfDyn[phase]
                    self.dynCt[phase] = staticCtPiv[caseNo] + 4 * cenLineInd * self.surgeVel[phase] * ((1 + cenLineInd * self.surgeVel[phase]) - 2 * self.discVelDyn[radius][phase])
                

            
caseNo = 0

p70Case0 = dynInflow(0)
p70Case1 = dynInflow(1)
p70Case2 = dynInflow(2)
p70Case3 = dynInflow(3)
p70Case4 = dynInflow(4)
p70Case5 = dynInflow(5)

defaultLocs = [0,1,3,4,5,6,8,9,10,11]
phases = [r"$18\degree$",r"$54\degree$",r"$90\degree$",r"$126\degree$",
          r"$162\degree$",r"$198\degree$",r"$234\degree$",r"$270\degree$",
          r"$306\degree$",r"$342\degree$"]

fig,ax = plt.subplots(ncols=2)
fig.set_size_inches(26/2.54, 10/2.54)
ax[0].plot(p70Case0.surgeVelTheory,1 - p70Case0.discVelDynTheory, label="case0")
ax[0].plot(p70Case1.surgeVelTheory,1 - p70Case0.discVelDynTheory, label="case1")
ax[0].plot(p70Case0.surgeVel,1 - p70Case0.discVelDyn[40], marker = "o",linestyle=":",color=ax[0].lines[0].get_color())
ax[0].plot(p70Case1.surgeVel,1 - p70Case1.discVelDyn[40], marker = "o",linestyle=":",color=ax[0].lines[1].get_color())
ax[0].legend()
ax[0].set_xlabel(r"$V_{sur}/V_{\infty}\ [-]$")
ax[0].set_ylabel(r"$v_D/V_{\infty}\ [-]$")
for x,y,z in zip(p70Case0.surgeVel[:10],1 - p70Case0.discVelDyn[40][:10],range(10)):
    ax[0].text(x,y,phases[z],color="k",fontsize=12)
ax[1].plot(p70Case3.surgeVelTheory,1 - p70Case3.discVelDynTheory, label="case3")
ax[1].plot(p70Case4.surgeVelTheory,1 - p70Case4.discVelDynTheory, label="case4")
ax[1].plot(p70Case3.surgeVel,1 - p70Case3.discVelDyn[40], marker = "o",linestyle=":",color=ax[1].lines[0].get_color())
ax[1].plot(p70Case4.surgeVel,1 - p70Case4.discVelDyn[40], marker = "o",linestyle=":",color=ax[1].lines[1].get_color())
ax[1].legend()
ax[1].set_xlabel(r"$V_{sur}/V_{\infty}\ [-]$")
for x,y,z in zip(p70Case3.surgeVel[:9],1 - p70Case3.discVelDyn[40][:9],range(9)):
    ax[1].text(x,y,phases[z],color="k",fontsize=12)
fig.savefig("dynInflowTorque24.png", dpi=300, bbox_inches="tight", pad_inches=0)

fig0,ax0 = plt.subplots(ncols=2)
fig0.set_size_inches(26/2.54, 10/2.54)
ax0[0].plot(p70Case1.surgeVelTheory,1 - p70Case1.discVelDynTheory, label="case1")
ax0[0].plot(p70Case2.surgeVelTheory,1 - p70Case2.discVelDynTheory, label="case2")
ax0[0].plot(p70Case3.surgeVelTheory,1 - p70Case3.discVelDynTheory, label="case3")
ax0[0].plot(p70Case1.surgeVel,1 - p70Case1.discVelDyn[40], marker = "x",markevery=defaultLocs,linestyle=":",color=ax0[0].lines[0].get_color())
ax0[0].plot(p70Case2.surgeVel,1 - p70Case2.discVelDyn[40], marker = "x",markevery=defaultLocs,linestyle=":",color=ax0[0].lines[1].get_color())
ax0[0].plot(p70Case3.surgeVel,1 - p70Case3.discVelDyn[40], marker = "x",markevery=defaultLocs,linestyle=":",color=ax0[0].lines[2].get_color())
ax0[0].plot(p70Case1.surgeVel[2],1 - p70Case1.discVelDyn[40][2],marker="o",markerfacecolor="g", color="k", linestyle="None", label = r"$90\degree phase$")
ax0[0].plot(p70Case1.surgeVel[7],1 - p70Case1.discVelDyn[40][7],marker="o",markerfacecolor="r", color="k", linestyle="None", label = r"$270\degree phase$")
ax0[0].plot(p70Case2.surgeVel[2],1 - p70Case2.discVelDyn[40][2],marker="o",markerfacecolor="g", color="k", linestyle="None")
ax0[0].plot(p70Case2.surgeVel[7],1 - p70Case2.discVelDyn[40][7],marker="o",markerfacecolor="r", color="k", linestyle="None")
ax0[0].plot(p70Case3.surgeVel[2],1 - p70Case3.discVelDyn[40][2],marker="o",markerfacecolor="g", color="k", linestyle="None")
ax0[0].plot(p70Case3.surgeVel[7],1 - p70Case3.discVelDyn[40][7],marker="o",markerfacecolor="r", color="k", linestyle="None")

fig1, ax1 = plt.subplots(nrows=3,ncols=2)
fig1.set_size_inches(26/2.54, 30/2.54)
ax1[0,0].set_title("Case0")
ax1[0,0].plot(p70Case0.surgeVelTheory,p70Case0.discVelDynTheory, label = "quasi-steady")
#ax1[0,0].plot(p70Case0.surgeVelTheory,p70Case0.discVelDynTheoryOld, label = "quasi-steady-old")
ax1[0,0].plot(p70Case0.surgeVel,p70Case0.discVelDyn[40],marker="x",markevery=defaultLocs,linestyle=":",color="k", label = "piv")
ax1[0,0].plot(p70Case0.surgeVel[2],p70Case0.discVelDyn[40][2],marker="o",markerfacecolor="g", color="k", linestyle="None", label = r"$90\degree phase$")
ax1[0,0].plot(p70Case0.surgeVel[7],p70Case0.discVelDyn[40][7],marker="o",markerfacecolor="r", linestyle="None", color="k", label = r"$270\degree phase$")
ax1[0,0].set_ylabel(r"$V_D/V_{\infty}\ [-]$")
ax1[0,0].legend()
ax1[0,1].set_title("Case1")
ax1[0,1].plot(p70Case1.surgeVelTheory,p70Case1.discVelDynTheory)
#ax1[0,1].plot(p70Case1.surgeVelTheory,p70Case1.discVelDynTheoryOld)
ax1[0,1].plot(p70Case1.surgeVel,p70Case1.discVelDyn[40],marker="x",markevery=defaultLocs,linestyle=":",color="k")
ax1[0,1].plot(p70Case1.surgeVel[2],p70Case1.discVelDyn[40][2],marker="o",markerfacecolor="g", color="k")
ax1[0,1].plot(p70Case1.surgeVel[7],p70Case1.discVelDyn[40][7],marker="o",markerfacecolor="r", color="k")
ax1[1,0].set_title("Case2")
ax1[1,0].plot(p70Case2.surgeVelTheory,p70Case2.discVelDynTheory)
#ax1[1,0].plot(p70Case1.surgeVelTheory,p70Case2.discVelDynTheoryOld)
ax1[1,0].plot(p70Case2.surgeVel,p70Case2.discVelDyn[40],marker="x",markevery=defaultLocs,linestyle=":",color="k")
ax1[1,0].plot(p70Case2.surgeVel[2],p70Case2.discVelDyn[40][2],marker="o",markerfacecolor="g", color="k")
ax1[1,0].plot(p70Case2.surgeVel[7],p70Case2.discVelDyn[40][7],marker="o",markerfacecolor="r", color="k")
ax1[1,0].set_ylabel(r"$V_D/V_{\infty}\ [-]$")
ax1[1,1].set_title("Case3")
ax1[1,1].plot(p70Case3.surgeVelTheory,p70Case3.discVelDynTheory)
#ax1[1,1].plot(p70Case3.surgeVelTheory,p70Case3.discVelDynTheoryOld)
ax1[1,1].plot(p70Case3.surgeVel,p70Case3.discVelDyn[40],marker="x",markevery=defaultLocs,linestyle=":",color="k")
ax1[1,1].plot(p70Case3.surgeVel[2],p70Case3.discVelDyn[40][2],marker="o",markerfacecolor="g", color="k")
ax1[1,1].plot(p70Case3.surgeVel[7],p70Case3.discVelDyn[40][7],marker="o",markerfacecolor="r", color="k")
ax1[2,0].set_title("Case4")
ax1[2,0].plot(p70Case4.surgeVelTheory,p70Case4.discVelDynTheory)
#ax1[2,0].plot(p70Case3.surgeVelTheory,p70Case4.discVelDynTheoryOld)
ax1[2,0].plot(p70Case4.surgeVel,p70Case4.discVelDyn[40],marker="x",markevery=defaultLocs,linestyle=":",color="k")
ax1[2,0].plot(p70Case4.surgeVel[2],p70Case4.discVelDyn[40][2],marker="o",markerfacecolor="g", color="k")
ax1[2,0].plot(p70Case4.surgeVel[7],p70Case4.discVelDyn[40][7],marker="o",markerfacecolor="r", color="k")
ax1[2,0].set_ylabel(r"$V_D/V_{\infty}\ [-]$")
ax1[2,0].set_xlabel(r"$V_{sur}/V_{\infty}\ [-]$")
ax1[2,1].set_title("Case5")
ax1[2,1].plot(p70Case5.surgeVelTheory,p70Case5.discVelDynTheory)
#ax1[2,1].plot(p70Case5.surgeVelTheory,p70Case5.discVelDynTheoryOld)
ax1[2,1].plot(p70Case5.surgeVel,p70Case5.discVelDyn[40],marker="x",markevery=defaultLocs,linestyle=":",color="k")
ax1[2,1].plot(p70Case5.surgeVel[2],p70Case5.discVelDyn[40][2],marker="o",markerfacecolor="g", color="k")
ax1[2,1].plot(p70Case5.surgeVel[7],p70Case5.discVelDyn[40][7],marker="o",markerfacecolor="r", color="k")
ax1[2,1].set_xlabel(r"$V_{sur}/V_{\infty}\ [-]$")

fig2, ax2 = plt.subplots(ncols=2)
fig2.set_size_inches(26/2.54,10/2.54)
ax2[0].plot(p70Case0.surgeVelTheory/max(p70Case0.surgeVelTheory),p70Case0.discVelDynTheory,label="Case0")
ax2[0].plot(p70Case1.surgeVelTheory/max(p70Case1.surgeVelTheory),p70Case1.discVelDynTheory,label="Case1")
ax2[0].plot(p70Case2.surgeVelTheory/max(p70Case2.surgeVelTheory),p70Case2.discVelDynTheory,label="Case2")
ax2[0].plot(p70Case3.surgeVelTheory/max(p70Case3.surgeVelTheory),p70Case3.discVelDynTheory,label="Case3")
ax2[0].plot(p70Case4.surgeVelTheory/max(p70Case4.surgeVelTheory),p70Case4.discVelDynTheory,label="Case4")
ax2[0].plot(p70Case5.surgeVelTheory/max(p70Case5.surgeVelTheory),p70Case5.discVelDynTheory,label="Case5")
ax2[0].legend()
ax2[0].grid()
ax2[0].set_ylabel(r"$V_D/V_{\infty}\ [-]$")
ax2[0].set_xlabel(r"$V_{sur}/\max (V_{sur})\ [-]$")
ax2[1].plot(p70Case0.surgeVel[:3]/max(p70Case0.surgeVel), p70Case0.discVelDyn[40][:3], color=ax2[0].lines[0].get_color(), linestyle = "--", label = r"$270\degree\ to\ 90\degree phases$")
ax2[1].plot(p70Case0.surgeVel[2:8]/max(p70Case0.surgeVel), p70Case0.discVelDyn[40][2:8], color=ax2[0].lines[0].get_color(), label = r"$90\degree\ to\ 270\degree phases$")
ax2[1].plot(p70Case0.surgeVel[7:]/max(p70Case0.surgeVel), p70Case0.discVelDyn[40][7:], color=ax2[0].lines[0].get_color(), linestyle = "--")
ax2[1].plot(p70Case1.surgeVel[:3]/max(p70Case1.surgeVel), p70Case1.discVelDyn[40][:3], color=ax2[0].lines[1].get_color(), linestyle = "--")
ax2[1].plot(p70Case1.surgeVel[2:8]/max(p70Case1.surgeVel), p70Case1.discVelDyn[40][2:8], color=ax2[0].lines[1].get_color())
ax2[1].plot(p70Case1.surgeVel[7:]/max(p70Case1.surgeVel), p70Case1.discVelDyn[40][7:], color=ax2[0].lines[1].get_color(), linestyle = "--")
ax2[1].plot(p70Case2.surgeVel[:3]/max(p70Case2.surgeVel), p70Case2.discVelDyn[40][:3], color=ax2[0].lines[2].get_color(), linestyle = "--")
ax2[1].plot(p70Case2.surgeVel[2:8]/max(p70Case2.surgeVel), p70Case2.discVelDyn[40][2:8], color=ax2[0].lines[2].get_color())
ax2[1].plot(p70Case2.surgeVel[7:]/max(p70Case2.surgeVel), p70Case2.discVelDyn[40][7:], color=ax2[0].lines[2].get_color(), linestyle = "--")
ax2[1].plot(p70Case3.surgeVel[:3]/max(p70Case3.surgeVel), p70Case3.discVelDyn[40][:3], color=ax2[0].lines[3].get_color(), linestyle = "--")
ax2[1].plot(p70Case3.surgeVel[2:8]/max(p70Case3.surgeVel), p70Case3.discVelDyn[40][2:8], color=ax2[0].lines[3].get_color())
ax2[1].plot(p70Case3.surgeVel[7:]/max(p70Case3.surgeVel), p70Case3.discVelDyn[40][7:], color=ax2[0].lines[3].get_color(), linestyle = "--")
ax2[1].plot(p70Case4.surgeVel[:3]/max(p70Case4.surgeVel), p70Case4.discVelDyn[40][:3], color=ax2[0].lines[4].get_color(), linestyle = "--")
ax2[1].plot(p70Case4.surgeVel[2:8]/max(p70Case4.surgeVel), p70Case4.discVelDyn[40][2:8], color=ax2[0].lines[4].get_color())
ax2[1].plot(p70Case4.surgeVel[7:]/max(p70Case4.surgeVel), p70Case4.discVelDyn[40][7:], color=ax2[0].lines[4].get_color(), linestyle = "--")
ax2[1].plot(p70Case5.surgeVel[:3]/max(p70Case5.surgeVel), p70Case5.discVelDyn[40][:3], color=ax2[0].lines[5].get_color(), linestyle = "--")
ax2[1].plot(p70Case5.surgeVel[2:8]/max(p70Case5.surgeVel), p70Case5.discVelDyn[40][2:8], color=ax2[0].lines[5].get_color())
ax2[1].plot(p70Case5.surgeVel[7:]/max(p70Case5.surgeVel), p70Case5.discVelDyn[40][7:], color=ax2[0].lines[5].get_color(), linestyle = "--")
ax2[1].set_xlabel(r"$V_{sur}/\max (V_{sur})\ [-]$")
ax2[1].grid()
ax2[1].legend()

fig3, ax3 = plt.subplots(nrows=3,ncols=2)
fig3.set_size_inches(26/2.54, 30/2.54)
ax3[0,0].set_title("Case0")
ax3[0,0].plot(p70Case0.phasesTheory/(2*np.pi),p70Case0.dynCtTheory,label="quasi-steady")
ax3[0,0].plot(p70Case0.phases/(2*np.pi),p70Case0.dynCt,label="piv-inflow")
ax3[0,0].plot(p70Case0.phases/(2*np.pi),p70Case0.dynCtPiv,label="piv-noca")
ax3[0,0].set_ylabel(r"$C_T\ [-]$")
ax3[0,0].grid()
ax3[0,0].legend()
ax3[0,1].set_title("Case1")
ax3[0,1].plot(p70Case1.phasesTheory/(2*np.pi),p70Case1.dynCtTheory,label="quasi-steady")
ax3[0,1].plot(p70Case1.phases/(2*np.pi),p70Case1.dynCt,label="piv-inflow")
ax3[0,1].grid()
ax3[0,1].plot(p70Case1.phases/(2*np.pi),p70Case1.dynCtPiv,label="piv-noca")
ax3[1,0].set_title("Case2")
ax3[1,0].plot(p70Case2.phasesTheory/(2*np.pi),p70Case2.dynCtTheory,label="quasi-steady")
ax3[1,0].plot(p70Case2.phases/(2*np.pi),p70Case2.dynCt,label="piv-inflow")
ax3[1,0].plot(p70Case2.phases/(2*np.pi),p70Case2.dynCtPiv,label="piv-noca")
ax3[1,0].set_ylabel(r"$C_T\ [-]$")
ax3[1,0].grid()
ax3[1,1].set_title("Case3")
ax3[1,1].plot(p70Case3.phasesTheory/(2*np.pi),p70Case3.dynCtTheory,label="quasi-steady")
ax3[1,1].plot(p70Case3.phases/(2*np.pi),p70Case3.dynCt,label="piv-inflow")
ax3[1,1].plot(p70Case3.phases/(2*np.pi),p70Case3.dynCtPiv,label="piv-noca")
ax3[1,1].grid()
ax3[2,0].set_title("Case4")
ax3[2,0].plot(p70Case4.phasesTheory/(2*np.pi),p70Case4.dynCtTheory,label="quasi-steady")
ax3[2,0].plot(p70Case4.phases/(2*np.pi),p70Case4.dynCt,label="piv-inflow")
ax3[2,0].plot(p70Case4.phases/(2*np.pi),p70Case4.dynCtPiv,label="piv-noca")
ax3[2,0].set_ylabel(r"$C_T\ [-]$")
ax3[2,0].set_xlabel(r"$t/T\ [-]$")
ax3[2,0].grid()
ax3[2,1].set_title("Case5")
ax3[2,1].plot(p70Case5.phasesTheory/(2*np.pi),p70Case5.dynCtTheory,label="quasi-steady")
ax3[2,1].plot(p70Case5.phases/(2*np.pi),p70Case5.dynCt,label="piv-inflow")
ax3[2,1].plot(p70Case5.phases/(2*np.pi),p70Case5.dynCtPiv,label="piv-noca")
ax3[2,1].set_xlabel(r"$t/T\ [-]$")
ax3[2,1].grid()

fig4,ax4 = plt.subplots()
ax4.plot(p70Case0.phases/(2*np.pi),p70Case0.dynCt,marker="x",label="case0")
ax4.plot(p70Case1.phases/(2*np.pi),p70Case1.dynCt,marker="x",label="case1")
ax4.plot(p70Case2.phases/(2*np.pi),p70Case2.dynCt,marker="x",label="case2")
ax4.plot(p70Case3.phases/(2*np.pi),p70Case3.dynCt,marker="x",label="case3")
ax4.plot(p70Case4.phases/(2*np.pi),p70Case4.dynCt,marker="x",label="case4")
ax4.plot(p70Case5.phases/(2*np.pi),p70Case5.dynCt,marker="x",label="case5")
ax4.axhline(0.56,color="k",linestyle="--",label="static")
ax4.set_ylabel(r"$C_T\ [-]$")
ax4.set_xlabel(r"$t/T\ [-]$")
plt.legend()
plt.grid()


            
        
        
        
            
            
            






