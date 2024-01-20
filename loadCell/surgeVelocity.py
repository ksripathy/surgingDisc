import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import sys

sys.path.append("..")

from src.miscTools import lagShift
from src.miscTools import symAvgArray
from piv.src.ioUtilsNoStitch import loadData

staticCT = 0.56
rho = 1.195
area = np.pi * (10e-2)**2

caseNo = 0
dynCaseRec = loadData(f"/home/ksripathy/phd/surgingDisc/piv/plots/results_p70Case{caseNo}.json")

caseDict = np.array([[1.0, 0.250, 2.387324, 5.0e-2],
                            [2.0, 0.125, 4.774648, 2.5e-2], 
                            [2.0, 0.250, 4.774648, 5.0e-2], 
                            [2.0, 0.375, 3.183099, 7.5e-2], 
                            [3.0, 0.250, 4.774648, 5.0e-2],
                            [2.7, 0.375, 4.297183, 7.5e-2],
                            [0.0, 0.000, 0.000000, 0.0e-2],
                            [0.0, 0.000, 0.000000, 0.0e-2]])


def discVel(caseNo):

    dynCase = loadData(f"/home/ksripathy/phd/surgingDisc/piv/plots/results_p70Case{caseNo}.json")

    if caseNo < 3:
        
        
        staticCase = loadData("/home/ksripathy/phd/surgingDisc/piv/plots/results_p70Case6.json")
        
    else:
        
        staticCase = loadData("/home/ksripathy/phd/surgingDisc/piv/plots/results_p70Case7.json")
        
    velStaticNdim = 1 - np.array(staticCase.cycleAxialInd["intp200"][0])
    vinfStatic = np.abs(staticCase.phaseVinf)
            
    vinfDyn = np.abs(np.array(dynCase.phaseVinf))
    qinfDyn = 0.5 * rho * vinfDyn**2
    phases = 2*np.pi * np.linspace(0.05, 1.15, 12)
    phasesTheory = 2*np.pi * np.linspace(0,1.15,200)

    velDynamicNdim = np.empty(shape=len(velStaticNdim), dtype=np.ndarray)


    velSurgeExp = np.empty(shape=len(velStaticNdim), dtype=np.ndarray)
    dynCT = np.zeros(len(phases))
    dynCTTheory = np.zeros(len(phasesTheory))
    
    freq = caseDict[caseNo,2]
    amp = caseDict[caseNo,3]
    omega = 2 * np.pi * freq

    for radius in range(len(velStaticNdim)):
        
        velDynamicNdim[radius] = np.zeros(12)
        velSurgeExp[radius] = np.zeros(12)
        
        indTheory = staticCase.cycleAxialInd["intp200"][0][40]
        velSurgeTheory = amp * omega * np.cos(phasesTheory)
        velSurgeTheoryNdim = velSurgeTheory/max(velSurgeTheory)
        velTheory = np.mean(vinfDyn) * (1 - indTheory) + (indTheory * velSurgeTheory)
        velTheoryNdim = velTheory/np.round(vinfStatic)
        dynCTTheory = staticCT + (np.mean(vinfDyn) + (indTheory * velSurgeTheory) - (2 * velTheory)) * (2 * rho * indTheory * velSurgeTheory / np.mean(qinfDyn))
        
        for phase in range(12):
            
            velDynamicNdim[radius][phase] = 1 - dynCase.cycleAxialInd["intp200"][phase][radius]
            #velSurgeExp[radius][phase] = vinfDyn[phase] * (velDynamicNdim[radius][phase] - velStaticNdim[radius])
            
            if radius == 40:
                
                ind = staticCase.cycleAxialInd["intp200"][0][radius]
                velSurge = amp * omega * np.cos(phases[phase])
                velExp = velDynamicNdim[radius][phase] * vinfDyn[phase]
                velTheory = vinfDyn[phase] * (1 - ind) + (ind * velSurge)
                #print("velExp:", velExp)
                dynCT[phase] = staticCT + (vinfDyn[phase] + (ind * velSurge) - (2 * velExp)) * (2 * rho * ind * velSurge / qinfDyn[phase])
                #dynCTTheory[phase] = staticCT + (vinfDyn[phase] + (ind * velSurge) - (2 * velTheory)) * (2 * rho * ind * velSurge / qinfDyn[phase])
            
    velSurgeSHM = amp * omega * np.cos(np.linspace(0,2*np.pi,200))
    velSurgeSHMNdim = velSurgeSHM/max(velSurgeSHM)
    velDiscTheoryNdim = ((velStaticNdim[40] * np.mean(vinfDyn)) + (1 - velStaticNdim[40]) * velSurgeSHM)/np.round(vinfStatic)

    velSurgeSHMReNdim = amp * omega * np.cos(phases)/max(velSurgeSHM)
    velDiscExpNdim = (velDynamicNdim[40] * vinfDyn)/np.round(vinfStatic)
    
    return (vinfStatic, vinfDyn, velSurgeSHMNdim, velSurgeSHMReNdim, velStaticNdim, velDiscTheoryNdim, velDiscExpNdim, phases, phasesTheory, velSurgeTheoryNdim, velTheoryNdim, dynCT, dynCTTheory)


vinfStatic, vinfDyn, velSurgeSHMNdim, velSurgeSHMReNdim, velStaticNdim, velDiscTheoryNdim, velDiscExpNdim, phases, phasesTheory, velSurgeTheoryNdim, velTheoryNdim, dynCTnl, dynCTTheorynl = discVel(caseNo)

fig1, ax1 = plt.subplots()
ax1.plot(velSurgeSHMNdim, velDiscTheoryNdim)
ax1.plot(velSurgeSHMReNdim, velDiscExpNdim)
ax1.plot(0,velStaticNdim[40]*np.mean(vinfDyn)/np.round(vinfStatic),marker="x")

fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()
fig4, ax4 = plt.subplots()
fig5, ax5 = plt.subplots()
fig6, ax6 = plt.subplots()

for case in range(6):
    
    vinfStatic, vinfDyn, velSurgeSHMNdim, velSurgeSHMReNdim, velStaticNdim, velDiscTheoryNdim, velDiscExpNdim, phases, phasesTheory, velSurgeTheoryNdim, velTheoryNdim, dynCT, dynCTTheory = discVel(case)
    ax2.plot(velSurgeSHMNdim, velDiscTheoryNdim)
    ax3.plot(velSurgeSHMReNdim[:3], velDiscExpNdim[:3], color=ax2.get_lines()[case].get_color(), linestyle = "--")
    ax3.plot(velSurgeSHMReNdim[2:8], velDiscExpNdim[2:8], color=ax2.get_lines()[case].get_color())
    ax3.plot(velSurgeSHMReNdim[7:], velDiscExpNdim[7:], color=ax2.get_lines()[case].get_color(), linestyle = "--")
    line = ax4.plot(phases/(2*np.pi), dynCT, label = f"case{case}")
    ax4.plot(phasesTheory/(2*np.pi), dynCTTheory, linestyle="--", color=line[0].get_color())
    line2 = ax5.plot(dynCT, velDiscExpNdim)
    ax5.plot(dynCTTheory, velTheoryNdim, color = line2[0].get_color(), linestyle="--")
    ax6.plot(velSurgeTheoryNdim, dynCTTheory)
    
ax4.legend()
    



