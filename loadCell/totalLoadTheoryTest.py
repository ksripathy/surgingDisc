import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import sys

sys.path.append("..")

from src.miscTools import lagShift
from src.miscTools import symAvgArray
from piv.src.ioUtilsNoStitch import loadData

p70Case6 = loadData("/home/ksripathy/phd/surgingDisc/piv/plots/results_p70Case6.json")
p70Case7 = loadData("/home/ksripathy/phd/surgingDisc/piv/plots/results_p70Case7.json")

caseDict = np.array([[1.0, 0.250, 2.387324, 5.0e-2],
                            [2.0, 0.125, 4.774648, 2.5e-2], 
                            [2.0, 0.250, 4.774648, 5.0e-2], 
                            [2.0, 0.375, 3.183099, 7.5e-2], 
                            [3.0, 0.250, 4.774648, 5.0e-2],
                            [2.7, 0.375, 4.297183, 7.5e-2],
                            [0.0, 0.000, 0.000000, 0.0e-2],
                            [0.0, 0.000, 0.000000, 0.0e-2]])

p70Mass = 37e-3
p70Ct = np.array([0.598, 0.598, 0.598, 0.563, 0.563, 0.563, 0.598, 0.563])

p70FstreamVelc = np.array([3.045224, 3.030344, 3.035440, 2.025460, 2.049157, 2.007377, 3.096507, 2.045414])
p70FstreamRho = np.array([1.196899, 1.198764, 1.198306, 1.194673, 1.192900, 1.189482, 1.194618, 1.195081])

time = np.array([4e-5*i for i in range(2*25000)])
area = 0.25 * np.pi * (20e-2)**2
caseNo = 3

if caseNo < 3:
    
    radialLocs = np.array(p70Case6.cycleAxialInd["spanLocs"][40:]) * 20e-2
    axialIndDistr = symAvgArray(p70Case6.cycleAxialInd["intp200"][0])
    axialIndV2 = axialIndDistr[0]
    
else:
    
    radialLocs = np.array(p70Case7.cycleAxialInd["spanLocs"][40:]) * 20e-2
    axialIndDistr = symAvgArray(p70Case7.cycleAxialInd["intp200"][0])
    axialIndV2 = axialIndDistr[0]

Ct = p70Ct[caseNo]
rho = p70FstreamRho[caseNo]
Vinf = p70FstreamVelc[caseNo]
dynPre = 0.5 * rho * Vinf**2
axialInd = 0.5 * (1 - np.sqrt(1 - Ct))
freq = caseDict[caseNo,2]
amp = caseDict[caseNo,3]
omega = 2 * np.pi * freq
qinf = 0.5*rho*Vinf**2

acclSHM = -amp * omega**2 * np.sin(omega * time)
accl1 = -amp * omega**2 * np.cos(omega * time)
projCrod = np.sqrt(45e-2**2 - (amp * np.sin(omega * time))**2)
accl2 = -amp**2 * omega**2 * (np.cos(omega * time)**2 - np.sin(omega * time)**2)/projCrod
accl3 = -amp**4 * omega**2 * np.sin(omega * time)**2 * np.cos(omega * time)**2 / projCrod**3

velSurgeSHM = amp * omega * np.cos(omega * time)
vel1 = -amp * omega * np.sin(omega * time)
vel2 = -amp**2 * omega * np.sin(omega * time) * np.cos(omega * time)/projCrod
velSurge = vel1 + vel2

inerLoadSHM = p70Mass * acclSHM
inerLoad = p70Mass * (accl1 + accl2 + accl3)

aero1 = Ct * dynPre * area
aeroLoadSHM = aero1 - 2 * rho * area * axialInd * velSurgeSHM * (axialInd * velSurgeSHM + Vinf * (1 - 2*axialInd))
aeroLoad = aero1 - 2 * rho * area * axialInd * velSurge * (axialInd * velSurge + Vinf * (1 - 2*axialInd))

aeroLoadSHMWen = aero1 - 2 * rho * area * (Vinf - velSurgeSHM * (1 - axialInd)) * velSurgeSHM * (1 - axialInd)



aeroLoadSHMAnnuli = np.empty(shape=40, dtype=np.ndarray)
aeroLoadAnnuli = np.empty(shape=40, dtype=np.ndarray)
aeroLoadStatic = 0

for i in range(len(aeroLoadAnnuli)):
    
    axialIndAnn = 0.5*(axialIndDistr[i] + axialIndDistr[i+1])
    areaAnn = np.pi * (radialLocs[i+1]**2 - radialLocs[i]**2)
    
    aero1Ann = 4 * axialIndAnn * (1 - axialIndAnn) * dynPre * areaAnn
    aeroLoadStatic = aero1Ann + aeroLoadStatic
    aeroLoadSHMAnnuli[i] = aero1Ann - 2 * rho * areaAnn * axialIndAnn * velSurgeSHM * (axialIndAnn * velSurgeSHM + Vinf * (1 - 2*axialIndAnn))
    aeroLoadAnnuli[i] = aero1Ann - 2 * rho * areaAnn * axialIndAnn * velSurge * (axialIndAnn * velSurge + Vinf * (1 - 2*axialIndAnn))
    
aeroLoadSHMRev = np.sum(aeroLoadSHMAnnuli)
aeroLoadRev = np.sum(aeroLoadAnnuli)
pivCt = aeroLoadStatic/(qinf * area)

peakIndices = signal.find_peaks(aeroLoadSHM)[0]
peakTimes = time[peakIndices]
        
surgePeriods = peakTimes[1:] - np.roll(peakTimes,1)[1:]

peakIndices2 = signal.find_peaks(inerLoadSHM)[0]
peakTimes2 = time[peakIndices2]
        
surgePeriods2 = peakTimes[1:] - np.roll(peakTimes,1)[1:]

fig1, ax1 = plt.subplots()
ax1.plot(time, aeroLoadSHM/(qinf * area))
#ax1.plot(time, aeroLoadSHMWen/(qinf * area))
ax1.plot(time, aeroLoadSHMRev/(qinf * area))
#ax1.plot(time, aeroLoadSHMAnnuli[0]/(qinf * area))
#ax1.plot(time, aeroLoadSHMAnnuli[9]/(qinf * area))
#ax1.plot(time, inerLoadSHM)
#ax1.plot(time, aeroLoadSHM + inerLoadSHM)

fig2, ax2 = plt.subplots()
ax2.plot(time, aeroLoad/(qinf * area))
ax2.plot(time, aeroLoadRev/(qinf * area))
#ax2.plot(time, inerLoad)
#ax2.plot(time, aeroLoad + inerLoad)

fig3,ax3 = plt.subplots()
ax3.plot(aeroLoadSHM/(qinf * area), velSurgeSHM/Vinf)

print(lagShift(aeroLoadSHMRev+inerLoadSHM, inerLoadSHM))
print(lagShift(aeroLoad+inerLoad, inerLoad))
 



