import os
import sys
import traceback
import warnings

def warn_with_traceback(message, category, filename, lineno, file=None, line=None):

    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))

# warnings.showwarning = warn_with_traceback

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

import numpy as np
import glob
import matplotlib.pyplot as plt
from matplotlib import ticker

from src.miscTools import freqsSignal
from src.miscTools import freqsSignal2
from src.miscTools import harmonicPassFilt
from src.miscTools import dftPlot
from src.recordedPhasesForces import recordedData
from src.aeroPhasesForces import loadsProcessing

def filePaths(discPor, kword, loadCell):
    
    totalLoadsPath = dataDir + f"/testMatrixTotal/curved/{loadCell}/p{discPor}Load{kword}.txt"
    totalLogPath = dataDir + f"/testMatrixTotal/curved/{loadCell}/p{discPor}Log{kword}.txt"

    inertialLoadsPath = dataDir + f"/testMatrixInertial/{loadCell}/p{discPor}Load{kword}.txt"
    inertialLogPath = dataDir + f"/testMatrixTotal/curved/{loadCell}/p{discPor}Log{kword}.txt"
    
    return totalLoadsPath,totalLogPath,inertialLoadsPath,inertialLogPath

discPor = "45"
loadCell="2N"
cases=["Case18s"]

totalObjs = []
inertialObjs = []
loadsProcObjs = []
aeroObjs = []

for case in cases:
    
    print(case)
    
    totalObj = recordedData(*filePaths(discPor,case,loadCell)[:2])
    totalObjs.append(totalObj)
    
    inertialObj = recordedData(*filePaths(discPor,case,loadCell)[2:])
    inertialObjs.append(inertialObj)
    
    loadsProcObj = loadsProcessing(totalObj,inertialObj)
    loadsProcObjs.append(loadsProcObj)
    
    aeroObj = loadsProcObj.aeroLoadsObj
    aeroObjs.append(aeroObj)
    
# objAttrb = totalObj.loadsData[totalObj.preSkipCycleStartIndices[1]:,0]
# freq, amp = freqsSignal(objAttrb - np.mean(objAttrb),2e3)
# fig2, ax2 = plt.subplots()
# ax2.stem(freq[:1000],amp[:1000])
# ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator(10))

# cycleNo = loadsProcObj.cycleQty - 1
cycleNo = int(0.5*(loadsProcObj.cycleQty - 1))
# cycleNo = 0

rawTotalLoads = totalObj.cycleLoads[cycleNo]
rawInertialLoads = inertialObj.cycleLoads[cycleNo]
rawTotalPhases = totalObj.cycleLoadPhases[cycleNo]
rawInertialPhases = inertialObj.cycleLoadPhases[cycleNo]
definedPhases = aeroObj.cycleLoadPhases[cycleNo]

totalObj.cycleMeanAndUcty(aeroObj.cumLoadPhases,aeroObj.cycleCumLoadPhases)
inertialObj.cycleMeanAndUcty(aeroObj.cumLoadPhases,aeroObj.cycleCumLoadPhases)

intpTotalLoads = totalObj.cycleIntpLoads[cycleNo]
intpInertialLoads = inertialObj.cycleIntpLoads[cycleNo]

totalLoadsMean = totalObj.cycleIntpLoadsMean
totalLoadsUcty = totalObj.cycleIntpLoadsUcty
totalLoadsStd = totalObj.cycleIntpLoadsStd

inertialLoadsStd = inertialObj.cycleIntpLoadsStd
inertialLoadsMean = inertialObj.cycleIntpLoadsMean
inertialLoadsUcty = inertialObj.cycleIntpLoadsUcty

loadsProcObj.computeAeroLoads()
intpAeroLoads = aeroObj.cycleIntpLoads[cycleNo]
aeroObj.cycleMeanAndUcty(aeroObj.cumLoadPhases,aeroObj.cycleCumLoadPhases)

aeroLoadsMean = aeroObj.cycleIntpLoadsMean
aeroLoadsUcty = aeroObj.cycleIntpLoadsUcty

#Interpolation test
fig1, ax1 = plt.subplots()
ax1.plot(rawTotalPhases, rawTotalLoads)
ax1.set_xlabel(r"$\phi [\degree]$")
ax1.set_ylabel(r"Force [N]")
ax1.set_ylim(-2,2)
# ax1.set_ylim(-0.2,0.25)
# ax1.plot(definedPhases, intpTotalLoads, color="tab:red")

fig3, ax3 = plt.subplots()
ax3.plot(rawInertialPhases, rawInertialLoads, color="tab:orange")
ax3.set_xlabel(r"$\phi [\degree]$")
ax3.set_ylabel(r"Force [N]")
ax3.set_ylim(-2,2)
# ax3.set_ylim(-0.20,0.25)

#Mean total and inertial load plot
fig4, ax4 = plt.subplots()
ax4.plot(definedPhases, totalLoadsMean)
ax4.fill_between(definedPhases, totalLoadsMean-totalLoadsUcty, totalLoadsMean+totalLoadsUcty, alpha=0.5)
ax4.plot(definedPhases, inertialLoadsMean)
ax4.fill_between(definedPhases, inertialLoadsMean-inertialLoadsUcty, inertialLoadsMean+inertialLoadsUcty, alpha=0.5)
ax4.set_xlabel(r"$\phi [\degree]$")
ax4.set_ylabel(r"Force [N]")
ax4.set_ylim(-2,2)
# ax4.set_ylim(-0.2,0.25)

#Aero loads
fig5, ax5 = plt.subplots()
ax5.plot(definedPhases, intpAeroLoads, "tab:green")
ax5.set_xlabel(r"$\phi [\degree]$")
ax5.set_ylabel(r"Force [N]")

#Mean aero laods and uncertainty
fig6, ax6 = plt.subplots()
ax6.plot(definedPhases, aeroLoadsMean, "tab:green")
# ax6.plot(definedPhases,totalLoadsMean-inertialLoadsMean, "tab:red")
# ax6.fill_between(definedPhases, totalLoadsMean-inertialLoadsMean-(totalLoadsUcty+inertialLoadsUcty), totalLoadsMean-inertialLoadsMean+(totalLoadsUcty+inertialLoadsUcty), alpha=0.5, color="tab:red")
ax6.fill_between(definedPhases, aeroLoadsMean-aeroLoadsUcty, aeroLoadsMean+aeroLoadsUcty, alpha=0.5, color="tab:green")
ax6.set_xlabel(r"$\phi [\degree]$")
ax6.set_ylabel(r"Force [N]")
ax6.set_ylim(-1,1)
# ax6.set_ylim(-0.25,0.35)

recordedPhases = inertialObj.loadsData[inertialObj.cycleStartIndices[0]:inertialObj.cycleStartIndices[-1],0]
surgeFreqs, surgeSPD = freqsSignal(recordedPhases - np.mean(recordedPhases),2e3)
surgeFreq = surgeFreqs[np.argmax(np.abs(surgeSPD))]
surgeFreqHarmonics = np.arange(1,21) * surgeFreq

unpackedLoads = getattr(totalObj, "intpLoads")
freqLimit = 5
# freqLimit = 2001
freq, amp = freqsSignal(unpackedLoads - np.mean(unpackedLoads),2e3)
fig7, ax7 = plt.subplots()
ax7.stem(freq[np.argwhere(freq<freqLimit)], amp[np.argwhere(freq<freqLimit)],markerfmt="tab:green",linefmt="tab:green")
ax7.vlines(surgeFreqHarmonics[np.argwhere(surgeFreqHarmonics<freqLimit)], ymin=min(amp[np.argwhere(freq<freqLimit)]), ymax=max(amp[np.argwhere(freq<freqLimit)]), color="k", linestyle="--")

fig8, ax8 = plt.subplots()
ax8.scatter(totalObj.loadsData[totalObj.cycleStartIndices[0]:totalObj.cycleStartIndices[-1],0], totalObj.loadsData[totalObj.cycleStartIndices[0]:totalObj.cycleStartIndices[-1],1])

freq2, amp2 = freqsSignal2(unpackedLoads - np.mean(unpackedLoads),2e3)

fig9, ax9 = plt.subplots()
ax9.stem(freq2[np.argwhere(freq2<freqLimit)], amp2[np.argwhere(freq2<freqLimit)],markerfmt="tab:green",linefmt="tab:green")
#ax9.vlines(surgeFreqHarmonics[np.argwhere(surgeFreqHarmonics<freqLimit)], ymin=min(amp2[np.argwhere(freq<freqLimit)]), ymax=max(amp2[np.argwhere(freq<freqLimit)]), color="k", linestyle="--")

cutOffFreqsIndxHalf = np.argwhere(np.logical_and(freq % surgeFreq < 5e-2,freq <2001))[3:]
cutOffFreqsIndx = np.concatenate((cutOffFreqsIndxHalf, np.flip(len(freq2) - cutOffFreqsIndxHalf)))
cutOffFreqs = freq2[cutOffFreqsIndx]

allIndx = np.arange(len(freq2))
complmntInd = np.setdiff1d(allIndx, cutOffFreqsIndx, assume_unique=True)

filtAmp2 = np.array(amp2)
filtAmp2[complmntInd] = 0

filtLoads = np.fft.ifft(filtAmp2) + np.mean(unpackedLoads)

cycleSlice = np.arange(cycleNo*(loadsProcObj.cycleSize),(cycleNo+1)*(loadsProcObj.cycleSize))

filtFreqs, filtForce = harmonicPassFilt(unpackedLoads, 2000, surgeFreq, 10, nonWholeCrit=0.02)

#Pre harmonic filter spectrum
fig11, ax11 = plt.subplots()
dftPlot(unpackedLoads, 2000, ax11, 30*surgeFreq, surgeFreq)

#Post harmonic filter spectrum
fig12, ax12 = plt.subplots()
dftPlot(filtForce, 2000, ax12, 30*surgeFreq, surgeFreq, color="tab:orange")

#Filterd loads plot
fig10, ax10 = plt.subplots()
ax10.plot(definedPhases, unpackedLoads[cycleSlice])
ax10.plot(definedPhases, filtForce[cycleSlice])












# filtAmp2Left = np.real(amp2[cutOffFreqsIndx])
# filtAmp2Right = np.real(amp2[-cutOffFreqsIndx])

