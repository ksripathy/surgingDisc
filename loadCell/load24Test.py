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

discPor = "70" #in percentage
kword = "RPM"

loadFileList2N = sorted(glob.glob(dataDir + f"/loadCellTest/2N/p{discPor}Load{kword}*"))
logFileList2N = sorted(glob.glob(dataDir + f"/loadCellTest/2N/p{discPor}Log{kword}*"))

loadFileList10N = sorted(glob.glob(dataDir + f"/loadCellTest/10N/p{discPor}Load{kword}*"))
logFileList10N = sorted(glob.glob(dataDir + f"/loadCellTest/10N/p{discPor}Log{kword}*"))

loadFileList20N = sorted(glob.glob(dataDir + f"/loadCellTest/20N/p{discPor}Load{kword}*"))
logFileList20N = sorted(glob.glob(dataDir + f"/loadCellTest/20N/p{discPor}Log{kword}*"))

loadFileList2NFlip = sorted(glob.glob(dataDir + f"/loadCellFlipTest/2N/flat/p{discPor}Load{kword}*"))
logFileList2NFlip = sorted(glob.glob(dataDir + f"/loadCellFlipTest/2N/flat/p{discPor}Log{kword}*"))

loadFileList10NFlip = sorted(glob.glob(dataDir + f"/loadCellFlipTest/10N/p{discPor}Load{kword}*"))
logFileList10NFlip = sorted(glob.glob(dataDir + f"/loadCellFlipTest/10N/p{discPor}Log{kword}*"))

loadFileList20NFlip = sorted(glob.glob(dataDir + f"/loadCellFlipTest/20N/flat/p{discPor}Load{kword}*"))
logFileList20NFlip = sorted(glob.glob(dataDir + f"/loadCellFlipTest/20N/flat/p{discPor}Log{kword}*"))

loadFileList2NFlipCurved = sorted(glob.glob(dataDir + f"/loadCellFlipTest/2N/curved/p{discPor}Load{kword}*"))
logFileList2NFlipCurved = sorted(glob.glob(dataDir + f"/loadCellFlipTest/2N/curved/p{discPor}Log{kword}*"))

loadFileList2NFlipCurved70 = sorted(glob.glob(dataDir + f"/loadCellFlipTest/2N/curved/p70Load{kword}*"))
logFileList2NFlipCurved70 = sorted(glob.glob(dataDir + f"/loadCellFlipTest/2N/curved/p70Log{kword}*"))

def CTRpm(loadFileList,logFileList):
    
    CTArr = np.empty(len(loadFileList),dtype=np.ndarray)
    surgePhaseArr = np.empty(len(loadFileList),dtype=np.ndarray)
    vInfArr = np.zeros(len(loadFileList))
    discArea = np.pi*0.1**2
    CTMean = np.zeros(len(CTArr))
    CTStd = np.zeros(len(CTArr))

    for i, (loadFile,logFile) in enumerate(zip(loadFileList,logFileList)):
        
        loadData = np.loadtxt(loadFile)
        logData = np.loadtxt(logFile)
        
        axialLoad = loadData[:,1]
        surgePhase = loadData[:,0]
        rhoInf = logData[10]
        vInf = logData[12]

        CT = axialLoad/(0.5*rhoInf*vInf**2*discArea)
        
        vInfArr[i] = vInf
        CTArr[i] = CT
        CTMean[i] = np.mean(CT)
        CTStd[i] = np.std(CT)
        surgePhaseArr[i] = surgePhase
            
        
    return vInfArr, CTArr, CTMean, CTStd
    
vInfArr2N, CTArr2N, CTMean2N, CTStd2N = CTRpm(loadFileList2N, logFileList2N)
vInfArr10N, CTArr10N, CTMean10N, CTStd10N = CTRpm(loadFileList10N, logFileList10N)
vInfArr20N, CTArr20N, CTMean20N, CTStd20N = CTRpm(loadFileList20N, logFileList20N)

vInfArr2NFlip, CTArr2NFlip, CTMean2NFlip, CTStd2NFlip = CTRpm(loadFileList2NFlip, logFileList2NFlip)
vInfArr10NFlip, CTArr10NFlip, CTMean10NFlip, CTStd10NFlip = CTRpm(loadFileList10NFlip, logFileList10NFlip)
vInfArr20NFlip, CTArr20NFlip, CTMean20NFlip, CTStd20NFlip = CTRpm(loadFileList20NFlip, logFileList20NFlip)

vInfArr2NFlipCurv, CTArr2NFlipCurv, CTMean2NFlipCurv, CTStd2NFlipCurv = CTRpm(loadFileList2NFlipCurved, logFileList2NFlipCurved)
vInfArr2NFlipCurv70, CTArr2NFlipCurv70, CTMean2NFlipCurv70, CTStd2NFlipCurv70 = CTRpm(loadFileList2NFlipCurved70, logFileList2NFlipCurved70)
    
"Dynamic data processing"
# freq, amp = freqsSignal(surgePhaseArr[caseNo] - np.mean(surgePhaseArr[caseNo]), 2e3)
# surgeFreq = freq[np.argmax(np.abs(amp))]

# freq2, amp2 = freqsSignal(CTArr[caseNo] - np.mean(CTArr[caseNo]), 2e3)
# surgeFreq2 = freq2[np.argmax(np.abs(amp2))]

# fig2, ax2 = plt.subplots()
# ax2.stem(freq[:300],amp[:300])

# fig3, ax3 = plt.subplots()
# ax3.stem(freq2[:300],amp2[:300])

#Static data processing
uc2N = 1.96 * CTStd2N / np.sqrt(len(CTArr2N[0]))
uc10N = 1.96 * CTStd10N / np.sqrt(len(CTArr10N[0]))
uc20N = 1.96 * CTStd20N / np.sqrt(len(CTArr20N[0]))

uc2NFlip = 1.96 * CTStd2NFlip / np.sqrt(len(CTArr2NFlip[0]))
uc10NFlip = 1.96 * CTStd10NFlip / np.sqrt(len(CTArr10NFlip[0]))
uc20NFlip = 1.96 * CTStd20NFlip / np.sqrt(len(CTArr20NFlip[0]))

uc2NFlipCurved = 1.96 * CTStd2NFlipCurv / np.sqrt(len(CTArr2NFlipCurv[0]))
uc2NFlipCurved70 = 1.96 * CTStd2NFlipCurv70 / np.sqrt(len(CTArr2NFlipCurv70[0]))

fig1, ax1 = plt.subplots()
ax1.errorbar(vInfArr2N, CTMean2N, xerr=0.1, yerr=0.2, label = "2N")
ax1.plot(vInfArr10N, CTMean10N, label = "10N")
ax1.plot(vInfArr20N, CTMean20N, label = "20N")
ax1.legend()
ax1.set_xlabel(r"$V_{\infty} [ms^{-1}]$")
ax1.set_ylabel(r"$C_T [-]$")
ax1.grid()

fig3, ax3 = plt.subplots()
ax3.plot(vInfArr2NFlip, CTMean2NFlip, label = "2N")
ax3.plot(vInfArr10NFlip, CTMean10NFlip, label = "10N")
ax3.plot(vInfArr20NFlip, CTMean20NFlip, label = "20N")
ax3.legend()
ax3.set_xlabel(r"$V_{\infty} [ms^{-1}]$")
ax3.set_ylabel(r"$C_T [-]$")
ax3.grid()

fig2, ax2 = plt.subplots()
ax2.errorbar(vInfArr2NFlipCurv,CTMean2NFlipCurv,xerr=0,yerr=uc2NFlipCurved, label="p45Curv")
ax2.grid()
ax2.set_xlabel(r"$V_{\infty} [ms^{-1}]$")
ax2.set_ylabel(r"$C_T [-]$")
#ax2.errorbar(vInfArr2NFlipCurv70,CTMean2NFlipCurv70,xerr=0,yerr=uc2NFlipCurved70, label="p70Curv")
#ax2.legend()

# ax2 = ax1.twinx()
# ax2.plot(vInfArr2N, CTStd2N, "r")
    
    


