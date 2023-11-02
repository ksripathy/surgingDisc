from nptdms import TdmsFile
import numpy as np
import scipy.signal as signal
import scipy.stats as stats
import copy

from miscTools import lagShift
from miscTools import lowPassFilt
from miscTools import cycleAvg
from miscTools import cycleAvgStd
from miscTools import nonCycleTrim
from miscTools import symCycleTrim
from miscTools import freqsSignal

class recordedForces:
    
    def __init__(self, inputFilePath, fstreamVelc, fstreamRho, discDia, approxSurgeFreq=0):
        
        self.generateMetaData(inputFilePath)
        
        dynamicPres = 0.5 * fstreamRho * fstreamVelc**2 * 0.25 * np.pi * discDia**2
        
        inputFile = TdmsFile.read(inputFilePath)
        self.forceSeries = inputFile.groups()[0].channels()[1][:]/dynamicPres
        self.totalSampleSize = len(self.forceSeries)
        
        self.timeSeries = np.array([4e-5*i for i in range(self.totalSampleSize)])
        self.deltaTime = self.timeSeries[1] - self.timeSeries[0]
        self.samplingFreq = 1/self.deltaTime
        self.totalTime = self.totalSampleSize/self.samplingFreq
        
        if approxSurgeFreq != 0:      
            self.surgeFreq = approxSurgeFreq
            
        else:    
            self.surgeFreq = self.computeSurgeFreq()
        
        self.filtForce = self.filtSignal()
        self.surgePeriods = self.computeSurgePeriod()
        
    def computeSurgeFreq(self):    
    
        freqs, amps = freqsSignal(self.forceSeries, self.samplingFreq)
        
        return freqs[np.argmax(np.abs(amps[1:])) + 1] #Amplitude peak at 0 freq is not considered 
    
    def filtSignal(self):
        
        cutOffFreq = np.ceil(self.surgeFreq)
        
        return lowPassFilt(self.forceSeries, self.samplingFreq, cutOffFreq)
    
    def computeSurgePeriod(self):
        
        peakIndices = signal.find_peaks(self.filtForce)[0]
        peakTimes = self.timeSeries[peakIndices]
        
        surgePeriods = peakTimes[1:] - np.roll(peakTimes,1)[1:]
        
        return surgePeriods
    
    def generateMetaData(self, inputFilePath):
        
        str1, str2 = inputFilePath.split("/")[-2:]
        
        str2 = str2.split(".")[0]
        self.uid = copy.deepcopy(str2)
        
        str2 = str2.replace("clean", "")
        
        if str1 == "dynamic":
            
            self.recordType = "Combined Loading"
            
        elif str1 == "clean":
            
            self.recordType = "Inertial Loading"
            
        elif str1 == "static":
            
            self.recordType = "Static Loading"
            
        self.porosity = str2[1:3]
        
        if str1 == "dynamic" or str1 == "clean":
        
            if "vmax1" not in str2:
                
                self.reducedSurgeFreq = float(str2[4])
                surgeAmp = float(str2[6:].replace("p", "."))
                self.reducedSurgeAmp = surgeAmp / 20 #Non dimensionalizing with disc dia
                
            else:
                
                self.reducedSurgeFreq = 2.7
                self.reducedSurgeAmp = 0.375
                
        elif str1 == "static":
            
            self.fstreamVelc = str2[-1]
                
class aeroForces:
    
    def __init__(self, windObj, noWindObj, trimCycles):
                
        self.surgePeriod = np.mean(windObj.surgePeriods[trimCycles:-trimCycles])
        self.surgeCycleSize = round(self.surgePeriod/windObj.deltaTime)
        self.surgeCycles = min(int(len(windObj.timeSeries) / self.surgeCycleSize), int(len(noWindObj.timeSeries) / self.surgeCycleSize)) #In case mismatch in length of signals
        self.trimSurgeCycles = self.surgeCycles - (2 * trimCycles)
        
        setattr(windObj, "forceAsymTrim", nonCycleTrim(windObj.forceSeries, self.surgeCycleSize, self.surgeCycles))
        setattr(windObj, "filtForceAsymTrim", nonCycleTrim(windObj.filtForce, self.surgeCycleSize, self.surgeCycles))
        setattr(noWindObj, "forceAsymTrim", nonCycleTrim(noWindObj.forceSeries, self.surgeCycleSize, self.surgeCycles))
        setattr(noWindObj, "filtForceAsymTrim", nonCycleTrim(noWindObj.filtForce, self.surgeCycleSize, self.surgeCycles))
        
        setattr(windObj, "forceTrim", symCycleTrim(windObj.forceAsymTrim, self.surgeCycleSize, trimCycles))
        setattr(windObj, "filtForceTrim", symCycleTrim(windObj.filtForceAsymTrim, self.surgeCycleSize, trimCycles))
        setattr(noWindObj, "forceTrim", symCycleTrim(noWindObj.forceAsymTrim, self.surgeCycleSize, trimCycles))
        setattr(noWindObj, "filtForceTrim", symCycleTrim(noWindObj.filtForceAsymTrim, self.surgeCycleSize, trimCycles))
        
        sigOffset = lagShift(windObj.filtForceTrim, noWindObj.filtForceTrim)[0]
        self.sigOffset = sigOffset
        self.sigOffsetRaw = lagShift(windObj.forceSeries, noWindObj.forceSeries)[0]
        setattr(noWindObj, "syncForceSeries", np.roll(noWindObj.forceSeries, sigOffset))
        setattr(noWindObj, "syncFiltForce", np.roll(noWindObj.filtForce, sigOffset))
        
        setattr(noWindObj, "syncForceAsymTrim", nonCycleTrim(noWindObj.syncForceSeries, self.surgeCycleSize, self.surgeCycles))
        setattr(noWindObj, "syncFiltForceAsymTrim", nonCycleTrim(noWindObj.syncFiltForce, self.surgeCycleSize, self.surgeCycles))
        
        setattr(noWindObj, "syncForceTrim", symCycleTrim(noWindObj.syncForceAsymTrim, self.surgeCycleSize, trimCycles))
        setattr(noWindObj, "syncFiltForceTrim", symCycleTrim(noWindObj.syncFiltForceAsymTrim, self.surgeCycleSize, trimCycles))
        
        timeAsymTrim = nonCycleTrim(windObj.timeSeries, self.surgeCycleSize, self.surgeCycles)
        self.timeSeries = symCycleTrim(timeAsymTrim, self.surgeCycleSize, trimCycles)
        
        self.forceSeries = windObj.forceTrim - noWindObj.syncForceTrim
        #Multi pass filtering for better noise removal
        self.filtForceV2temp = lowPassFilt(self.forceSeries, windObj.samplingFreq, np.ceil(windObj.surgeFreq))
        self.filtForceV2 = lowPassFilt(self.filtForceV2temp, windObj.samplingFreq, np.ceil(windObj.surgeFreq))
        #Filtering aero data for better plot visualisation 
        self.filtForce = lowPassFilt(windObj.filtForceTrim - noWindObj.syncFiltForceTrim, windObj.samplingFreq, np.ceil(windObj.surgeFreq))
        self.filtForceV3 = lowPassFilt(windObj.filtForceTrim - noWindObj.syncFiltForceTrim - np.mean(noWindObj.syncForceTrim), windObj.samplingFreq, np.ceil(windObj.surgeFreq))
        
        self.cycleNdimTime = np.linspace(0,1,num=self.surgeCycleSize)
        self.cycleNdimTimeRe = np.linspace(0,1,20000)
        self.cycleTime = self.genCycles(self.timeSeries)
        self.cycleForce = self.genCycles(self.forceSeries)
        self.cycleFiltForce = self.genCycles(self.filtForce)
        self.cycleFiltForceV2 = self.genCycles(self.filtForceV2)
        self.cycleFiltForceV3 = self.genCycles(self.filtForceV3)
        
        setattr(windObj, "cycleForce", self.genCycles(windObj.forceTrim))
        setattr(windObj, "cycleFiltForce", self.genCycles(windObj.filtForceTrim))
        setattr(noWindObj, "cycleSyncForce", self.genCycles(noWindObj.syncForceTrim))
        setattr(noWindObj, "cycleSyncFiltForce", self.genCycles(noWindObj.syncFiltForceTrim))
        
        self.cycleAvgForce, self.cycleStdForce = cycleAvgStd(self.cycleForce)
        self.cycleAvgFiltForce = cycleAvg(self.cycleFiltForce)
        self.cycleAvgFiltForceRe = np.interp(self.cycleNdimTimeRe, self.cycleNdimTime, self.cycleAvgFiltForce)
        self.cycleAvgFiltForceV2 = cycleAvg(self.cycleFiltForceV2)
        self.cycleAvgFiltForceV3 = cycleAvg(self.cycleFiltForceV3)
        setattr(windObj, "cycleAvgForce", cycleAvg(windObj.cycleForce))
        setattr(windObj, "cycleAvgFiltForce", cycleAvg(windObj.cycleFiltForce))
        setattr(noWindObj, "cycleAvgSyncForce", cycleAvg(noWindObj.cycleSyncForce))
        setattr(noWindObj, "cycleAvgSyncFiltForce", cycleAvg(noWindObj.cycleSyncFiltForce))
        
    def genCycles(self, sigTrim):
        
        cycleSig = np.empty(self.trimSurgeCycles, dtype=np.ndarray)
        
        for cycle in range(self.trimSurgeCycles):
            
            cycleIter = np.arange(cycle * self.surgeCycleSize, (cycle+1) * self.surgeCycleSize)
            cycleSig[cycle] = sigTrim[cycleIter]
            
        return cycleSig
        
        
        
        
        