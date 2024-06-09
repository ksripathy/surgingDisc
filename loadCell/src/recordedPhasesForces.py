'''Libraries for processing loads data obtained using the surging discs tests conducted in 2024'''

import numpy as np
from operator import itemgetter
from scipy.interpolate import CubicSpline

from src.miscTools import cycleAvg
from src.miscTools import cycleAvgStd
from src.miscTools import lowPassFilt
from src.miscTools import freqsSignal

class recordedData:
    
    def __init__(self, loadsPath, logPath, skipCycles=10):
        
        self.loadsData = np.loadtxt(loadsPath)
        logData = np.loadtxt(logPath)
        self.skipCycles = skipCycles
        
        self.rhoInf = logData[10]
        self.fstreamVel = logData[12]
        
        self.cycleStartIndices, self.cycleLoadPhases, self.cycleLoads, self.cumLoadPhases, self.preSkipCycleStartIndices = self.cycleSplit()
        
        self.loadsDataClean, self.cumLoadPhasesClean = self.cleanRawData()
        
        self.dftFreqs, surgeSPD = freqsSignal(self.loadsDataClean[:,0] - np.mean(self.loadsDataClean[:,0]),2e3)
        self.surgeFreq = self.dftFreqs[np.argmax(np.abs(surgeSPD))]
        
    def cycleSplit(self):
        
        preSkipCycleInitialIndices = [0]
        cycleInitialIndices = []
        cycleCount = 1
        cumSurgePhases = np.array(self.loadsData[:,0])
        
        for i in range(1,np.shape(self.loadsData)[0]):
            
            if self.loadsData[i-1,0] - self.loadsData[i,0] > 300:
                
                cycleCount = cycleCount + 1
                cumSurgePhases[i:] = cumSurgePhases[i:] + 360
                
                if cycleCount <= self.skipCycles:
                    
                    preSkipCycleInitialIndices.append(i)
                
                else:
                    
                    cycleInitialIndices.append(i)      
        
        cycleQty = len(cycleInitialIndices) - 1 #Removing last cycle since it is expected to be incomplete
        
        cycleSurgePhases = np.empty(cycleQty, dtype = np.ndarray)
        cycleLoads = np.empty(cycleQty, dtype = np.ndarray)
        
        for i in range(cycleQty):
            
            cycleSurgePhases[i] = self.loadsData[cycleInitialIndices[i]:cycleInitialIndices[i+1],0]
            cycleLoads[i] = self.loadsData[cycleInitialIndices[i]:cycleInitialIndices[i+1],1]
            
            
        return cycleInitialIndices, cycleSurgePhases, cycleLoads, cumSurgePhases, preSkipCycleInitialIndices
    
    def cycleSizes(self):
        
        cycleSizes = []
        
        for i in range(len(self.cycleLoadPhases)):
    
            cycleSizes.append(len(self.cycleLoadPhases[i]))
            
        return cycleSizes
    
    def cleanRawData(self): #Trim incomplete cycles in the start and end of data. Remove duplicate phases if any
        
        beginTrim = self.cycleStartIndices[0]
        endTrim = self.cycleStartIndices[-1]
        
        loadsDataTrim = self.loadsData[beginTrim:endTrim,:]
        cumLoadPhasesTrim = self.cumLoadPhases[beginTrim:endTrim]
        
        duplicateIndices = np.argwhere(np.diff(cumLoadPhasesTrim) <= 0)
        
        print("Number of duplicate entries deleted:", len(duplicateIndices))
        
        loadsDataClean = np.delete(loadsDataTrim, duplicateIndices, axis=0)
        cumLoadPhasesClean = np.delete(cumLoadPhasesTrim, duplicateIndices)
        
        return loadsDataClean, cumLoadPhasesClean
    
    def resampleLoads(self,definedPhases):
        
        recPhases = self.cumLoadPhasesClean
        recLoads = self.loadsDataClean[:,1]
        
        intpLoads = np.interp(definedPhases, recPhases, recLoads)
        return intpLoads    
    
    def cycleMeanAndUcty(self,definedCumLoadPhasesClean,definedCycleCumLoadPhases): 
        
        cycleQty = len(definedCycleCumLoadPhases)
        
        if not hasattr(self,"cycleIntpLoads"): #Exception for aeroLoadsObj. For aeroLoadsObj this attribute will be generated in aeroPhasesForces
            self.cycleIntpLoads = np.empty(cycleQty,dtype=np.ndarray)
        
            for i in range(cycleQty):
                
                definedPhases = definedCycleCumLoadPhases[i]
                self.cycleIntpLoads[i] = self.resampleLoads(definedPhases)
                
            self.intpLoads = self.resampleLoads(definedCumLoadPhasesClean)
            
        self.cycleIntpLoadsMean, self.cycleIntpLoadsStd = cycleAvgStd(self.cycleIntpLoads[:-1]) #Last cycle not considered in averaging to remove potential artifacting from boundary effects
        self.cycleIntpLoadsUcty = 1.96 * self.cycleIntpLoadsStd / np.sqrt(cycleQty-1)
        
    def nonDimLoads(self):
        
        cycleQty = len(self.cycleIntpLoads)
        
        ndimLoads = np.empty(cycleQty, dtype=np.ndarray)
        
        for i in range(cycleQty):
            
            ndimLoads[i] = self.cycleIntpLoads[i]/(0.5 * self.rhoInf * self.fstreamVel**2 * np.pi * 0.10**2)
            
        ndimLoadsMean = self.cycleIntpLoadsMean/(0.5 * self.rhoInf * self.fstreamVel**2 * np.pi * 0.10**2)
        ndimLoadsUcty = self.cycleIntpLoadsUcty/(0.5 * self.rhoInf * self.fstreamVel**2 * np.pi * 0.10**2)
        
        self.ndimCycleIntpLoads = ndimLoads
        self.ndimCycleIntpLoadsMean = ndimLoadsMean
        self.ndimCycleIntpLoadsUcty = ndimLoadsUcty
            
            
        
        
        
        
            
        