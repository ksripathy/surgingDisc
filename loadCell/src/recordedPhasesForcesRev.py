'''Libraries for processing loads data obtained using the surging discs tests conducted in 2024'''

import numpy as np
from operator import itemgetter
from scipy.interpolate import CubicSpline

from src.miscTools import cycleAvg
from src.miscTools import cycleAvgStd
from src.miscTools import lowPassFilt
from src.miscTools import freqsSignal
from src.miscTools import harmonicPassFilt

class recordedData:
    
    def __init__(self, loadsPath, logPath):
        
        self.data = np.loadtxt(loadsPath)
        logData = np.loadtxt(logPath)
        
        self.rhoInf = logData[10]
        self.fstreamVel = logData[12]
        self.area = np.pi * 0.1**2
        
        self.dataExtd, self.cycleData = self.cycleSplit(self.data)
        self.dataExtdClean, self.cycleDataClean = self.cleanRawData(self.dataExtd, self.cycleData, skipBeginCycles=1)
        
        self.dftFreqs, surgeSPD = freqsSignal(self.dataExtdClean[:,0] - np.mean(self.dataExtdClean[:,0]),2e3)
        self.surgeFreq = self.dftFreqs[np.argmax(np.abs(surgeSPD))]
        
    def cycleSplit(self, data, isDataExtd=False):
        
        cycleInitialIndices = [0]
        cumSurgePhases = np.array(data[:,0])
        
        cycleCount = 1
        
        for i in range(1,np.shape(data)[0]):
            
            if data[i-1,0] - data[i,0] > 300:
                    
                    cycleInitialIndices.append(i)
                    cumSurgePhases[i:] = cumSurgePhases[i:] + 360
                    cycleCount = cycleCount + 1
        
        if isDataExtd:
            
            cumSurgePhases = data[:,1] #Replace cumSurgePhases in order to prevent addition by 360
            dataExtd = data
            
        else:
        
            dataExtd = np.insert(data, 1, cumSurgePhases, axis=-1)
            
        lastIndex = np.shape(data)[0]
        cycleInitialIndices.append(lastIndex)
        
        cycleQty = cycleCount
        # print("Cycle Qty:", cycleQty)
        
        cycleData = np.empty((cycleQty,7),dtype = np.ndarray)
        #[cycleInitIndices, cycleSize, phase, cumPhase, loads0, loads1, loads2]
        
        for i in range(cycleQty):
            
            # print(i)
            
            cycleData[i,0] = cycleInitialIndices[i]
            cycleData[i,2] = dataExtd[cycleInitialIndices[i]:cycleInitialIndices[i+1],0]
            cycleData[i,3] = cumSurgePhases[cycleInitialIndices[i]:cycleInitialIndices[i+1]]
            cycleData[i,4] = dataExtd[cycleInitialIndices[i]:cycleInitialIndices[i+1],2]
            cycleData[i,5] = dataExtd[cycleInitialIndices[i]:cycleInitialIndices[i+1],3]
            cycleData[i,6] = dataExtd[cycleInitialIndices[i]:cycleInitialIndices[i+1],4]
            cycleData[i,1] = len(cycleData[i,2])
               
        return dataExtd, cycleData
    
    def cleanRawData(self, dataExtd, cycleData, skipBeginCycles, skipEndCycles=1): #Trim incomplete cycles in the start and end of data. Remove duplicate phases if any
        
        beginTrim = cycleData[skipBeginCycles,0]
        endTrim = cycleData[-skipEndCycles,0]
        dataExtdTrim = dataExtd[beginTrim:endTrim,:]
        self.dataExtdTrim = dataExtdTrim
        
        duplicateIndices = np.argwhere(np.diff(dataExtdTrim[:,1]) <= 0)
        self.duplicateIndices = duplicateIndices
        print("Number of duplicate entries deleted:", len(duplicateIndices))
        
        dataExtdClean = np.delete(dataExtdTrim,duplicateIndices, axis=0)
        cycleDataClean = self.cycleSplit(dataExtdClean, isDataExtd=True)[1]
        
        return dataExtdClean, cycleDataClean 
    
    def resampleData(self,dataExtdClean,cycleDataClean,cycleQty,resampledPhases):
        

        cycleSize = len(resampledPhases)
        
        dataExtdIntp = np.zeros((cycleQty * cycleSize,5),dtype=np.ndarray)
        
        definedPhases = np.tile(resampledPhases,reps=cycleQty)
        
        definedCumPhases = np.array(definedPhases)
        for i in range(cycleQty * cycleSize):
            
            if i % cycleSize == 0:
                
                definedCumPhases[i:] = definedCumPhases[i:] + 360
                
        self.definedCumPhases = definedCumPhases
        
        dataExtdIntp[:,0] = definedPhases
        dataExtdIntp[:,1] = definedCumPhases
        dataExtdIntp[:,2] = np.interp(definedCumPhases, dataExtdClean[:,1], dataExtdClean[:,2])
        dataExtdIntp[:,3] = np.interp(definedCumPhases, dataExtdClean[:,1], dataExtdClean[:,3])
        dataExtdIntp[:,4] = np.interp(definedCumPhases, dataExtdClean[:,1], dataExtdClean[:,4])
        
        cycleDataIntp = self.cycleSplit(dataExtdIntp, isDataExtd=True)[1]
        
        return dataExtdIntp, cycleDataIntp
                
    def cycleAveraging(self, cycleDataIntp):
        
        cycleDataAvg = np.empty(4,dtype=np.ndarray)
        cycleDataStd = np.empty(4,dtype=np.ndarray)
        
        cycleDataAvg[0], cycleDataStd[0] = cycleAvgStd(cycleDataIntp[:,2])
        cycleDataAvg[1], cycleDataStd[1] = cycleAvgStd(cycleDataIntp[:,4])
        cycleDataAvg[2], cycleDataStd[2] = cycleAvgStd(cycleDataIntp[:,5])
        cycleDataAvg[3], cycleDataStd[3] = cycleAvgStd(cycleDataIntp[:,6])
        
        return cycleDataAvg, cycleDataStd
        
    def nonDimensionalize(self, dataExtd):
        
        nonDimData = np.empty(np.shape(dataExtd), dtype=np.ndarray)
        
        nonDimData[:,0] = dataExtd[:,0]
        nonDimData[:,1] = dataExtd[:,1]/360
        nonDimData[:,2] = dataExtd[:,2]/(0.5 * self.rhoInf * self.fstreamVel**2 * self.area)
        nonDimData[:,3] = dataExtd[:,3]/(0.5 * self.rhoInf * self.fstreamVel**2 * self.area)
        nonDimData[:,4] = dataExtd[:,4]/(0.5 * self.rhoInf * self.fstreamVel**2 * self.area)
        
        return nonDimData
    
    def filterData(self, dataExtd):
        
        dftFreqs, surgeSPD = freqsSignal(dataExtd[:,0] - np.mean(dataExtd[:,0]),2e3)
        surgeFreq = dftFreqs[np.argmax(np.abs(surgeSPD))]
        
        filtdDataExtd = np.array(dataExtd)
        
        filtdDataExtd[:,2] = harmonicPassFilt(dataExtd[:,2],2e3,surgeFreq,harmonicsQty=self.maxHarmonic,nonWholeCrit=0.01)[1]
        filtdDataExtd[:,3] = harmonicPassFilt(dataExtd[:,3],2e3,surgeFreq,harmonicsQty=self.maxHarmonic,nonWholeCrit=0.01)[1]
        filtdDataExtd[:,4] = harmonicPassFilt(dataExtd[:,4],2e3,surgeFreq,harmonicsQty=self.maxHarmonic,nonWholeCrit=0.01)[1]
        
        cycleFiltdData = self.cycleSplit(filtdDataExtd, isDataExtd=True)[1]
        
        return filtdDataExtd, cycleFiltdData
        
            
            
        
        
        
        
            
        