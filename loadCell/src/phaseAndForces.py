'''Libraries for processing loads data obtained using the surging discs tests conducted in 2024'''

import numpy as np
from scipy.interpolate import CubicSpline
from operator import itemgetter

from src.miscTools import cycleAvg
from src.miscTools import cycleAvgStd
from src.miscTools import lowPassFilt
from src.miscTools import highPassFilt
from src.miscTools import freqsSignal
from src.miscTools import bandPassFilt

class recordedData:
    
    def __init__(self, totalLoadsPath, totalLogPath, inertialLoadsPath, inertialLogPath, skipcycles=10):
        
        totalLoadsData = np.loadtxt(totalLoadsPath)
        totalLogData = np.loadtxt(totalLogPath)
        
        totalRhoInf = totalLogData[10]
        self.fstreamVel = totalLogData[12]
            
        inertialLoadsData = np.loadtxt(inertialLoadsPath)
        inertialLogData = np.loadtxt(inertialLogPath)
        
        inertialRhoInf = inertialLogData[10]
        
        cycleIndices1, self.totalCycleLoadPhases, totalCycleLoads = self.cycleSplit(totalLoadsData,skipcycles)[:3]
        cycleIndices2, self.inertialCycleLoadPhases, inertialCycleLoads = self.cycleSplit(inertialLoadsData,skipcycles)[:3]
        
        cycleIndices3, cumTotalLoadPhases = itemgetter(0,3)(self.cycleSplit(totalLoadsData,0))
        cycleIndices4, cumInertialLoadPhases = itemgetter(0,3)(self.cycleSplit(inertialLoadsData,0))
        
        self.cumTotalPhases = cumTotalLoadPhases
        self.cumInertialPhases = cumInertialLoadPhases
        
        print("cycleIndices1:", cycleIndices1)
        print("Average diff cycleIndices1:", np.mean(np.diff(cycleIndices1)))
        print("cycleIndices2:", cycleIndices2)
        print("Average diff cycleIndices2:", np.mean(np.diff(cycleIndices2)))
        
        #Equalize the number of cycles between total and inertial load
        totalCycleQty = len(self.totalCycleLoadPhases)
        inertialCycleQty = len(self.inertialCycleLoadPhases)
        deltaCycleQty = totalCycleQty - inertialCycleQty
        
        print("Delta Cycle Qty:", deltaCycleQty)
        
        if deltaCycleQty > 0:
            
            self.totalCycleLoadPhases = np.delete(self.totalCycleLoadPhases,np.arange(-deltaCycleQty,0))
            totalCycleLoads = np.delete(totalCycleLoads,np.arange(-deltaCycleQty,0))
            
            #Equalize the number of samples in the trimmed total and inertial loads
            self.surgeCycleSize, deltaSampleQty = self.deltaInertialTotalSamples(self.totalCycleLoadPhases,self.inertialCycleLoadPhases)
            
            totalLoadsIndices = np.arange(cycleIndices1[0],cycleIndices1[-1-deltaCycleQty]+deltaSampleQty)
            inertialLoadsIndices = np.arange(cycleIndices2[0],cycleIndices2[-1])
            
            totalLoadsDataTrim = totalLoadsData[cycleIndices1[0]:cycleIndices1[-1-deltaCycleQty]+deltaSampleQty,:]
            inertialLoadsDataTrim = inertialLoadsData[cycleIndices2[0]:cycleIndices2[-1],:]
            
        elif deltaCycleQty < 0:
            
            self.inertialCycleLoadPhases = np.delete(self.inertialCycleLoadPhases,np.arange(deltaCycleQty,0))
            inertialCycleLoads = np.delete(inertialCycleLoads,np.arange(deltaCycleQty,0))
            
            self.surgeCycleSize, deltaSampleQty = self.deltaInertialTotalSamples(self.totalCycleLoadPhases,self.inertialCycleLoadPhases)
            
            totalLoadsIndices = np.arange(cycleIndices1[0],cycleIndices1[-1]+deltaSampleQty)
            inertialLoadsIndices = np.arange(cycleIndices2[0],cycleIndices2[-1+deltaCycleQty])
            
            #This will ensure that totalLoadsData does not run out of samples when deltaSampleQty is huge
            # if deltaSampleQty > len(totalLoadsData[:-1,1]) - cycleIndices1[-1]:
                
            #     self.totalCycleLoadPhases = np.delete(self.totalCycleLoadPhases,-1)
            #     totalCycleLoads = np.delete(totalCycleLoads,-1)
            #     self.inertialCycleLoadPhases = np.delete(self.inertialCycleLoadPhases,-1)
            #     inertialCycleLoads = np.delete(inertialCycleLoads,-1)
                
            #     self.surgeCycleSize, deltaSampleQty = self.deltaInertialTotalSamples(self.totalCycleLoadPhases,self.inertialCycleLoadPhases)
                
            #     totalLoadsIndices = np.arange(cycleIndices1[0],cycleIndices1[-2]+deltaSampleQty)
            #     inertialLoadsIndices = np.arange(cycleIndices2[0],cycleIndices2[-2+deltaCycleQty])
            
            totalLoadsDataTrim = totalLoadsData[cycleIndices1[0]:cycleIndices1[-1]+deltaSampleQty,:]
            inertialLoadsDataTrim = inertialLoadsData[cycleIndices2[0]:cycleIndices2[-1+deltaCycleQty],:]
            
        else:
            
            self.surgeCycleSize, deltaSampleQty = self.deltaInertialTotalSamples(self.totalCycleLoadPhases,self.inertialCycleLoadPhases)
            
            totalLoadsIndices = np.arange(cycleIndices1[0],cycleIndices1[-1]+deltaSampleQty)
            inertialLoadsIndices = np.arange(cycleIndices2[0],cycleIndices2[-1])
            
            #This will ensure that totalLoadsData does not run out of samples when deltaSampleQty is huge
            # if deltaSampleQty > len(totalLoadsData[:-1,1]) - cycleIndices1[-1]:
                
            #     self.totalCycleLoadPhases = np.delete(self.totalCycleLoadPhases,-1)
            #     totalCycleLoads = np.delete(totalCycleLoads,-1)
            #     self.inertialCycleLoadPhases = np.delete(self.inertialCycleLoadPhases,-1)
            #     inertialCycleLoads = np.delete(inertialCycleLoads,-1)
                
            #     self.surgeCycleSize, deltaSampleQty = self.deltaInertialTotalSamples(self.totalCycleLoadPhases,self.inertialCycleLoadPhases)
                
            #     totalLoadsIndices = np.arange(cycleIndices1[0],cycleIndices1[-2]+deltaSampleQty)
            #     inertialLoadsIndices = np.arange(cycleIndices2[0],cycleIndices2[-2])
            
            totalLoadsDataTrim = totalLoadsData[cycleIndices1[0]:cycleIndices1[-1]+deltaSampleQty,:]
            inertialLoadsDataTrim = inertialLoadsData[cycleIndices2[0]:cycleIndices2[-1],:]
            
        self.totalLoadPhases = totalLoadsDataTrim[:,0]
        self.inertialLoadPhases = inertialLoadsDataTrim[:,0]
        totalLoads = totalLoadsDataTrim[:,1]
        inertialLoads = inertialLoadsDataTrim[:,1]
        
        self.totalCycleCt = totalCycleLoads/(0.5 * totalRhoInf * self.fstreamVel**2 * np.pi * 10e-2**2)
        self.inertialCycleCt = inertialCycleLoads/(0.5 * totalRhoInf * self.fstreamVel**2 * np.pi * 10e-2**2)
        self.totalCt = totalLoads/(0.5 * totalRhoInf * self.fstreamVel**2 * np.pi * 10e-2**2)
        self.inertialCt = inertialLoads/(0.5 * totalRhoInf * self.fstreamVel**2 * np.pi * 10e-2**2)
        
        self.aeroCt = self.totalCt - self.inertialCt
        
        self.trimSurgeCycles = len(self.inertialCycleCt) - 1 #Reduced by 1 to take care off potential issues from rounding off surge cycle size
        
        self.cyclePhases = self.genCycles(self.inertialLoadPhases)
        self.cycleTotal = self.genCycles(self.totalCt)
        self.cycleInertial = self.genCycles(self.inertialCt)
        self.cycleAero = self.genCycles(self.aeroCt)
        
        self.cycleAvgTotal = cycleAvg(self.cycleTotal)
        self.cycleAvgInertial = cycleAvg(self.cycleInertial)
        self.cycleAvgAero = cycleAvg(self.cycleAero)
        
        freq, amp = freqsSignal(inertialLoadsData[:,0] - np.mean(inertialLoadsData[:,0]), 2e3)
        surgeFreq = freq[np.argmax(np.abs(amp))]
        
        filtTotalLoadsData = np.array(totalLoadsData)
        filtInertialLoadsData = np.array(inertialLoadsData)
        
        filtTotalLoadsData[:,1] = lowPassFilt(totalLoadsData[:,1],2e3,np.ceil(surgeFreq))
        filtInertialLoadsData[:,1] = lowPassFilt(inertialLoadsData[:,1],2e3,np.ceil(surgeFreq))
        
        # self.filtTotal = filtTotalLoadsData[:,1]
        # self.filtInertial = filtInertialLoadsData[:,1]
        
        # filtTotalLoadsIntp = CubicSpline(cumTotalLoadPhases[cycleIndices3[1]:],filtTotalLoadsData[cycleIndices3[1]:,1]) #First cycle of recorded data does not have strictly increasing x values
        # filtInertialLoadsIntp = CubicSpline(cumInertialLoadPhases[cycleIndices4[1]:],filtInertialLoadsData[cycleIndices4[1]:,1])
        
        # self.totalIntp = filtTotalLoadsIntp
        # self.inertialIntp = filtInertialLoadsIntp 
        
        filtTotalLoadsDataTrim = filtTotalLoadsData[totalLoadsIndices,:]
        filtInertialLoadsDataTrim = filtInertialLoadsData[inertialLoadsIndices,:]
        
        # self.userDefinedCyclePhases = np.linspace(0,360,self.surgeCycleSize.astype(int),False)
        
        # cloneUserDefinedPhases = np.tile(self.userDefinedCyclePhases,self.trimSurgeCycles)
        # self.phasesTrim = 360*skipcycles*np.ones((self.surgeCycleSize*self.trimSurgeCycles).astype(int))
        
        # for i,phase in enumerate(self.phasesTrim):
            
        #     self.phasesTrim[i] = cloneUserDefinedPhases[i] + self.phasesTrim[i]
            
        #     if self.surgeCycleSize % i+1 == 0:
                
        #         self.phasesTrim[i:] =  self.phasesTrim[i:] + 360
        
        
        
        # filtTotalLoads = filtTotalLoadsIntp(self.phasesTrim)
        # filtInertialLoads = filtInertialLoadsIntp(self.phasesTrim)
        
        # self.filtTotalCt = filtTotalLoads/(0.5 * totalRhoInf * self.fstreamVel**2 * np.pi * 10e-2**2)
        # self.filtInertialCt = filtInertialLoads/(0.5 * totalRhoInf * self.fstreamVel**2 * np.pi * 10e-2**2)
        # self.filtAeroCt = self.filtTotalCt - self.filtInertialCt
        # self.cycleAvgFiltAeroCt, self.cycleStdFiltAeroCt = cycleAvgStd(self.genCycles(self.filtAeroCt))
        # self.filtAeroCtUnc = 1.96 * self.cycleAvgFiltAeroCt/np.sqrt(self.trimSurgeCycles)
        
        
        # filtTotalLoadsDataTrim = np.array(totalLoadsDataTrim)
        # filtInertialLoadsDataTrim = np.array(inertialLoadsDataTrim)
        
        # filtTotalLoadsDataTrim[:,1] = lowPassFilt(totalLoadsDataTrim[:,1],2e3,np.ceil(surgeFreq))
        # filtInertialLoadsDataTrim[:,1] = lowPassFilt(inertialLoadsDataTrim[:,1],2e3,np.ceil(surgeFreq))
        
        self.filtTotalCt = filtTotalLoadsDataTrim[:,1]/(0.5 * totalRhoInf * self.fstreamVel**2 * np.pi * 10e-2**2)
        self.filtInertialCt = filtInertialLoadsDataTrim[:,1]/(0.5 * totalRhoInf * self.fstreamVel**2 * np.pi * 10e-2**2)
        self.filtAeroCt = self.filtTotalCt - self.filtInertialCt
        
        self.cycleAvgFiltAeroCt, self.cycleStdFiltAeroCt = cycleAvgStd(self.genCycles(self.filtAeroCt))
        self.filtAeroCtUnc = 1.96 * self.cycleAvgFiltAeroCt/np.sqrt(self.trimSurgeCycles)
        self.cycleAvgPhases = np.linspace(0,360,self.surgeCycleSize.astype(int))
        
    def cycleSplit(self, loadsData, skipCycles=10):
        
        cycleInitialIndices = []
        cycleCount = 1
        cumSurgePhases = np.array(loadsData[:,0])
        
        for i in range(1,np.shape(loadsData)[0]):
            
            if loadsData[i-1,0] - loadsData[i,0] > 300:
                
                cycleCount = cycleCount + 1
                cumSurgePhases[i:] = cumSurgePhases[i:] + 360
                
                if cycleCount > skipCycles:
                    
                    cycleInitialIndices.append(i)      
        
        cycleQty = len(cycleInitialIndices) - 1 #Removing last cycle since it is expected to be incomplete
        
        cycleSurgePhases = np.empty(cycleQty, dtype = np.ndarray)
        cycleLoads = np.empty(cycleQty, dtype = np.ndarray)
        
        for i in range(cycleQty):
            
            cycleSurgePhases[i] = loadsData[cycleInitialIndices[i]:cycleInitialIndices[i+1],0]
            cycleLoads[i] = loadsData[cycleInitialIndices[i]:cycleInitialIndices[i+1],1]
            
            
        return cycleInitialIndices, cycleSurgePhases, cycleLoads, cumSurgePhases 
    
    def deltaInertialTotalSamples(self, totalCycleLoadPhases, inertialCycleLoadPhases):
        
        totalCycleSizes = []
        inertialCycleSizes = []
        
        for i in range(len(totalCycleLoadPhases)):
    
            totalCycleSizes.append(len(totalCycleLoadPhases[i]))
            inertialCycleSizes.append(len(inertialCycleLoadPhases[i]))
            
        print("Total samples Qty:", sum(totalCycleSizes))
        print("Inertial samples Qty:", sum(inertialCycleSizes))
            
        return np.round(np.mean(inertialCycleSizes)), sum(inertialCycleSizes) - sum(totalCycleSizes)
        
    def genCycles(self, sigTrim):
        
        cycleSig = np.empty(self.trimSurgeCycles, dtype=np.ndarray)
        
        for cycle in range(self.trimSurgeCycles):
            
            cycleIter = np.arange(cycle * self.surgeCycleSize, (cycle+1) * self.surgeCycleSize).astype(int)
            cycleSig[cycle] = sigTrim[cycleIter]
            
        return cycleSig
        
                    
            
