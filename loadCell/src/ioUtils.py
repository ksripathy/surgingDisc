import numpy as np
import glob
from scipy.optimize import minimize

from src.forcesRevised import recordedForces
from src.forcesRevised import aeroForces
from src.miscTools import lagShift
from src.glauerts import glauertsCorrectionSingleDOF
from miscTools import nonCycleTrim
from miscTools import symCycleTrim

class caseData:
    
    def __init__(self, tdmsFilePath):
        
        #Lookup table for case metaData. redFreq, redAmp, freq, amp
        caseDict = np.array([[1.0, 0.250, 2.387324, 5.0e-2],
                            [2.0, 0.125, 4.774648, 2.5e-2], 
                            [2.0, 0.250, 4.774648, 5.0e-2], 
                            [2.0, 0.375, 3.183099, 7.5e-2], 
                            [3.0, 0.250, 4.774648, 5.0e-2],
                            [2.7, 0.375, 4.297183, 7.5e-2],
                            [0.0, 0.000, 0.000000, 0.0e-2],
                            [0.0, 0.000, 0.000000, 0.0e-2]])

        p45Mass = 55e-3
        p70Mass = 37e-3
        
        p45Ct = np.array([1.12, 1.02]) #[case6, case7]
        p70Ct = np.array([0.598, 0.563])

        p45FstreamVelc = np.array([2.988726, 3.057711, 2.998582, 2.050731, 2.068581, 2.031872, 3.092307, 2.052653]) #[case0, case1, ..., case7]
        p45FstreamRho = np.array([1.188489, 1.186424, 1.185136, 1.186961, 1.186961, 1.201072, 1.194394, 1.194716])

        p70FstreamVelc = np.array([3.045224, 3.030344, 3.035440, 2.025460, 2.049157, 2.007377, 3.096507, 2.045414])
        p70FstreamRho = np.array([1.196899, 1.198764, 1.198306, 1.194673, 1.192900, 1.189482, 1.194618, 1.195081])

        pivPhaseTimes = np.arange(0.05,1,10)
        
        #Generate case metadata
        self.case = tdmsFilePath.split("/")[-1].split(".")[0]
        self.discPor = float(self.case[1:3])/100
        self.caseNo = int(self.case[7])
        self.redFreq = caseDict[self.caseNo, 0]
        self.redAmp = caseDict[self.caseNo, 1]
        self.vmax = self.redFreq * self.redAmp
        self.approxFreq = caseDict[self.caseNo, 2]
        self.amp = caseDict[self.caseNo, 3]
        self.discDia = 20e-2
        
        if self.discPor == 0.45:
            
            self.discMass = p45Mass
            self.fstreamVel = p45FstreamVelc[self.caseNo]
            self.fstreamRho = p45FstreamRho[self.caseNo]
            
            if round(p45FstreamVelc[self.caseNo]) == 3:
                
                self.staticCt = p45Ct[0]
                
            elif round(p45FstreamVelc[self.caseNo]) == 2:
                
                self.staticCt = p45Ct[1]
            
        elif self.discPor == 0.70:
            
            self.discMass = p70Mass
            self.fstreamVel = p70FstreamVelc[self.caseNo]
            self.fstreamRho = p70FstreamRho[self.caseNo]
            
            if round(p70FstreamVelc[self.caseNo]) == 3:
                
                self.staticCt = p70Ct[0]
                
            elif round(p70FstreamVelc[self.caseNo]) == 2:
                
                self.staticCt = p70Ct[1]
            
        
        self.axialInd = glauertsCorrectionSingleDOF(self.staticCt)
        self.dynPre = 0.5 * self.fstreamRho * self.fstreamVel**2
        self.discArea = 0.25 * np.pi * self.discDia**2
        self.ctNdimTerm = self.dynPre * self.discArea
        
        self.windObj = recordedForces(tdmsFilePath, self.fstreamVel, self.fstreamRho, self.discDia, self.approxFreq)
        #self.freq = self.windObj.computeSurgeFreq(manualMode=True) #accurate estimation of observed surge frequency
        self.surgePeriod = np.mean(self.windObj.surgePeriods[2:-2])
        self.freq = 1/self.surgePeriod
        
        if self.redFreq != 0.0:
            
            #self.aeroMeanLoadTheory = self.staticCt * self.ctNdimTerm * (1 + self.vmax**2)
            
            #print("Old combined load avg:", np.mean(self.windObj.forceSeries))
            
            #Subtract mean load theory from windObj since the combined loading measured in experiment will contain the mean aero load
            #self.windObj.forceSeries = self.windObj.forceSeries - self.aeroMeanLoadTheory
            #self.windObj.filtForce = self.windObj.filtForce - self.aeroMeanLoadTheory
            
            '''print("Theory mean aero load:", self.aeroMeanLoadTheory)
            print("New combined load avg:", np.mean(self.windObj.forceSeries))'''
            
            #Construct theoretical inertial load
            self.theoryInertialLoad()
            
            '''noWindLoadAmp = self.discMass * self.amp * self.freq**2
            aeroLoadAmp = 2 * self.staticCt * self.vmax * self.ctNdimTerm #Kiran's model for disc aero ct in surge
            #aeroLoadAmp = self.staticCt * self.vmax * self.ctNdimTerm #Disc surge ct from clem's results
            
            phaseAngle = 0
            #phaseAngle = np.arctan(aeroLoadAmp/noWindLoadAmp)
            #print("phase angle:", phaseAngle)
            #print("static CT:", self.staticCt)'''
            
            phaseOffset = self.computeOffset()
            print("phase offset:", phaseOffset)
            
            self.aeroObj = aeroForces(self.windObj, self.noWindObj, 2, phaseOffset)
            
            self.windCt = self.windObj.cycleAvgFiltForce/self.ctNdimTerm
            self.noWindCt = self.noWindObj.cycleAvgSyncForce/self.ctNdimTerm
            self.aeroCt = self.aeroObj.cycleAvgFiltForceV2/self.ctNdimTerm
            
        else:
            
            self.aeroCt = np.mean(self.windObj.forceSeries)/self.ctNdimTerm
        
    def theoryInertialLoad(self):
        
        class recForcesContainer(object):
            
            pass
        
        #Dummy object to store theoretical inertial loading
        self.noWindObj = recForcesContainer()
        
        self.noWindObj.timeSeries = self.windObj.timeSeries
    
        angFreq = 2 * np.pi * self.freq
        self.acclSurgeSHM = -self.amp * angFreq**2 * np.sin(angFreq * self.noWindObj.timeSeries)
        
        self.velSurgeSHM = self.amp * angFreq * np.cos(angFreq * self.noWindObj.timeSeries)
        
        self.noWindObj.forceSeries = self.discMass * self.acclSurgeSHM
        self.noWindObj.loadAmp = self.discMass * self.amp * angFreq**2
        self.noWindObj.filtForce = self.noWindObj.forceSeries #Theoretical construction does not require low pass filter
        
    def theoryInertialLoadV2(self):
        
        class recForcesContainer(object):
            
            pass
        
        #Dummy object to store theoretical inertial loading
        self.noWindObj = recForcesContainer()
        
        
        self.noWindObj.timeSeries = self.windObj.timeSeries
        angFreq = 2 * np.pi * self.freq
        
        projCrod = 45e-2**2 - (self.amp*np.sin(angFreq * self.noWindObj.timeSeries))**2
        a1 = self.amp**2 * (np.cos(angFreq * self.noWindObj.timeSeries)**2 - np.sin(angFreq * self.noWindObj.timeSeries)**2)
        a2 = self.amp**4 * np.sin(angFreq * self.noWindObj.timeSeries)**2 * np.cos(angFreq * self.noWindObj.timeSeries)**2
        a3 = np.sqrt(projCrod)
        accl = angFreq**2 * (self.amp*np.cos(angFreq * self.noWindObj.timeSeries) + (a1/a3) + (a2/a3**3))
        
        self.noWindObj.forceSeries = self.discMass * accl
        self.noWindObj.loadAmp = self.discMass * self.amp * angFreq**2
        self.noWindObj.filtForce = self.noWindObj.forceSeries #Theoretical construction does not require low pass filter
        
    def computeOffset(self): #Computes offset of theoretical inertial load wrt theoretical combined load
        
        surgeCycleSize = round(self.surgePeriod/self.windObj.deltaTime)
        surgeCycles = min(int(len(self.windObj.timeSeries) / surgeCycleSize), int(len(self.noWindObj.timeSeries) / surgeCycleSize)) #In case mismatch in length of signals
        
        inerLoadSHM = self.noWindObj.forceSeries
        inerLoadSHMAsymTrim = nonCycleTrim(inerLoadSHM, surgeCycleSize, surgeCycles)
        inerLoadSHMTrim = symCycleTrim(inerLoadSHMAsymTrim, surgeCycleSize, 2)
        
        aero1 = self.staticCt * self.ctNdimTerm
        aeroLoadSHM = aero1 - 2 * self.fstreamRho * self.discArea * self.axialInd * self.velSurgeSHM * (self.axialInd * self.velSurgeSHM + self.fstreamVel * (1 - 2*self.axialInd))
        aeroLoadSHMAsymTrim = nonCycleTrim(aeroLoadSHM, surgeCycleSize, surgeCycles)
        aeroLoadSHMTrim = symCycleTrim(aeroLoadSHMAsymTrim, surgeCycleSize, 2)
        
        combinedLoadSHM = inerLoadSHM + aeroLoadSHM
        combinedLoadSHMTrim = inerLoadSHMTrim + aeroLoadSHMTrim
        
        return lagShift(combinedLoadSHMTrim, inerLoadSHMTrim)[0]
        
        
        
    