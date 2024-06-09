import numpy as np
from copy import deepcopy

class loadsProcessing:
    
    def __init__(self, totalLoadsObj, inertialLoadsObj):
        
        self.totalLoadsObj = totalLoadsObj
        self.inertialLoadsObj = inertialLoadsObj
        
        #Initatie aero loads obj
        self.aeroLoadsObj = self.initAeroLoadsObj(totalLoadsObj)
        
        self.totalLoadsObj.cycleLoadsSizes = self.totalLoadsObj.cycleSizes()
        self.inertialLoadsObj.cycleLoadsSizes = self.inertialLoadsObj.cycleSizes()
        
        self.cycleQty = min(len(self.totalLoadsObj.cycleSizes()),len(self.inertialLoadsObj.cycleSizes()))
        self.cycleSize = min(np.min(self.totalLoadsObj.cycleSizes()),np.min(self.totalLoadsObj.cycleSizes()))
        
        self.definedPhases = np.linspace(0,360,self.cycleSize,endpoint=False)
        
        self.aeroLoadsObj.loadPhases = np.tile(self.definedPhases,reps=self.cycleQty+totalLoadsObj.skipCycles-1)
        self.aeroLoadsObj.cumLoadPhases = np.array(self.aeroLoadsObj.loadPhases)
        for i in range(len(self.aeroLoadsObj.cumLoadPhases)):
            
            if i % self.cycleSize == 0:
                
                self.aeroLoadsObj.cumLoadPhases[i:] = self.aeroLoadsObj.cumLoadPhases[i:] + 360
        
        
        self.aeroLoadsObj.cycleLoadPhases = np.empty(self.cycleQty, dtype=np.ndarray)
        self.aeroLoadsObj.cycleCumLoadPhases = np.empty(self.cycleQty, dtype=np.ndarray)
        for i in range(self.cycleQty):
            
            self.aeroLoadsObj.cycleLoadPhases[i] = self.definedPhases
            self.aeroLoadsObj.cycleCumLoadPhases[i] =  self.definedPhases + (totalLoadsObj.skipCycles*360) + (i*360)
            
        
            
    def initAeroLoadsObj(self, recLoadsObj):
        
        aeroLoadsObj = deepcopy(recLoadsObj)
        del aeroLoadsObj.cycleStartIndices
        del aeroLoadsObj.cycleLoadPhases
        del aeroLoadsObj.cycleLoads
        del aeroLoadsObj.cumLoadPhases
        aeroLoadsObj.skipCycles = 0
        
        return aeroLoadsObj
    
    def computeAeroLoads(self):
        
        totalLoads = self.totalLoadsObj.cycleIntpLoads
        inertialLoads = self.inertialLoadsObj.cycleIntpLoads
        
        aeroLoads = np.empty(self.cycleQty, dtype=np.ndarray)
        
        for i in range(self.cycleQty):
            
            aeroLoads[i] = totalLoads[i] - inertialLoads[i]
            
        self.aeroLoadsObj.cycleIntpLoads = aeroLoads
        self.aeroLoadsObj.intpLoads = self.totalLoadsObj.intpLoads - self.inertialLoadsObj.intpLoads