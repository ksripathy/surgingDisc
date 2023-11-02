import numpy as np
from copy import deepcopy
from scipy.interpolate import interpn
from scipy.stats import linregress
from numpy.ma import masked_array

class pivFrames2D:
    
    def __init__(self, tecData):
        
        self.posX = tecData[:,0]
        self.posY = tecData[:,1]
        self.velX = tecData[:,2]
        self.velY = tecData[:,3]
        self.vel = np.sqrt(self.velX**2 + self.velY**2)
        self.vortZ = tecData[:,9]
        self.stdVelX = tecData[:,-3]
        self.stdVelY = tecData[:,-2]
        self.isValid = tecData[:,-1].astype(np.bool_)
        
        #Calculate number of rows and cols in FoV
        
        for i in range(len(self.posY)):
            
            if self.posY[i] != self.posY[i+1]:
                
                self.colQty = i+1
                break
            
        self.rowQty = int(len(self.posY)/self.colQty)
        
        self.gridPosX = self.posX.reshape((self.rowQty, self.colQty))
        self.gridPosY = self.posY.reshape((self.rowQty, self.colQty))
        self.gridVelX = self.velX.reshape((self.rowQty, self.colQty))
        self.gridVelY = self.velY.reshape((self.rowQty, self.colQty))
        self.gridVel = self.vel.reshape((self.rowQty, self.colQty))
        self.gridVortZ = self.vortZ.reshape((self.rowQty, self.colQty))
        self.gridStdVelX = self.stdVelX.reshape((self.rowQty, self.colQty))
        self.gridStdVelY = self.stdVelY.reshape((self.rowQty, self.colQty))
        self.gridIsValid = self.isValid.reshape((self.rowQty, self.colQty))
        
    def frameOffset(self, xOset=0, yOset=0):
        
        frameCopy = deepcopy(self)
    
        frameCopy.gridVelX = np.roll(self.gridVelX, (xOset, yOset), axis=(1,0))
        frameCopy.gridVelY = np.roll(self.gridVelY, (xOset, yOset), axis=(1,0))
        frameCopy.gridVel = np.roll(self.gridVel, (xOset, yOset), axis=(1,0))
        frameCopy.gridVortZ = np.roll(self.gridVortZ, (xOset, yOset), axis=(1,0))
        frameCopy.gridStdVelX = np.roll(self.gridStdVelX, (xOset, yOset), axis=(1,0))
        frameCopy.gridStdVelY = np.roll(self.gridStdVelY, (xOset, yOset), axis=(1,0))
        frameCopy.gridIsValid = np.roll(self.gridIsValid, (xOset, yOset), axis=(1,0))
        
        frameCopy.velX = frameCopy.gridVelX.flatten()
        frameCopy.velY = frameCopy.gridVelY.flatten()
        frameCopy.vel = frameCopy.gridVel.flatten()
        frameCopy.vortZ = frameCopy.gridVortZ.flatten()
        frameCopy.stdVelX = frameCopy.gridStdVelX.flatten()
        frameCopy.stdVelY = frameCopy.gridStdVelY.flatten()
        frameCopy.isValid = frameCopy.gridIsValid.flatten()
        
        return frameCopy
    
    def combineFrame(self, otherFrameObj, overlapSmoothing=False):
        
        frameCopy = deepcopy(self)
        #List of column indices corresponding to overlap region
        overlapPosX = []
        
        #Identify xmin and xmax overlap
        for j in range(self.colQty):
                
            if self.gridIsValid[10,j] and otherFrameObj.gridIsValid[10,j]:
                
                overlapPosX.append(j)
                
        frameCopy.overlapIndicesX = overlapPosX
                            
        minOverlapX = overlapPosX[0]
        maxOverlapX = overlapPosX[-1]
        midOverlapX = int(0.5 * (minOverlapX + maxOverlapX))
        
        for i in range(self.rowQty):
                for j in range(self.colQty):
                            
                    if self.gridIsValid[i,j] and otherFrameObj.gridIsValid[i,j]:
                        
                        #Overlap smoothing by associating half of overlap region to values away from FoV edge
                        if overlapSmoothing:
                        
                            if j < midOverlapX:
                                
                                frameCopy.gridVelX[i,j] = otherFrameObj.gridVelX[i,j]
                                frameCopy.gridVelY[i,j] = otherFrameObj.gridVelY[i,j]
                                frameCopy.gridVel[i,j] = otherFrameObj.gridVel[i,j]
                                frameCopy.gridVortZ[i,j] = otherFrameObj.gridVortZ[i,j]
                                frameCopy.gridStdVelX[i,j] = otherFrameObj.gridStdVelX[i,j]
                                frameCopy.gridStdVelY[i,j] = otherFrameObj.gridStdVelY[i,j]
                                frameCopy.gridIsValid[i,j] = otherFrameObj.gridIsValid[i,j]
                                
                            else:
                                
                                frameCopy.gridVelX[i,j] = self.gridVelX[i,j]
                                frameCopy.gridVelY[i,j] = self.gridVelY[i,j]
                                frameCopy.gridVel[i,j] = self.gridVel[i,j]
                                frameCopy.gridVortZ[i,j] = self.gridVortZ[i,j]
                                frameCopy.gridStdVelX[i,j] = self.gridStdVelX[i,j]
                                frameCopy.gridStdVelY[i,j] = self.gridStdVelY[i,j]
                           
                        else:    
                        
                            frameCopy.gridVelX[i,j] = 0.5 * (self.gridVelX[i,j] + otherFrameObj.gridVelX[i,j])
                            frameCopy.gridVelY[i,j] = 0.5 * (self.gridVelY[i,j] + otherFrameObj.gridVelY[i,j])
                            frameCopy.gridVel[i,j] = 0.5 * (self.gridVel[i,j] + otherFrameObj.gridVel[i,j])
                            frameCopy.gridVortZ[i,j] = 0.5 * (self.gridVortZ[i,j] + otherFrameObj.gridVortZ[i,j])
                            frameCopy.gridStdVelX[i,j] = 0.5 * (self.gridStdVelX[i,j] + otherFrameObj.gridStdVelX[i,j])
                            frameCopy.gridStdVelY[i,j] = 0.5 * (self.gridStdVelY[i,j] + otherFrameObj.gridStdVelY[i,j])
                        
                    elif self.gridIsValid[i,j] and ~otherFrameObj.gridIsValid[i,j]:
                        
                        frameCopy.gridVelX[i,j] = self.gridVelX[i,j]
                        frameCopy.gridVelY[i,j] = self.gridVelY[i,j]
                        frameCopy.gridVel[i,j] = self.gridVel[i,j]
                        frameCopy.gridVortZ[i,j] = self.gridVortZ[i,j]
                        frameCopy.gridStdVelX[i,j] = self.gridStdVelX[i,j]
                        frameCopy.gridStdVelY[i,j] = self.gridStdVelY[i,j]
                        
                    elif ~self.gridIsValid[i,j] and otherFrameObj.gridIsValid[i,j]:
                        
                        frameCopy.gridVelX[i,j] = otherFrameObj.gridVelX[i,j]
                        frameCopy.gridVelY[i,j] = otherFrameObj.gridVelY[i,j]
                        frameCopy.gridVel[i,j] = otherFrameObj.gridVel[i,j]
                        frameCopy.gridVortZ[i,j] = otherFrameObj.gridVortZ[i,j]
                        frameCopy.gridStdVelX[i,j] = otherFrameObj.gridStdVelX[i,j]
                        frameCopy.gridStdVelY[i,j] = otherFrameObj.gridStdVelY[i,j]
                        frameCopy.gridIsValid[i,j] = otherFrameObj.gridIsValid[i,j]
                        
        frameCopy.velX = frameCopy.gridVelX.flatten()
        frameCopy.velY = frameCopy.gridVelY.flatten()
        frameCopy.vel = frameCopy.gridVel.flatten()
        frameCopy.vortZ = frameCopy.gridVortZ.flatten()
        frameCopy.stdVelX = frameCopy.gridStdVelX.flatten()
        frameCopy.stdVelY = frameCopy.gridStdVelY.flatten()
        frameCopy.isValid = frameCopy.gridIsValid.flatten()
        
        return frameCopy
    
    def intpVel(self, x, y):
        
        intpArgX = self.gridPosX[0,:]
        intpArgY = self.gridPosY[:,0]
        
        intpVelX = interpn((intpArgY, intpArgX), self.gridVelX, [y,x], method="nearest")
        intpVelY = interpn((intpArgY, intpArgX), self.gridVelY, [y,x], method="nearest")
        
        return intpVelX[0], intpVelY[0]
    
    def velDistrAxial(self, spanLoc, sampleQty, minX = 0, maxX = 0, trimFoV=0):
        
        if minX == 0 and maxX == 0:
        
            minX = min(self.posX) + trimFoV
            maxX = max(self.posX) - trimFoV
        
        axialLocs = np.linspace(minX, maxX, sampleQty)
        
        velX = np.zeros(sampleQty)
        velY = np.zeros(sampleQty)
        
        for i in range(sampleQty):
        
            velX[i], velY[i] = self.intpVel(axialLocs[i], spanLoc)
            
        return axialLocs, velX, velY
    
    def velDistrSpan(self, axialLoc, sampleQty, minY = 0, maxY = 0, trimFoV=0):
        
        if minY == 0 and maxY == 0:
        
            minY = min(self.posY) + trimFoV
            maxY = max(self.posY) - trimFoV
        
        spanLocs = np.linspace(minY, maxY, sampleQty)
        
        velX = np.zeros(sampleQty)
        velY = np.zeros(sampleQty)
        
        for i in range(sampleQty):
        
            velX[i], velY[i] = self.intpVel(axialLoc, spanLocs[i])
            
        return spanLocs, velX, velY
    
    def velDisc(self, spanLoc, discLoc=0):
        
        axialLocs, velXDistr, velYDistr = self.velDistrAxial(spanLoc, 1500, trimFoV=10)
        
        #Find index corresponding to first instance of zero
        for i in range(len(axialLocs)):
            
            if velXDistr[i] == 0:
                firstZeroIndex = i
                
                break
            
        #Find index corresponding to last instance of zero
        revVelXDistr = np.flip(velXDistr)
        
        for i in range(len(axialLocs)):
            
            if revVelXDistr[i] == 0:
                lastZeroIndex = len(axialLocs) - 1 - i
                
                break
            
        #print("First zero index", firstZeroIndex)
        #print("Last zero index", lastZeroIndex)
            
        #Shrink distribution to 10 mm upstream and downstream of last and first zero respectively
        #minAxialLocShrink = axialLocs[firstZeroIndex] - (axialLocs[lastZeroIndex] - axialLocs[firstZeroIndex])
        #maxAxialLocShrink = axialLocs[lastZeroIndex] + (axialLocs[lastZeroIndex] - axialLocs[firstZeroIndex])
        
        #axialLocsShrink, velXDistrShrink, velYDistrShrink = self.velDistrAxial(spanLoc, int(6 * (axialLocs[lastZeroIndex] - axialLocs[firstZeroIndex])), minX = minAxialLocShrink, maxX = maxAxialLocShrink)
        
        #Shrink distribution to 50 mm upstream and downstream of last and first zero respectively
        minAxialLocShrink = axialLocs[firstZeroIndex] - 50
        maxAxialLocShrink = axialLocs[lastZeroIndex] + 50
        
        axialLocsShrink, velXDistrShrink, velYDistrShrink = self.velDistrAxial(spanLoc, int(2 * (maxAxialLocShrink - minAxialLocShrink)), minX = minAxialLocShrink, maxX = maxAxialLocShrink)
        
        #Linear fit
        linFitVelXDistrShrink = linregress(axialLocsShrink, velXDistrShrink)
        intpVelXDistrShrink = linFitVelXDistrShrink.intercept + linFitVelXDistrShrink.slope * axialLocsShrink
        
        #Isolate upstream part from above distribution
        revVelXDistrShrink = np.flip(velXDistrShrink)
        
        for i in range(len(axialLocsShrink)):
            
            if revVelXDistrShrink[i] == 0:
                lastZeroIndexShrink = len(axialLocsShrink) - 1 - i
                
                break
            
        #Linear Fit using upstream velocity alone
        linFitUpVelXDistrShrink = linregress(axialLocsShrink[lastZeroIndexShrink+1:], velXDistrShrink[lastZeroIndexShrink+1:])
        intpUpVelXDistrShrink = linFitUpVelXDistrShrink.intercept + linFitUpVelXDistrShrink.slope * axialLocsShrink
        
        #figDiscVel, axDiscVel = plt.subplots()
        #axDiscVel.plot(axialLocsShrink, velXDistrShrink)
        #axDiscVel.plot(axialLocsShrink, intpVelXDistrShrink)
        #axDiscVel.plot(axialLocsShrink, intpUpVelXDistrShrink)       
        
        #If axial location of disc is not specified use the below location
        if discLoc == 0:
            discLoc = axialLocs[lastZeroIndex]
        
        return discLoc, linFitUpVelXDistrShrink.intercept + linFitUpVelXDistrShrink.slope * discLoc
        #return np.interp(axialLocs[lastZeroIndex]-10, axialLocsShrink, velXDistrShrink)
        
    def recFstreamVelc(self):
        
        spanLocs, velX, velY = self.velDistrSpan(300, 1000, trimFoV=5)
        
        return np.mean(velX)
    
    def discAxialInd(self):
        
        discLoc60, discVel60 = self.velDisc(60.0)
        discLoc72, discVel72 = self.velDisc(72.0, discLoc60)
        discLoc84, discVel84 = self.velDisc(84.0)
        
        
        indLoc60 = 1 - discVel60 / self.recFstreamVelc()
        indLoc72 = 1 - discVel72 / self.recFstreamVelc()
        indLoc84 = 1 - discVel84 / self.recFstreamVelc()
        
        return indLoc60, indLoc72, indLoc84
    
    def discAxialIndv2(self, spanLoc, discLoc):
        
        return 1 - self.velDisc(spanLoc, discLoc)[1]/self.recFstreamVelc()