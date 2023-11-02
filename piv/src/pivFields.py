import numpy as np
from scipy.interpolate import interpn
from scipy.stats import linregress 
import matplotlib.pyplot as plt
from pivFrames import pivFrames2D
from miscTools import nearestValueIndex

class planarPIVField():
    
    def __init__(self, filePath):
    
        with open(filePath, 'r') as tecFile:
            tecFileLines = tecFile.read().split("\n")
            
        framesFirstLineLocs = []
            
        #Find number of frames in the tecFile
        for i,line in enumerate(tecFileLines):
            
            if "ZONE" in line.split():
                framesFirstLineLocs.append(i)
                
        frameQty = len(framesFirstLineLocs)
        
        if frameQty == 1:
            self.frame0 = pivFrames2D(np.loadtxt(filePath, skiprows=4))
            
        elif frameQty > 1:
            
            for i in range(frameQty - 1):
                setattr(self,f"frame{i}",pivFrames2D(np.loadtxt(filePath, skiprows=framesFirstLineLocs[i]+2, max_rows = framesFirstLineLocs[i+1] - 4 - 1)))#Max rows is the number of rows in tecFile corresponding to a frame
                
            #Storing final frame    
            setattr(self,f"frame{frameQty-1}",pivFrames2D(np.loadtxt(filePath, skiprows=framesFirstLineLocs[-1]+2)))
            
        self.frameAxialLocs = self.frame0.gridPosX[0,:]
        self.frameSpanLocs = self.frame0.gridPosY[:,0]
        self.deltaPosX = np.round(self.frame0.gridPosX[0,1] - self.frame0.gridPosX[0,0],3)
        self.deltaPosY = np.round(self.frame0.gridPosY[0,0] - self.frame0.gridPosX[1,0],3)
            
    def combine2Frames(self, frame0OffsetParm=0, overlapSmoothing=False):
        
        if frame0OffsetParm!=0:
            offsetFrame0 = self.frame0.frameOffset(frame0OffsetParm)
            
        else:
            offsetFrame0 = self.frame0
            
        self.combinedFrames = self.frame1.combineFrame(offsetFrame0, overlapSmoothing)
        self.overlapAxialIndices = self.combinedFrames.overlapIndicesX
        
    def maskBoundaryIndex(self, frameNo, spanIndex, calbEdgeArtifactTrim = 0):
        
        #With non-zero calbEdgeArtifactTrim, the mask boundary can be better detected by the code, if for zero value the mask boundary is not detected
        #This is required when frame.gridIsValid = 0 at the left most edge of the frame in special circumstances
        #For P70Case3,4,5 use non-zero calbEdgeArtifactTrim for better mask boundary detection
        
        if frameNo == -1:
            frameObj = getattr(self,"combinedFrames")
            
        else:
            frameObj = getattr(self,f"frame{frameNo}")
        
        startAxialIter = self.overlapAxialIndices[0] + 1 + calbEdgeArtifactTrim
        endAxialIter = len(frameObj.gridPosX[spanIndex, :]) + 1
        
        for axialIndex in range(startAxialIter, endAxialIter):
            
            if bool(frameObj.gridIsValid[spanIndex, axialIndex]) - bool(frameObj.gridIsValid[spanIndex, axialIndex-1]) == 1:
                ret = axialIndex
                break
            
        return ret
    
    #Provide odd number for distribution size to obtain induction at disc centre
    def discVelXDistr(self, frameNo, discCentreLoc, discAxialLoc, intpWinSize=10, distrSize=21):
        
        '''Returns an array of interpolated axial velocity at different spanwise locations.
        Velocity is interpolated from the upstream mask edge corresponding to each spanwise location'''
        
        if frameNo == -1:
            frameObj = getattr(self,"combinedFrames")
            
        else:
            frameObj = getattr(self,f"frame{frameNo}")
            
        spanLocs = np.linspace(discCentreLoc - 100, discCentreLoc + 100, distrSize)
        velXDistr = np.empty(distrSize, dtype=np.ndarray)
        
        for i in range(distrSize):
            
            nearestSpanIndex = nearestValueIndex(self.frameSpanLocs, spanLocs[i])
            startAxialIndex = self.maskBoundaryIndex(frameNo, nearestSpanIndex)
            
            if intpWinSize == 0:
                
                velXDistr[i] = frameObj.intpVel(frameObj.gridPosX[nearestSpanIndex,startAxialIndex], spanLocs[i])[0]
                
            else:
            
                if i == int(0.5 * (distrSize + 1) - 1): #At disc centre
            
                    endAxialIndex = startAxialIndex + intpWinSize
                        
                else: #Non disc centre. Here dditional 10 samples in window for the additional field between disc and mounting screw 
                    
                    endAxialIndex = startAxialIndex + intpWinSize + 10
                    
                intpVelX = np.zeros(endAxialIndex - startAxialIndex + 1)
                
                for j,axialIndex in enumerate(range(startAxialIndex, endAxialIndex+1)):
                    
                    intpVelX[j] = frameObj.intpVel(frameObj.gridPosX[nearestSpanIndex,axialIndex], spanLocs[i])[0]
                    
                #intpObj = linregress(frameObj.gridPosX[nearestSpanIndex,startAxialIndex:endAxialIndex+1],intpVelX)
                #velXDistr[i] = intpObj.intercept + intpObj.slope * discAxialLoc
                
                polyFit5 = np.polyfit(frameObj.gridPosX[nearestSpanIndex,startAxialIndex:endAxialIndex+1], intpVelX, deg=3)
                velXDistr[i] = np.poly1d(polyFit5)(discAxialLoc)
            
        return velXDistr
    
    def discVelXDistrV2(self, frameNo, discCentreLoc, discAxialLoc, intpWinSize=10, distrSize=21):
        
        '''Returns an array of interpolated axial velocity at different spanwise locations.
        Velocity is interpolated from the upstream mask edge corresponding to the disc centre location'''
        
        if frameNo == -1:
            frameObj = getattr(self,"combinedFrames")
            
        else:
            frameObj = getattr(self,f"frame{frameNo}")
            
        spanLocs = np.linspace(discCentreLoc - 100, discCentreLoc + 100, distrSize)
        velXDistr = np.empty(distrSize, dtype=np.ndarray)
        
        nearestSpanIndex = nearestValueIndex(self.frameSpanLocs, discCentreLoc)
        startAxialIndex = self.maskBoundaryIndex(frameNo, nearestSpanIndex, calbEdgeArtifactTrim=5)
        endAxialIndex = startAxialIndex + intpWinSize
        
        for i in range(distrSize):
            
            if intpWinSize == 0:
                
                velXDistr[i] = frameObj.intpVel(frameObj.gridPosX[nearestSpanIndex,startAxialIndex], spanLocs[i])[0]
                
            else:
            
                intpVelX = np.zeros(endAxialIndex - startAxialIndex + 1)
                
                for j,axialIndex in enumerate(range(startAxialIndex, endAxialIndex+1)):
                    
                    intpVelX[j] = frameObj.intpVel(frameObj.gridPosX[nearestSpanIndex,axialIndex], spanLocs[i])[0]
                    
                #intpObj = linregress(frameObj.gridPosX[nearestSpanIndex,startAxialIndex:endAxialIndex+1],intpVelX)
                #velXDistr[i] = intpObj.intercept + intpObj.slope * discAxialLoc
                
                polyFit5 = np.polyfit(frameObj.gridPosX[nearestSpanIndex,startAxialIndex:endAxialIndex+1], intpVelX, deg=3)
                velXDistr[i] = np.poly1d(polyFit5)(discAxialLoc)
            
        return velXDistr
    
    def discAxialIndDistr(self, frameNo, discCentreLoc, discAxialLoc, intpWinSize=10, distrSize=21):

        if frameNo == -1:
            frameObj = getattr(self,"combinedFrames")
            
        else:
            frameObj = getattr(self,f"frame{frameNo}")
            
        fstreamVelc = frameObj.recFstreamVelc()
        
        discVelXDistr = self.discVelXDistr(frameNo, discCentreLoc, discAxialLoc, intpWinSize, distrSize)
        axialIndDistr = 1 - discVelXDistr/fstreamVelc
        
        return axialIndDistr
    
    def discAxialIndDistrV2(self, frameNo, discCentreLoc, discAxialLoc, intpWinSize=10, distrSize=21):

        if frameNo == -1:
            frameObj = getattr(self,"combinedFrames")
            
        else:
            frameObj = getattr(self,f"frame{frameNo}")
            
        fstreamVelc = frameObj.recFstreamVelc()
        
        discVelXDistr = self.discVelXDistrV2(frameNo, discCentreLoc, discAxialLoc, intpWinSize, distrSize)
        axialIndDistr = 1 - discVelXDistr/fstreamVelc
        
        return axialIndDistr
    
    def axialIndSpanDistr(self, frameNo, spanLoc, axialLoc):
        
        if frameNo == -1:
            frameObj = getattr(self,"combinedFrames")
            
        else:
            frameObj = getattr(self,f"frame{frameNo}")
            
        axialVelXDistr = frameObj.velDistrSpan(axialLoc, 21, spanLoc-100, spanLoc+100)[1]
        fstreamVelc = frameObj.recFstreamVelc()
        
        axialIndSpanDistr = 1 - axialVelXDistr/fstreamVelc
        
        return axialIndSpanDistr
              
    #Todo: change the reference for the frame object in remaining methods from self to self.combinedFrame
    def intpVel(self, x, y):
        
        intpArgX = self.gridPosX[0,:]
        intpArgY = self.gridPosY[:,0]
        
        intpVelX = interpn((intpArgY, intpArgX), self.gridVelX, [y,x], method="linear")
        intpVelY = interpn((intpArgY, intpArgX), self.gridVelY, [y,x], method="linear")
        
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
        
        