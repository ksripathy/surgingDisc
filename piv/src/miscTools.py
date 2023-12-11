import numpy as np
from scipy import signal

def movingAverage(a, n=2):
    
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

#Function for obtaining array of attribute of multiple objects of same class organized as an array
def attributeArray(objectArray, attribute, args=None):
    '''Returns the requested attribute of all the objects of numpy array'''
    
    if args == None:#Requested attribute is a variable
        
        return np.array([getattr(obj, attribute) for obj in objectArray])
    
    else:#Requested attribute is a method
        
        return np.array([getattr(obj, attribute)(*args) for obj in objectArray])
    
#Function to obtain the nearest value in a numpy array
def nearestValueIndex(arr, arg):
    
    deltaArr = np.abs(arr - arg)
    nearestIndex = np.argpartition(deltaArr,1)[0]
    
    return nearestIndex
    
#Function to resize axes bbox in a matplotlib figure
def axesBboxSize(bboxWidth, bboxHeight, axesObj):
    
    leftBorder = axesObj.figure.subplotpars.left
    rightBorder = axesObj.figure.subplotpars.right
    topBorder = axesObj.figure.subplotpars.top
    bottomBorder = axesObj.figure.subplotpars.bottom
    
    figWidth = float(bboxWidth)/(rightBorder - leftBorder)
    figHeight = float(bboxHeight)/(topBorder - bottomBorder)
    
    axesObj.figure.set_size_inches(figWidth, figHeight)
    
#Function for finding time shift between signals
def lagShift(sig1, sig2):
    
    corr = signal.correlate(sig1, sig2)/((len(sig1)-1) * np.std(sig1) * np.std(sig2))
    lags = signal.correlation_lags(sig1.size, sig2.size)
    
    return lags[np.argmax(corr)], max(corr)

#Low pass filter function
def lowPassFilt(sig, sampleFreq, cutOffFreq):
    
    Wn = cutOffFreq/(0.5 * sampleFreq)
    b, a = signal.butter(2, Wn)
    
    return signal.filtfilt(b, a, sig)

#Function to trim the end of a signal with incomplete period
def nonCycleTrim(sig, cycleSize, surgeCycleQty):
    
    return sig[:cycleSize*surgeCycleQty]

#Function to trim the initial and final cycles of a signal
def symCycleTrim(sig, cycleSize, trimCycleQty):

    if trimCycleQty != 0:
        
        return sig[cycleSize*trimCycleQty:-cycleSize*trimCycleQty]
    
    else:
        
        return sig
        
#Function to calculate frequeuncies of a signal
def freqsSignal(sig, samplingFreq):
    
    sampleQty = len(sig)
    signalDuration = sampleQty/samplingFreq
    freqHarmonics = np.arange(sampleQty)
    
    twoSideFreqRange = freqHarmonics/signalDuration
    freqRange = twoSideFreqRange[:sampleQty//2]
    
    normSpectralPower = np.fft.fft(sig)/sampleQty
    normSpectralPower = normSpectralPower[:sampleQty//2]
    
    return freqRange, normSpectralPower

#Cycle averaging
def cycleAvg(cyclesArray):
    
    cycleSize = len(cyclesArray[0])
    
    cycleQty = len(cyclesArray)
        
    res = np.zeros(cycleSize)
    
    for i in range(cycleSize):
        
        for cycle in range(cycleQty):
            
            res[i] = res[i] + cyclesArray[cycle][i]
            
        res[i] = res[i]/cycleQty
        
    return res

#Cycle standard deviation
def cycleAvgStd(cyclesArray):
    
    cycleSize = len(cyclesArray[0])
    
    cycleQty = len(cyclesArray)
        
    resAvg = np.zeros(cycleSize)
    resStd = np.zeros(cycleSize)
    
    for i in range(cycleSize):
        
        temp = np.zeros(cycleQty)
        
        for cycle in range(cycleQty):
            
            temp[cycle] = cyclesArray[cycle][i]
            
        resAvg[i] = np.mean(temp)
        resStd[i] = np.std(temp)
        
    return resAvg, resStd

#Update function for cycle data animation
def animFrameUpdate(frame, *fargs):
    
    plotObj = fargs[0]
    filtToggle = fargs[1]
    
    xData = plotObj.aeroObj.cycleTime[frame]
    
    if filtToggle:
    
        y1Data = plotObj.windObj.cycleFiltForce[frame]
        y2Data = plotObj.noWindObj.cycleSyncFiltForce[frame]
        y3Data = plotObj.aeroObj.cycleFiltForce[frame]
        
        minYLim = min(min(y1Data),min(y2Data))
        maxYLim = max(max(y1Data),max(y2Data))
        
        plotObj.filtAnimFig.axes[0].set_xlim(xData[0], xData[-1])
        plotObj.filtAnimFig.axes[0].set_ylim(minYLim, maxYLim)
        
        plotObj.filtAnimFig.axes[0].lines[0].set_data(xData, y1Data)
        plotObj.filtAnimFig.axes[0].lines[1].set_data(xData, y2Data)
        plotObj.filtAnimFig.axes[0].lines[2].set_data(xData, y3Data)
        
        plotObj.filtAnimFig.axes[0].lines[3].set_xdata(xData)
        plotObj.filtAnimFig.axes[0].lines[4].set_xdata(xData)
        plotObj.filtAnimFig.axes[0].lines[5].set_xdata(xData)
        
        return tuple(plotObj.filtAnimFig.axes[0].lines)
        
    else:
        
        y1Data = plotObj.windObj.cycleForce[frame]
        y2Data = plotObj.noWindObj.cycleSyncForce[frame]
        y3Data = plotObj.aeroObj.cycleForce[frame]
        
        minYLim = min(min(y1Data),min(y2Data),min(y3Data))
        maxYLim = max(max(y1Data),max(y2Data),max(y3Data))
        
        plotObj.rawAnimFig.axes[0].set_xlim(xData[0], xData[-1])
        plotObj.rawAnimFig.axes[0].set_ylim(minYLim, maxYLim)
        
        plotObj.rawAnimFig.axes[0].lines[0].set_data(xData, y1Data)
        plotObj.rawAnimFig.axes[0].lines[1].set_data(xData, y2Data)
        plotObj.rawAnimFig.axes[0].lines[2].set_data(xData, y3Data)
        
        plotObj.rawAnimFig.axes[0].lines[3].set_xdata(xData)
        plotObj.rawAnimFig.axes[0].lines[4].set_xdata(xData)
        plotObj.rawAnimFig.axes[0].lines[5].set_xdata(xData)
        
        return tuple(plotObj.rawAnimFig.axes[0].lines)  

#Update function for vortRing animation
def animVortUpdate(frame, *fargs):
        
        plotObj = fargs[0]
        
        funcReturnList = []
        
        for j in range(plotObj.VMObj.mesh.rotor.bladeQuantity):
            
            bladeObj = getattr(plotObj.VMObj.mesh, f"annulusBlade{j}")
            
            #Retreive coordinate locations to plot
            QChordPos = attributeArray(bladeObj.QChordHistory[frame], 'controlPointLoc')
            rootVortPos = attributeArray(bladeObj.rootVortHistory[frame], 'controlPointLoc')
            tipVortPos = tipVortPos = attributeArray(bladeObj.tipVortHistory[frame], 'controlPointLoc')
            
            plotObj.QChordLines[j].set_xdata(QChordPos[:,0])
            plotObj.QChordLines[j].set_ydata(QChordPos[:,1])
            plotObj.QChordLines[j].set_3d_properties(QChordPos[:,2])
            funcReturnList.append(plotObj.QChordLines[j])
            
            plotObj.rootVortLines[j].set_xdata(rootVortPos[:,0])
            plotObj.rootVortLines[j].set_ydata(rootVortPos[:,1])
            plotObj.rootVortLines[j].set_3d_properties(rootVortPos[:,2])
            funcReturnList.append(plotObj.rootVortLines[j])
            
            plotObj.tipVortLines[j].set_xdata(tipVortPos[:,0])
            plotObj.tipVortLines[j].set_ydata(tipVortPos[:,1])
            plotObj.tipVortLines[j].set_3d_properties(tipVortPos[:,2])
            funcReturnList.append(plotObj.tipVortLines[j])
            
        return tuple(funcReturnList)

            
            
            
            
            
    
        
