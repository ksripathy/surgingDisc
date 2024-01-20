import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy import signal

data = np.loadtxt("case0DispData.csv", delimiter=",", skiprows=1)
time = data[31:,0]*1e-3
disp = data[31:,1]*1e-3

timeRe = (time[-1] - time[0]) * np.linspace(0,1,len(disp)*2)
dispRe = signal.resample(disp, len(disp)*2)

surgePeriod = 1/4.297183
ndimTime = time/surgePeriod

def mergeCycles(x, y, period): 
    
    ndimX = x/period
    
    cycleQty = int(ndimX[-1])
    ndimXMerge = ndimX - ndimX.astype(int)
    ndimXMergeSort = np.sort(ndimXMerge)
    sortedIndices = np.argsort(ndimXMerge)
    yMergeSort = y[sortedIndices]
    
    return ndimXMergeSort, yMergeSort

timeMerge, dispMerge = mergeCycles(time, disp, surgePeriod)
timeReMerge, dispReMerge = mergeCycles(timeRe, dispRe, surgePeriod)    

dispTheory = 5e-2*np.sin(2*np.pi*timeMerge) + 0.5*(max(dispMerge) + min(dispMerge))

deltaTime = time[-1] - time[-2]

vel = np.diff(dispMerge)/deltaTime
accl = np.diff(vel)/deltaTime

#timeMerged, dispMerged = mergeCycles(time, disp, surgePeriod)
fig1, ax1 = plt.subplots()
ax1.plot(time, disp, marker='x')

fig2,ax2 = plt.subplots()
ax2.plot(timeMerge, dispMerge)
#ax2.plot(timeReMerge, dispReMerge)
ax2.plot(timeMerge, dispTheory)
#ax2.plot(ndimTimeMergeSort[:-1], vel)
#ax2.plot(ndimTimeMergeSort[:-2], accl)
xIndex = 0

'''for x,y in zip(ndimTimeMergeSort,dispMergeSort):
    
    plt.text(x,y,str(xIndex),color="red",fontsize=12)
    xIndex = xIndex + 1'''

#omega = 2*np.pi*
length = 6
radius = 2
x0 = 0.5 * np.pi

def testFunc(x):
    return np.abs(radius*np.cos(x) + np.sqrt(length**2 - (radius**2 * (np.sin(x)**2))) - length)

meanX = minimize(testFunc, x0)

timeMerge = 2*timeMerge

"Theory test-2"
amp = 12.5e-2
dispTheory1 = amp*np.cos(2*np.pi*timeMerge)
acclTheory1 = -amp*(2*np.pi/surgePeriod)**2*np.cos(2*np.pi*timeMerge)

projCrod = 45e-2**2 - (amp*np.sin(2*np.pi*timeMerge))**2
dispTheory2 = amp*np.cos(2*np.pi*timeMerge) + np.sqrt(projCrod) - 45e-2
a1 = amp**2 * (np.cos(2*np.pi*timeMerge)**2 - np.sin(2*np.pi*timeMerge)**2)
a2 = amp**4 * np.sin(2*np.pi*timeMerge)**2 * np.cos(2*np.pi*timeMerge)**2
a3 = np.sqrt(projCrod)
acclTheory2 = -(2*np.pi/surgePeriod)**2 * (amp*np.cos(2*np.pi*timeMerge) + (a1/a3) + (a2/a3**3))

fig3,ax3 = plt.subplots()
ax3.plot(timeMerge, acclTheory1*55e-3/0.1667)
#ax2.plot(timeReMerge, dispReMerge)
ax3.plot(timeMerge, acclTheory2*55e-3/0.1667)
