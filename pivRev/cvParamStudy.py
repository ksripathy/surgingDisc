#%%
import os
import sys
import pickle
import math

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy import signal as sig
from scipy.ndimage import median_filter
from matplotlib import colors

from src.ioUtilsPickle import pickleData

discPor = "45"
caseID = "05"
caseObj=pickleData(dataDir+"/pickledFilledNdimData",discPor,caseID)
lowerCV, upperCV = caseObj.circWzParamStudy()

#%%
phaseID = "00"

omegaFluxLowerCV = caseObj.getCntrlVolArrData(lowerCV,"omegaFlux",plotPhaseID=phaseID)
omegaFluxUpperCV = caseObj.getCntrlVolArrData(upperCV,"omegaFlux",plotPhaseID=phaseID)

if int(caseID) < 6: 
    gammaLowerCV = caseObj.getCntrlVolArrData(lowerCV,"gamma",plotPhaseID=phaseID)
    gammaUpperCV = caseObj.getCntrlVolArrData(upperCV,"gamma",plotPhaseID=phaseID)
    
    DGamma_DTLowerCV = caseObj.getCntrlVolArrData(lowerCV,"DGamma_DT",plotPhaseID=phaseID)
    DGamma_DTUpperCV = caseObj.getCntrlVolArrData(upperCV,"DGamma_DT",plotPhaseID=phaseID)
    
    curlLoadLowerCV = DGamma_DTLowerCV - omegaFluxLowerCV
    curlLoadUpperCV = DGamma_DTUpperCV - omegaFluxUpperCV
    
# %%
fig, ax = plt.subplots()
ax.hist(omegaFluxLowerCV.flatten(),density=True,bins=100)

# %% Vorticity study

from src.miscTools import curl2D

Vx = caseObj.phase00.gridVx * caseObj.phase00.VInf
Vy = caseObj.phase00.gridVy * caseObj.phase00.VInf
x = caseObj.gridX * caseObj.discDia
y = caseObj.gridY * caseObj.discDia

Wz = curl2D(x,y,Vx,Vy)

# %%
