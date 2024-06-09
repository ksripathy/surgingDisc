import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt

from scipy.io import savemat

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

with open(dataDir + "/pickledSajad/p70Case07.pickle","rb") as handle:
    caseObj = pickle.load(handle)
    savemat(dataDir + "/pickledSajad/p70Case07.mat", caseObj.__dict__)
    
    
