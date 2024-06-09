#%%
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pickle


#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

from src.ioUtils import caseData

p45Case6 = caseData(dataDir,"45","06")
p45Case6Container = p45Case6.packFilledData()

#%%

with open(dataDir + f"/pickledDataRev/{p45Case6.uid}" + ".pickle", "wb") as handle:
    
    pickle.dump(p45Case6Container, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
# with open(dataDir + f"/pickledData/{p45Case6.uid}" + ".pickle", "rb") as handle:
    
#     obj = pickle.load(handle)

#%%
fig, ax = plt.subplots()
ax.plot(p45Case6.yLine,p45Case6.indVxcIntpArr)
ax.plot(p45Case6.yLine,p45Case6.indVxcRawArr)

fig2, ax2 = plt.subplots()
ax2.plot(p45Case6.yLine,p45Case6.indVxuIntpArr)
ax2.plot(p45Case6.yLine,p45Case6.indVxuRawArr)

fig3, ax3 = plt.subplots()
ax3.plot(p45Case6.yLine,p45Case6.indVycIntpArr)
ax3.plot(p45Case6.yLine,p45Case6.indVycRawArr)

fig4, ax4 = plt.subplots()
ax4.plot(p45Case6.yLine,p45Case6.indVyuIntpArr)
ax4.plot(p45Case6.yLine,p45Case6.indVyuRawArr)
# %%
