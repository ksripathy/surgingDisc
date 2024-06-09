import os
import sys
import pickle

'''File to obtain velocity, circulation and plots
from pickled piv fields data generated using genData.py'''

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

from src.ioUtilsPickle import pickleData

# p45Case0 = caseData(dataDir,"45","00")
# with open(dataDir + f"/pickledData/{p45Case0.uid}" + ".pickle", "wb") as handle:
#     pickle.dump(p45Case0.pickleContainer, handle, protocol=pickle.HIGHEST_PROTOCOL)
# del p45Case0
    
discPors = ["45","70"]
cases = ["00","01","02","03","04","05","06","07"]

for por in discPors:
    
    for Case in cases:
        
        print(f"\nProcessing P{por}Case{Case}...")
        
        caseObj = pickleData(dataDir+"/pickledFilledNdimData",por,Case)
        with open(dataDir + f"/pickledPlotNdimDataSegs/{caseObj.uid}" + ".pickle", "wb") as handle:
            pickle.dump(caseObj.pickleContainer, handle, protocol=pickle.HIGHEST_PROTOCOL)
        del caseObj
        
        print(f"Processed P{por}Case{Case}!")
