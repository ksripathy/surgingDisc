import os
import sys
import pickle

'''File to generate PIV fields data and interpolator from tec file'''

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

from src.ioUtils import caseData

# p45Case0 = caseData(dataDir,"45","00")
# with open(dataDir + f"/pickledData/{p45Case0.uid}" + ".pickle", "wb") as handle:
#     pickle.dump(p45Case0.pickleContainer, handle, protocol=pickle.HIGHEST_PROTOCOL)
# del p45Case0
    
discPors = ["45","70"]
# cases = ["00","01","02","03","04","05","06","07"]
cases = ["06","07"]

for por in discPors:
    
    for Case in cases:
        
        caseObj = caseData(dataDir,por,Case)
        # caseObj = pickleData(dataDir+"/pickledRawDataIntpr",por,Case)
        with open(dataDir + f"/pickledSajad/{caseObj.uid}" + ".pickle", "wb") as handle:
            pickle.dump(caseObj.pickleContainer, handle, protocol=pickle.HIGHEST_PROTOCOL)
        del caseObj