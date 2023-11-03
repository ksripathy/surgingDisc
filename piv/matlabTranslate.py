from scipy.io import savemat
import numpy as np

testDic = {"a" : np.arange(5), "b" : np.zeros((5,5)), "c" : {"d" : np.arange(5), "e" :  np.zeros((5,5))}}

savemat("testMlabFile.mat", testDic)
