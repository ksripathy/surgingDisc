import numpy as np
import matplotlib.pyplot as plt
from src.miscTools import freqsSignal
from src.miscTools import lagShift

x = np.linspace(0,6*np.pi,150)
y1 = np.cos(x)
y2 = np.sin(x)

y = 20*y1 + y2
freq = freqsSignal(y,50)[0][1]
freq1 = freqsSignal(y1,50)[0][1]
freq2 = freqsSignal(y2,50)[0][1]
lag = lagShift(y2,y)[0]

fig,ax = plt.subplots()
ax.plot(x,y1)
ax.plot(x,y2)
ax.plot(x,y)