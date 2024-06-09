import numpy as np
import matplotlib.pyplot as plt

sampleRate = 2000
N = 10

phi = np.linspace(0,20*np.pi,sampleRate*N,endpoint=False)
A1 = 1
f1 = 2
f1Harmonics = np.arange(1,21)*f1
x1 = A1*np.sin(f1*phi)

A2 = 1
f2 = 5
f2Harmonics = np.arange(1,21)*f2
x2 = A2*np.sin(f2*phi)

x = x1 + x2

freqs = np.arange(len(phi))/N

dftx1 = np.fft.fft(x1)
dftx2 = np.fft.fft(x2)
dftx = np.fft.fft(x)

f1Limit = 5
fig1, ax1 = plt.subplots()
ax1.stem(freqs[np.argwhere(freqs<f1Limit)],dftx1[np.argwhere(freqs<f1Limit)])
ax1.vlines(f1Harmonics[np.argwhere(f1Harmonics<f1Limit)], ymin=min(dftx1[np.argwhere(freqs<f1Limit)]), ymax=max(dftx1[np.argwhere(freqs<f1Limit)]), color="k", linestyle="--")

f2Limit = 25
fig2, ax2 = plt.subplots()
ax2.stem(freqs[np.argwhere(freqs<f2Limit)],dftx2[np.argwhere(freqs<f2Limit)], markerfmt="tab:orange", linefmt="tab:orange")
ax2.vlines(f2Harmonics[np.argwhere(f2Harmonics<f2Limit)], ymin=min(dftx2[np.argwhere(freqs<f2Limit)]), ymax=max(dftx2[np.argwhere(freqs<f2Limit)]), color="k", linestyle="--")

fLimit = 20
fig3, ax3 = plt.subplots()
ax3.stem(freqs[np.argwhere(freqs<fLimit)],dftx[np.argwhere(freqs<fLimit)], markerfmt="tab:green", linefmt="tab:green")
ax3.vlines(f1Harmonics[np.argwhere(f1Harmonics<fLimit)], ymin=min(dftx[np.argwhere(freqs<fLimit)]), ymax=max(dftx[np.argwhere(freqs<fLimit)]), color="k", linestyle="dotted")
ax3.vlines(f2Harmonics[np.argwhere(f2Harmonics<fLimit)], ymin=min(dftx[np.argwhere(freqs<fLimit)]), ymax=max(dftx[np.argwhere(freqs<fLimit)]), color="k", linestyle="--")

dftx1Filt = np.zeros(len(dftx1),dtype=complex)
dftx1Filt[20] = dftx1[20]
dftx1Filt[-20] = dftx1[-20]

dftx2Filt = np.zeros(len(dftx2),dtype=complex)
dftx2Filt[50] = dftx2[50]
dftx2Filt[-50] = dftx2[-50]

dftxFilt = np.zeros(len(dftx),dtype=complex)
dftxFilt[20] = dftx[20]
dftxFilt[-20] = dftx[-20]
dftxFilt[50] = dftx[50]
dftxFilt[-50] = dftx[-50]

idftx1 = np.fft.ifft(dftx1Filt)
idftx2 = np.fft.ifft(dftx2Filt)
idftx = np.fft.ifft(dftxFilt)

cycleNo = 0
phiPlot = phi[:sampleRate]
cycleSlice = np.arange(cycleNo*sampleRate,(cycleNo+1)*sampleRate)

fig4, ax4 = plt.subplots()
ax4.plot(phiPlot, x1[cycleSlice])
ax4.plot(phiPlot,idftx1[cycleSlice],"k--")

fig5,ax5 = plt.subplots()
ax5.plot(phiPlot, x2[cycleSlice], "tab:orange")
ax5.plot(phiPlot,idftx2[cycleSlice], "k--")

fig6, ax6 = plt.subplots()
ax6.plot(phiPlot, x[cycleSlice], "tab:green")
ax6.plot(phiPlot, idftx[cycleSlice], "k--")




