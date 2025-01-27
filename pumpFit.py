

from cProfile import label
import scipy.optimize as spo

import numpy as np
import matplotlib.pyplot as plt

def pumpCurve(Q,a,b,n):
    """_summary_

    Args:
        Q (_type_): flow rate (m3/s)
        a (_type_): Resistance variable (pressure/flowrate)
        b (_type_): Dead Head (pressure at zero flowrate)
        n (_type_): Flowrate Exponent (q**n)
    """
    
    return b + a * Q**n
    # return b + a * Q **2
    
data = [(0,50), (20/3600,47), (80/3600,27), (110/3600,15), (114/3600,10)]#, (115/3600,0)]
# data = [(0,50), (20,47), (80,27), (110,15), (114,10), (115,0)]

xdata = np.asarray(data,dtype=np.float64)[:,0]
ydata = np.asarray(data,dtype=np.float64)[:,1]

p0 = [-1, ydata[0], 2]

popt, pconv = spo.curve_fit(pumpCurve,xdata,ydata, p0)
print(f"a: {popt[0]}\nb: {popt[1]}\nn: {popt[2]}")
xmodel = np.linspace(xdata[0], xdata[-1])
ymodel = pumpCurve(xmodel,*popt)
# ymodel = pumpCurve(xmodel, -popt[0], -popt[1], popt[2])

fig, ax = plt.subplots()

ax.plot(xdata,ydata,'b-',label='data')
ax.plot(xmodel,ymodel,'r',label='model')
ax.set(xlabel='Q',ylabel='H',title='Head vs. Flowrate')
ax.grid()

plt.show()