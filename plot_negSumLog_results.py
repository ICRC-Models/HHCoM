#/usr/bin/python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

filename = '/gscratch/csde/carajb/HHCoM/Params/negSumLogL_calib_07June19.dat'
negSumLogData = np.loadtxt(filename)
fig = plt.figure()

plt.plot(negSumLogData)
plt.xlabel('Run #')
plt.ylabel('negSumLog')
plt.show()
