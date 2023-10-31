import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

paths = glob.glob("*/*-dvdl-*")
average_dvdl = []
lam = []
eq_time = 0 #ps

for path in paths:
    time, buff, acid = np.loadtxt(path, comments=["#","@"], unpack=True)
    lam.append(float(path.partition('/')[0]))
    average_dvdl.append(acid)

np.savetxt("avg_dvdl.dat", list(zip(lam, average_dvdl)), fmt = '%.4e')
