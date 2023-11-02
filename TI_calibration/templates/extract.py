########### Import modules ###########

import numpy as np
import glob

########### Extract and average raw CpHMD data ###########

paths = glob.glob("*/*-dvdl-*")
average_dvdl = []
lam = []
eq_time = 1_000 #ps

for path in paths:
    time, buff, acid = np.loadtxt(path, comments=["#","@"], unpack=True)
    lam.append(float(path.partition('/')[0]))
    eq_steps = np.where(time == eq_time)[0][0]
    average_dvdl.append(acid[eq_steps:])

np.savetxt("avg_dvdl.dat", list(zip(lam, average_dvdl)), fmt = '%.4e')
