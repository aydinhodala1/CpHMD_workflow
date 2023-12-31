import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def hill(pHs, n, pka):
	return 1 / (10 **(n* (pka - pHs))+1)

paths = glob.glob("*/*-coord-*")
average_lam = []
pH = []
eq_time = 5_000 #ps

for path in paths:
    time,buff,*acids = np.loadtxt(path, comments=["#","@"], unpack=True)
    rounded = []
    eq_steps = np.where(time == eq_time)[0][0]
    for acid in acids:
        for value in acid[eq_steps:]:
            if value > 0.9:
                rounded.append(1)
            elif value < 0.1:
                rounded.append(0)
    pH.append(float(path.partition('/')[0]))
    average_lam.append(np.mean(rounded))

np.savetxt("avg_lam.dat", list(zip(pH, average_lam)), fmt = '%.4e')

popt, pcov = curve_fit(hill, pH, average_lam)

xs = np.linspace(min(pH),max(pH),5000)

plt.plot(pH, average_lam, '.',label = "Simulation data",zorder= 3)
plt.plot(xs, hill(xs, popt[0], popt[1]), label = "Fitted curve", zorder = 2)
plt.legend()
plt.savefig("titration_curve.pdf")
plt.savefig("titration_curve.png")
