import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def hill(pHs, n, pka):
	return 1 / (10 **(n* (pka - pHs))+1)

paths = glob.glob("../*/avg_lam.dat")
average_lam = []
pH = []

for path in paths:
    ph_temp, avg_lam_temp = np.loadtxt(path, unpack = True)
    pH = np.append(pH, ph_temp)
    average_lam = np.append(average_lam, avg_lam_temp)

popt, pcov = curve_fit(hill, pH, average_lam)

print(popt)

xs = np.linspace(min(pH),max(pH),5000)

plt.plot(pH, average_lam, '.',label = "Simulation data",zorder= 3)
plt.plot(xs, hill(xs, popt[0], popt[1]), label = "Fitted curve", zorder = 2)
plt.legend()
plt.savefig("titration_curve.png")
plt.savefig("titration_curve.pdf")