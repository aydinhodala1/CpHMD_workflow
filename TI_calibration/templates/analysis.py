import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def poly(x, As):
    rAs = As[::-1]
    x = np.array(x)
    y = np.zeros(len(x))
    for n, a in enumerate(rAs):
        y += a * x ** n
    return y

avg_dvdl = []
avg_dvdl_err = []
lambdas = []
num_terms = 7

paths = glob.glob("*/avg_dvdl.dat")

for path in paths:
    lam, dvdl = np.loadtxt(path, comments = ['#','@'], unpack = True)
    avg_dvdl.append(dvdl)
    lambdas.append(lam)

popt, pcov = curve_fit(poly,lambdas, avg_dvdl, p0 = np.zeros(num_terms))

xs = np.linspace(-.1,1.1, 1201)

plt.plot(lambdas, avg_dvdl,'.',label = "TI simulation data", zorder = 2)
plt.plot(xs, poly(xs, popt), label = "Fitted polynomial", zorder = 1)
plt.legend()
plt.set_xlabel("$\lambda$")
plt.set_ylabel("$dV/d\lambda$ /kJmol$^{-1}$")

plt.savefig("calibration_curve.pdf")
plt.savefig("calibration_curve.png")

np.savetxt("raw_coefficients.dat", popt, newline = " ", fmt = "%.4e")
