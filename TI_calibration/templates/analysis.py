########### Import modules ###########

import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

########### Define functions ###########

def poly(x: float, As: list) -> list:
    '''
    Evaluates polynomial with coefficients As at points x. Accepts numpy arrays

    Inputs:
    x (float): Point at which to evaluate polynomial
    As (list): Coefficients for polynomial

    Outputs:
    y (float): Value of polynomial at point x
    '''
    rAs = As[::-1]
    x = np.array(x)
    y = np.zeros(len(x))
    for n, a in enumerate(rAs):
        y += a * x ** n
    return y

#Define number of terms for polynomial
avg_dvdl = []
lambdas = []
num_terms = 7

########### Read in data files ###########

paths = glob.glob("*/avg_dvdl.dat")

for path in paths:
    lam, dvdl = np.loadtxt(path, comments = ['#','@'], unpack = True)
    avg_dvdl.append(dvdl)
    lambdas.append(lam)

########### Fit polynomial to dvdl data ###########

popt, pcov = curve_fit(poly,lambdas, avg_dvdl, p0 = np.zeros(num_terms))

#Plot calibration curve
xs = np.linspace(-.1,1.1, 1201)

plt.plot(lambdas, avg_dvdl,'.',label = "TI simulation data", zorder = 2)
plt.plot(xs, poly(xs, popt), label = "Fitted polynomial", zorder = 1)
plt.legend()
plt.set_xlabel("$\lambda$")
plt.set_ylabel("$dV/d\lambda$ /kJmol$^{-1}$")

plt.savefig("calibration_curve.pdf")
plt.savefig("calibration_curve.png")

np.savetxt("raw_coefficients.dat", popt, newline = " ", fmt = "%.4e")
