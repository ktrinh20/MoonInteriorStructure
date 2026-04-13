'''
Compute structure of a 2-layer icy moon
    Two knowns: M and MoI
    Two unknowns: core radius, core density

    note: MoI = 0.3547 +/- 0.0024 (Casajus et al., 2021) or
                0.3475 +/- 0.0026 (Anderson et al. 1998)

Code clean-up by Claude Opus 4.6
                
Author: Kevin T. Trinh
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

## CONSTRAINTS
M = 4.8e22   # total mass [kg]
MoI = 0.3547 # moment of inertia factor
R = 1561e3   # total radius [m]

## ASSUMPTIONS
r2 = R # outer radius of rocky core shell [m]

## FUNCTIONS
def calc_rho1(r1, r2, rho2):
    return ((3*M/4/np.pi) - (r2**3 - r1**3)*rho2) / r1**3

def calc_MOI(r1, r2, rho1, rho2):
    return 8*np.pi/15/M/R/R * (rho1*r1**5 + rho2*(r2**5 - r1**5))

def getStructure(rho2):
    '''Calculate core size and density given a water shell density.'''
    def moi_residual(r1):
        rho1 = calc_rho1(r1, r2, rho2)
        if rho1 <= 0 or rho1 <= rho2:
            return np.nan
        return calc_MOI(r1, r2, rho1, rho2) - MoI

    test_r1_arr = np.linspace(1e3, r2 - 1e3, 50)
    residuals = np.array([moi_residual(r1) for r1 in test_r1_arr])
    valid = ~np.isnan(residuals)
    sign_changes = np.where(np.diff(np.sign(residuals[valid])))[0]
    if len(sign_changes) == 0:
        return np.nan, np.nan
    idx = sign_changes[0]
    valid_r1 = test_r1_arr[valid]
    r1 = brentq(moi_residual, valid_r1[idx], valid_r1[idx + 1])
    rho1 = calc_rho1(r1, r2, rho2)
    return r1, rho1


## GET SOLUTION
n = 50
rho2_arr = np.linspace(917, 1000, n) # ocean-ice shell density [kg/m^3]

r1_arr = np.nan * np.zeros(n)
rho1_arr = np.nan * np.zeros(n)

for i in range(n):
    r1_arr[i], rho1_arr[i] = getStructure(rho2_arr[i])

plt.figure()
plt.plot(r1_arr/1e3, rho1_arr)
plt.xlabel('Rocky core radius (km)')
plt.ylabel(r'Rocky core density ($kg/m^3$)')
