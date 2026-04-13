'''
Compute structure of a 3-layer icy moon
    Two knowns: M and MoI
    Two unknowns: core radius, silicate mantle outer radius

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
r3 = R       # ocean-ice shell outer radius [m]
rho3 = 950   # ocean-ice shell density [kg/m^3]
rho2 = 3300  # silicate mantle density [kg/m^3]

## FUNCTIONS
def calc_r1(r2, r3, rho1, rho2, rho3):
    val = (3*M/4/np.pi - (r3**3 - r2**3)*rho3 - r2**3*rho2) / (rho1 - rho2)
    if val < 0:
        return np.nan  # unphysical: no room for a core of this density
    return val**(1/3)

def calc_MOI(r1, r2, r3, rho1, rho2, rho3):
    return 8*np.pi/15/M/R/R * (r1**5*rho1 + (r2**5 - r1**5)*rho2 + (r3**5 - r2**5)*rho3)

def getStructure(r3, rho1, rho2, rho3):
    '''Calculate an interior structure by assuming layer densities.'''
    def moi_residual(r2):
        r1 = calc_r1(r2, r3, rho1, rho2, rho3)
        if np.isnan(r1) or r1 <= 0 or r1 >= r2:
            return np.nan
        return calc_MOI(r1, r2, r3, rho1, rho2, rho3) - MoI

    test_r2_arr = np.linspace(1e3, r3 - 1e3, 2000)
    residuals = np.array([moi_residual(r2) for r2 in test_r2_arr])
    valid = ~np.isnan(residuals)
    sign_changes = np.where(np.diff(np.sign(residuals[valid])))[0]
    if len(sign_changes) == 0:
        return np.nan, np.nan
    idx = sign_changes[0]
    valid_r2 = test_r2_arr[valid]
    r2 = brentq(moi_residual, valid_r2[idx], valid_r2[idx + 1])
    r1 = calc_r1(r2, r3, rho1, rho2, rho3)
    return r1, r2


## GET SOLUTION
n = 100
rho1_arr = np.linspace(4840, 8000, n) # metal core density [kg/m^3]
r1_arr = np.nan * np.zeros(n)
r2_arr = np.nan * np.zeros(n)

for i in range(n):
    r1_arr[i], r2_arr[i] = getStructure(r3, rho1_arr[i], rho2, rho3)

plt.figure()
plt.plot(r1_arr/1e3, rho1_arr)
plt.xlabel('Metal core radius (km)')
plt.ylabel(r'Metal core density ($kg/m^3$)')
