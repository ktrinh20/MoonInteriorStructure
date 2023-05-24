'''
Compute structure of a 3-layer icy moon
    Need a few equations:
    - hydrostatic equilibrium
    - mass conservation
    
    Assume we know layer densities. Want to know ocean-ice shell thickness!
    Two knowns: M and MoI
    Two unknowns: core size, ocean-ice shell thickness
    
    note: MoI = 0.3547 +/- 0.0024 (Casajus et al., 2021) or 
                0.3475 +/- 0.0026 (Anderson et al. 1998)
'''
from math import pi
import numpy as np

M = 4.8e22  # mass of sphere

err = 0.0024     # error in MoI
MoI = 0.3547
R1 = 1565e3 # total radius [m]
rhoa = M / (4*pi/3*R1**3)   # bulk (average) density
rho1 = 1000 # ocean-ice shell density [kg / m^3]
rho2 = 3300     # silicate mantle density
rho3 = 5150     # core density

maxR3 = 1500e3   # max core radius tested [m]
R3 = 1e3
dr3 = 5e2   # increase in tested core size
m_arr = [] # make sure that mass makes sense [kg]
R3_arr = []   # store tested metallic core radius [m]
R2_arr = []   # store silicate mantle radius [m]
h_arr = []  # store ocean-ice shell thickness (R1 - R2_arr)
MoI_arr = []    # store calculated MoI
while R3 <= maxR3:
    # compute ocean-ice shell thickness to satisfy mass constraint
    R3_arr.append(R3)
    m3 = rho3 * 4*pi/3*R3**3    # metallic core mass [kg]
    r2 = ((3*(M-m3)/(4*pi) + rho2*R3**3 - rho1*R1**3)/(rho2-rho1))**(1/3)
    R2_arr.append(r2)
    h_arr.append(R1 - r2)
    
    # check mass
    m_arr.append(4*pi/3 * (rho3*R3**3 + rho2*(r2**3-R3**3) + rho1*(R1**3-r2**3)))
    
    # compute MoI using Eq. 6 in Schubert et al. (2009)
    moi = (2/5)*(rho1/rhoa + ((rho3-rho2)/rhoa)*(R3/R1)**5 + 
                 ((rho2-rho1)/rhoa)*(r2/R1)**5)
    MoI_arr.append(moi)
    R3 += dr3
    
# find structure that leads to lowest error in MoI
MoI_arr = np.asarray(MoI_arr)
diff_arr = np.absolute(MoI_arr - MoI)
index = diff_arr.argmin()
within_error = abs(MoI_arr[index] - MoI) <= err
print("MoI error = ", (MoI_arr[index] - MoI) / MoI * 100, " %")
print("MoI result = ", MoI_arr[index])
print("Within MoI error? ", within_error)
print("Mass error = ", (m_arr[index] - M)/M * 100, " %")
print("Core radius = ", R3_arr[index] / 1e3, " km")
print("Core mass fraction = ", rho3*4*pi/3*R3_arr[index]**3 / M, " %")
print("Ocean-ice shell thickness = ", h_arr[index] / 1e3, " km")
print("H2O mass fraction = ", rho1*4*pi/3*(R1**3-R2_arr[index]**3)/M, " %")
