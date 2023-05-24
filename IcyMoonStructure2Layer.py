'''
Compute structure of a 2-layer icy moon
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
R = 1565e3 # total radius [m]
rhow_arr = np.linspace(950, 1050, 101) # density of water [kg m^-3]
rhoc_arr = np.linspace(3000, 3800, 551) # density of rock-metal core [kg m^-3]
solution = [0, 0, 0, 0]
for rhow in rhow_arr:
    for rhoc in rhoc_arr:
        Rc = ((3*M/4/pi - rhow*R**3)/(rhoc-rhow))**(1/3)
        moi = (8*pi/15) * (rhow*(R**5-Rc**5)+rhoc*Rc**5) /(M*R**2)
        if (abs(MoI-moi) < abs(MoI-solution[3])):
            solution = [rhow, rhoc, Rc, moi]

print("Water-shell density = ", solution[0], " kg m^-3")
print("Rock-metal density = ", solution[1], " kg m^-3")
print("Rock-metal radius = ", solution[2]/1e3, " km")
print("Water-shell thickness = ", (R - Rc)/1e3, " km")
print("Output MoI = ", solution[3])
        
        


