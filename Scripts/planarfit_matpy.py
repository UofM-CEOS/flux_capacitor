import numpy as np
from oct2py import octave
from fluxer.eddycov.flux import (planarfit, rotate_vectors)

octave.addpath("~/Scripts/Projects/CEOS/Flux_Capacitor/Scripts/EddyCalc")

np.random.seed(123)
uvw = np.random.randn(10, 3)

# Compare with EddyCalc implementation
k, b = octave.getPlanarFitCoeffs(uvw[:, np.newaxis, 0],
                                 uvw[:, np.newaxis, 1],
                                 uvw[:, np.newaxis, 2])
kpy, bpy = planarfit(uvw)
print np.allclose(k, kpy)       # True
print np.allclose(b.ravel(), bpy)  # True

# Now compare calculated theta, phi
uvw_rot, theta, phi = octave.rotateWindVector(uvw, "PF", k)
uvw_rot_py, (phi_py, theta_py) = rotate_vectors(uvw, k_vector=kpy)
print np.allclose(uvw_rot, uvw_rot_py)  # True
print np.allclose(theta, theta_py)      # True
print np.allclose(phi, phi_py)          # True

# Rui's version
phi_rui, theta_rui = octave.get_tilt_angles(uvw)
print np.allclose(phi, phi_rui)  # False
print np.allclose(theta, theta_rui)  # True
