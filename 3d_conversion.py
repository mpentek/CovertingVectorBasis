# with the kind help from Thomas Oberbichler

import numpy as np

# co-variant basis
g1 = np.array([2.2, 0.0, -2.5])
g2 = np.array([0.0, 1.95, 0.45])
g3 = np.array([-4.5, 0.23, 5.2])

'''
g11 = g1 * g1
g12 = g1 * g2
g13 = g1 * g3
g22 = g2 * g2
g23 = g2 * g3
g33 = g3 * g3

det = 2 * g12 * g13 * g23 - g11 * g23 * g23 - g12 * g12 * g33 + g11 * g22 * g33 - g13 * g13 * g22

G11 = (g22 * g33 - g23 * g23) / det
G12 = (g13 * g23 - g12 * g33) / det
G13 = (g12 * g23 - g13 * g22) / det
G22 = (g11 * g33 - g13 * g13) / det
G23 = (g12 * g13 - g11 * g23) / det
G33 = (g11 * g22 - g12 * g12) / det

G1 = G11 * g1 + G12 * g2 + G13 * g3
G2 = G12 * g1 + G22 * g2 + G23 * g3
G3 = G13 * g1 + G23 * g2 + G33 * g3
'''

# contra-variant basis
G1 = np.cross(g2, g3)/np.dot(g1,np.cross(g2,g3))
G2 = np.cross(g3, g1)/np.dot(g2,np.cross(g3,g1))
G3 = np.cross(g1, g2)/np.dot(g3,np.cross(g1,g2))

# check prints
print(np.dot(g1, G1))
print(np.dot(g2, G2))
print(np.dot(g3, G3))

# in cartesian basis
Pc = np.array([1.25, 1.4, 3.56])

# co-variant component defined by contra-variant basis
P1 = np.dot(Pc, G1)
P2 = np.dot(Pc, G2)
P3 = np.dot(Pc, G3)

# contra-variant component defined by co-variant basis
p1 = np.dot(Pc, g1)
p2 = np.dot(Pc, g2)
p3 = np.dot(Pc, g3)

# resultant 1
Pres1 = p1*G1 + p2*G2 + p3*G3

# resultant 2
# this is needed: co-variant basis with the contra-variant components
Pres2 = P1*g1 + P2*g2 + P3*g3
print(Pres2)

np.testing.assert_allclose(Pres1, Pres2)
np.testing.assert_allclose(Pres2, Pc)
