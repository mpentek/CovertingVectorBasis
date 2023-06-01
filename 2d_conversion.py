# re-creating the example from Thomas Oberbichler
# https://observablehq.com/@oberbichler/covariant-and-contravariant

import numpy as np

# co-variant basis
g1 = np.array([2.2, 0.0])
g2 = np.array([0.0, 1.95])

g11 = np.dot(g1,g1)
g12 = np.dot(g1,g2)
g21 = np.dot(g2,g1)
g22 = np.dot(g2,g2)

gij = np.array([[g11, g12],[g21, g22]])
Gij = np.linalg.inv(gij)

# contra-variant basis
G1 = np.matmul(Gij, g1)
G2 = np.matmul(Gij, g2)

# check prints
print(np.dot(g1, G1))
print(np.dot(g1, G2))
print(np.dot(g2, G1))
print(np.dot(g2, G2))

# in cartesian basis
Pc = np.array([1.25, 1.4])

# co-variant component defined by contra-variant basis
P1 = np.dot(Pc, G1)
P2 = np.dot(Pc, G2)

# contra-variant component defined by co-variant basis
p1 = np.dot(Pc, g1)
p2 = np.dot(Pc, g2)

# resultant 1
Pres1 = p1*G1 + p2*G2

# resultant 2
# this is needed: co-variant basis with the contra-variant components
Pres2 = P1*g1 + P2*g2
print(Pres2)

np.testing.assert_allclose(Pres1, Pres2)
np.testing.assert_allclose(Pres2, Pc)