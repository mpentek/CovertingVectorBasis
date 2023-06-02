import numpy as np

'''
See the repository for more information: https://github.com/mpentek/CovertingVectorBasis/tree/main


Print from MainKratosStatic_RactionOut.py -> Hypar0DegStruct
1
[2.9999999993, -2.9999999987, 4.999999999]
[3](37911.2,-37912,20037.1)

675
[2.9999999987, 2.9999999994, 3.0]
[3](37726.8,37727.3,-19864.1)

676
[-2.9999999993, -2.9999999987, 3.0]
[3](-37726.5,-37726.1,-19863.6)

765
[-2.9999999993, 2.9999999987, 4.999999999]
[3](-37911.6,37910.8,20036.8)
'''

locations = [
    {"corner" : 1,
    "id"    : 676,
    "origin_point": np.array([-3.0, -3.0, 3.0]), 
    "end_points": [np.array([-2.292893, -2.292893, 0.0]), # compression member
                   np.array([-4.931852, -2.482362, 0.0]), # tension member 1
                   np.array([-2.482362, -4.931852, 0.0])], # tension member 2
    "current_components": np.array([-37726.5, -37726.1, -19863.6]),
    "areas" : [1.96e-3, 1.77e-4, 1.77e-4]},
    {"corner" : 2,
    "id"    : 1,
    "origin_point": np.array([3.0, -3.0, 5.0]), 
    "end_points": [np.array([2.292893, -2.292893, 0.0]), # compression member
                   np.array([3, -6, 0.0]), # tension member 1
                   np.array([6, -3, 0.0])], # tension member 2
    "current_components": np.array([37911.2, -37912, 20037.1]),
    "areas" : [1.96e-3, 1.77e-4, 1.77e-4]},
    {"corner" : 3,
    "id"    : 675,
    "origin_point": np.array([3.0, 3.0, 3.0]), 
    "end_points": [np.array([2.292893, 2.292893, 0.0]), # compression member
                   np.array([4.931852, 2.482362, 0.0]), # tension member 1
                   np.array([2.482362, 4.931852, 0.0])], # tension member 2
    "current_components": np.array([37726.8,37727.3,-19864.1]),
    "areas" : [1.96e-3, 1.77e-4, 1.77e-4]},
    {"corner" : 3,
    "id"    : 765,
    "origin_point": np.array([-3.0, 3.0, 5.0]), 
    "end_points": [np.array([-2.292893, 2.292893, 0.0]), # compression member
                   np.array([-3.0, 6.0, 0.0]), # tension member 1
                   np.array([-6.0, 3.0, 0.0])], # tension member 2
    "current_components": np.array([-37911.6, 37910.8, 20036.8]),
    "areas" : [1.96e-3, 1.77e-4, 1.77e-4]}]

# the current is cartesian
current_coord_syst = np.array([[1.0, 0.0, 0.0],
                               [0.0, 1.0, 0.0],
                               [0.0, 0.0, 1.0]])

for loc in locations:
    
    new_coord_syst = np.empty([3, 3])
    for idx, ep in enumerate(loc["end_points"]):
        direction = np.subtract(ep, loc["origin_point"])
        norm = np.linalg.norm(direction)
        if norm > 1e-10:
            direction = direction/norm
        new_coord_syst[idx, :] = direction

    # this is the co-variant basis
    g1 = new_coord_syst[0, :]
    g2 = new_coord_syst[1, :]
    g3 = new_coord_syst[2, :]
    
    # now we determine the contra-variant basis for this co-variant basis
    G1 = np.cross(g2, g3)/np.dot(g1,np.cross(g2,g3))
    G2 = np.cross(g3, g1)/np.dot(g2,np.cross(g3,g1))
    G3 = np.cross(g1, g2)/np.dot(g3,np.cross(g1,g2))
    
    # original components in cartesian coord syst
    orig_force = loc["current_components"]
    
    # the co-variant components are
    [P1, P2, P3] = [np.dot(orig_force, G1), np.dot(orig_force, G2), np.dot(orig_force, G3)]
    # these are the force 
    print('Forces')
    print(P1, P2, P3)
    # these are the prestress values
    print('Prestres')
    print(', '.join([str(p/a) for p,a in zip([P1, P2, P3], loc["areas"])]))
    print()

    # loc["new_components"] = np.matmul(np.matmul(np.linalg.inv(new_coord_syst), current_coord_syst), loc["current_components"])

    # vec_res_current = np.matmul(current_coord_syst, loc["current_components"])
    # vec_res_new = np.matmul(new_coord_syst, loc["new_components"])

    # np.testing.assert_allclose(vec_res_current, vec_res_new)
    
    # print()
    
'''
The results:

Forces
-73889.51109963686 54060.33610210658 54059.74731808234
Prestres
-37698730.15287595, 305425627.6955174, 305422301.23210365

Forces
-101473.49450019744 46336.930351203206 46335.37543069791
Prestres
-51772191.07152931, 261790566.95595032, 261781782.09433848

Forces
-73890.92462700642 54060.78001729857 54061.51599732882
Prestres
-37699451.3403094, 305428135.6909524, 305432293.7702193

Forces
-101472.36246645762 46334.9030942086 46336.45801471389
Prestres
-51771613.50329471, 261779113.5266023, 261787898.38821408

'''
