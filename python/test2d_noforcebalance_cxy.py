import matlab.engine
eng = matlab.engine.connect_matlab()

import wrenchStampingLib as ws
import numpy as np
from numpy import newaxis



#import contact_modes as cm


# Parameters
kFrictionH = 0.7
kFrictionE = 0.3

kContactForce = 100
kObjWeight = 5

kCharacteristicLength = 0.15

##
## Geometrical Problem definition
##
#
#             -----------
#          h2|          | h1
#        --------------------
#        |        Y         |
#        |        ^         |
#     e2 |        | O       | e1
#    =============|---> X ===========
#
kW = 1 # object width
kH = 0.4 # object height

# list of contact points and contact normals
p_W_e = np.array(([kW/2, -kH/2],
                  [-kW/2, -kH/2],
                  [-kW/2, -kH/2],
                  [-kW/2, kH/2])).T
n_W_e = np.array(([0, 1],
                  [0, 1],
                  [1,0],
                  [1,0])).T
p_H_h = np.array(([[0], [0]]))
n_H_h = -np.array(([[0], [1]]))

CP_W_G = np.array([[0],[0]])

# R_WH = -np.eye(2)
# p_WH = np.array(([[0], [kH/2]]))
R_WH = np.array(([0,1],[-1,0]))
p_WH = np.array(([[kW/2], [kH/4]]))

##
## Geometrical Pre-processing
##

kNumSlidingPlanes = 1 # for 2D planar problems
jacs = eng.preProcessing(matlab.double([kFrictionE]),
        matlab.double([kFrictionH]),
        matlab.double([kNumSlidingPlanes]),
        matlab.double([kObjWeight]),
        matlab.double(p_W_e.tolist()),
        matlab.double(n_W_e.tolist()),
        matlab.double(p_H_h.tolist()),
        matlab.double(n_H_h.tolist()),
        matlab.double(R_WH.tolist()),
        matlab.double(p_WH.tolist()),
        matlab.double(CP_W_G.tolist()), nargout=9)
# print('jacs:')
# print(np.array(jacs))

# read outputs
N_e = np.asarray(jacs[0])
T_e = np.asarray(jacs[1])
N_h = np.asarray(jacs[2])
T_h = np.asarray(jacs[3])
eCone_allFix = np.asarray(jacs[4])
hCone_allFix = np.asarray(jacs[6])
F_G = np.asarray(jacs[8])

J_e = np.vstack((N_e, T_e))
J_h = np.vstack((N_h, T_h))

# [e_modes, h_modes] = sharedGraspModeEnumeration(CP_W_e, CN_W_e, CP_H_h, CN_H_h);
modes = eng.sharedGraspModeEnumeration(matlab.double(p_W_e.tolist()),
        matlab.double(n_W_e.tolist()),
        matlab.double(p_H_h.tolist()),
        matlab.double(n_H_h.tolist()), nargout = 2);
e_modes = np.asarray(modes[0]).astype('int32').T
h_modes = np.asarray(modes[1]).astype('int32').T

##
## Goal
##

# Palm Pivot
G = np.array([0., 1., 0., 0., 0., 0.]);
G = G[newaxis, :]
b_G = np.array([[ 0.0056]]);

# e_mode_goal = np.array([[2, 2]]).astype('int32').T # sf
e_mode_goal = np.array([[0,2,0,1]]).astype('int32').T
h_mode_goal = np.array([[1]]).astype('int32').T # ff

results = ws.wrenchSpaceAnalysis_2d(J_e, J_h, eCone_allFix, hCone_allFix, F_G,
    kContactForce, kFrictionE, kFrictionH, kCharacteristicLength,
    G, b_G, e_modes, h_modes, e_mode_goal, h_mode_goal)

stability_margin = results[0]
print ('stability_margin = ')
print (stability_margin)
print ('results = ')
print (results)
