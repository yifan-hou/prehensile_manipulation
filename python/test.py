import matlab.engine
eng = matlab.engine.connect_matlab()
import numpy as np
from numpy import newaxis


import contact_modes as cm
import example as em

# Parameters
kFrictionH = 0.5
kFrictionE = 0.25
kContactForce = 15
kObjWeight = 10
kCharacteristicLength = 0.15

##
## Geometrical Problem definition
##

kW = 0.0435 # object width
kH = 0.02 # object height

# list of contact points and contact normals
p_W_e = np.array(([kW/2, kW/2, 0],
                  [kW/2, -kW/2, 0],
                  [-kW/2, kW/2, 0],
                  [-kW/2, -kW/2, 0])).T
n_W_e = np.array(([0, 0, 1],
                  [0, 0, 1],
                  [0, 0, 1],
                  [0, 0, 1])).T
p_H_h = np.array(([kW/2, kW/2, 0],
                  [kW/2, -kW/2, 0],
                  [-kW/2, kW/2, 0],
                  [-kW/2, -kW/2, 0])).T
n_H_h = np.array(([0, 0, -1],
                  [0, 0, -1],
                  [0, 0, -1],
                  [0, 0, -1])).T

CP_W_G = np.array([0, 0, kH/2]);
CP_W_G = CP_W_G[:, newaxis]


R_WH = np.eye(3)
p_WH = np.array(([0, 0, kH]))
p_WH = p_WH[:, newaxis]

##
## Optional
##
G = np.array([1., 0., 0., 0., 0., 0., 0, 0, 0, 0, 0, 0]);
G = G[newaxis, :]
b_G = np.array([0.1]);

e_mode_goal = np.array([0]);
h_mode_goal = np.array([0]);


##
## Geometrical Pre-processing
##

kNumSlidingPlanes = 4
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
        matlab.double(CP_W_G.tolist()), nargout=7)
# print('jacs:')
# print(np.array(jacs))

# read outputs
N_e = np.asarray(jacs[0])
T_e = np.asarray(jacs[1])
N_h = np.asarray(jacs[2])
T_h = np.asarray(jacs[3])
eCone_allFix = np.asarray(jacs[4])
hCone_allFix = np.asarray(jacs[5])
F_G = np.asarray(jacs[6])

b_e = np.zeros((N_e.shape[0], 1))
t_e = np.zeros((T_e.shape[0], 1))
b_h = np.zeros((N_h.shape[0], 1))
t_h = np.zeros((T_h.shape[0], 1))

J_e = np.vstack((N_e, T_e))
J_h = np.vstack((N_h, T_h))

e_modes, cs_lattice, info = cm.enum_sliding_sticking_3d_proj(N_e, b_e, T_e, t_e)
# divide into cs modes and sliding modes
kNumContactsE = p_W_e.shape[1];
e_cs_modes = np.zeros((len(e_modes), kNumContactsE));
e_ss_modes = [0]*len(e_modes);
for i in range(len(e_modes)):
      e_cs_modes[i,:] = e_modes[i][0, 0:kNumContactsE]
      e_ss_modes[i] = e_modes[i][:, kNumContactsE:].copy()
e_cs_modes = np.array(e_cs_modes)

kNumContactsH = p_H_h.shape[1];
h_cs_modes = np.zeros((1, kNumContactsH)).astype('int32');
h_ss_modes = [np.zeros((1, kNumContactsH*kNumSlidingPlanes)).astype('int32')];

# cs mode: 1 separation, 0 contact
# ss mode: each element can be 1(positive), -1(negative) or 0(on the plane). If all elements for a contaft are 0, this is a sticking contact
# # debug
# e_cs_modes = np.array([[1, 1, 0, 1]]);
# e_ss_modes = [np.array([[0, 0, 0, 0, 1, 1, 0, 0]])];
# h_cs_modes = np.array([[0, 1, 1, 1]]);
# h_ss_modes = [np.array([[1, 1, 0, 0, 0, 0, 0, 0]])];


# # Test modeCleaning()
# s_modes = em.modeCleaning(h_cs_modes, h_ss_modes, 4)
# print(s_modes)
# num_s = 0
# for i in range(len(s_modes)):
#   num_s += s_modes[i].shape[0]
# print('final s modes:')
# print(num_s)

em.wrenchSpaceAnalysis(J_e, J_h, eCone_allFix, hCone_allFix, F_G,
    kContactForce, kCharacteristicLength, kNumSlidingPlanes,
    e_cs_modes, e_ss_modes, h_cs_modes, h_ss_modes, G, b_G, e_mode_goal, h_mode_goal)