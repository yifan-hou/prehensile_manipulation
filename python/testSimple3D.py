import matlab.engine
eng = matlab.engine.connect_matlab()
import numpy as np
from numpy import newaxis

import contact_modes as cm
import wrenchStampingLib as ws

# Parameters
kFrictionH = 1.0
kFrictionE = 1.0
kContactForce = 15.0
kObjWeight = 10.0
kCharacteristicLength = 0.15
print_level = 1

##
## Geometrical Problem definition
##
#         +Y
#          |
# --1--------------2-> + X
#          |
#          |
#
#         +Z
#          |
#   1______|______2
#    \     |     /
#     \    |    /
#      \   |   /
#       \  |  /
#        \ | /
#         \|/
#----------0---------> + X

kW = 0.0435 # object width
kH = 0.0435 # object height

# list of contact points and contact normals
p_W_e = np.array(([0, kW/2, 0],
                  [0, -kW/2, 0])).T
# p_W_e = p_W_e[:, newaxis]

n_W_e = np.array(([0, 0, 1],
                  [0, 0, 1])).T
# n_W_e = n_W_e[:, newaxis]

p_H_h = np.array(([-kW/2, 0, 0],
                  [ kW/2, 0, 0])).T
n_H_h = np.array(([0, 0, -1],
                  [0, 0, -1])).T

CP_W_G = np.array([0, 0, kH/2]);
CP_W_G = CP_W_G[:, newaxis]

R_WH = np.eye(3)
p_WH = np.array(([0, 0, kH]))
p_WH = p_WH[:, newaxis]


##
## Geometrical Pre-processing
##

kNumSlidingPlanes = 2
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
eTCone_allFix = np.asarray(jacs[5])
hCone_allFix = np.asarray(jacs[6])
hTCone_allFix = np.asarray(jacs[7])
F_G = np.asarray(jacs[8])

b_e = np.zeros((N_e.shape[0], 1))
t_e = np.zeros((eTCone_allFix.shape[0], 1))
b_h = np.zeros((N_h.shape[0], 1))
t_h = np.zeros((hTCone_allFix.shape[0], 1))

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
      e_ss_modes[i] = e_ss_modes[i].astype('int32')
e_cs_modes = np.array(e_cs_modes)

kNumContactsH = p_H_h.shape[1];
e_cs_modes = e_cs_modes.astype('int32')
h_cs_modes = np.zeros((1, kNumContactsH)).astype('int32')
h_ss_modes = [np.zeros((1, kNumContactsH*kNumSlidingPlanes)).astype('int32')]

# cs mode: 1 separation, 0 contact
# ss mode: each element can be 1(positive), -1(negative) or 0(on the plane). If all elements for a contaft are 0, this is a sticking contact
# v: [v_HO, v_HH]
G = np.array([0., 0., 0., 0., 1., 0., 0, 0, 0, 0, 0, 0]);
G = G[newaxis, :]
b_G = np.array([0.1]);

e_cs_modes_goal = np.array([[0, 0]]).astype('int32');
e_ss_modes_goal = [np.array([[0, 0, 0, 0]]).astype('int32')];
h_cs_modes_goal = np.array([[0, 0]]).astype('int32');
h_ss_modes_goal = [np.array([[0, 0, 0, 0]]).astype('int32')];

# e_cs_modes_goal = np.array([[0, 0]]).astype('int32');
# e_ss_modes_goal = [np.array([[0, 0, 0, 0, 0, 0, 0, 0]]).astype('int32')];
# h_cs_modes_goal = np.array([[0, 0]]).astype('int32');
# h_ss_modes_goal = [np.array([[0, 0, 0, 0, 0, 0, 0, 0]]).astype('int32')];

# save arguments to files
np.save('data/J_e.npy', J_e)
np.save('data/J_h.npy', J_h)
np.save('data/eCone_allFix.npy', eCone_allFix)
np.save('data/hCone_allFix.npy', hCone_allFix)
np.save('data/F_G.npy', F_G)
np.save('data/kContactForce.npy', kContactForce)
np.save('data/kFrictionE.npy', kFrictionE)
np.save('data/kFrictionH.npy', kFrictionH)
np.save('data/kCharacteristicLength.npy', kCharacteristicLength)
np.save('data/kNumSlidingPlanes.npy', kNumSlidingPlanes)
np.save('data/e_cs_modes.npy', e_cs_modes)
np.save('data/h_cs_modes.npy', h_cs_modes)
np.save('data/G.npy', G)
np.save('data/b_G.npy', b_G)
np.save('data/e_cs_modes_goal.npy', e_cs_modes_goal)
np.save('data/h_cs_modes_goal.npy', h_cs_modes_goal)
np.save('data/print_level.npy', print_level)

e_ss_modes_rows = []
for i in range(0,len(e_ss_modes)):
  e_ss_modes_rows.append(e_ss_modes[i].shape[0])
e_ss_modes_rows = np.array(e_ss_modes_rows).astype('int32')
e_ss_modes_rows = e_ss_modes_rows[newaxis, :]
e_ss_modes_data = np.concatenate(e_ss_modes, axis=0)
np.save('data/e_ss_modes_data.npy', e_ss_modes_data)
np.save('data/e_ss_modes_rows.npy', e_ss_modes_rows)

h_ss_modes_rows = []
for i in range(0,len(h_ss_modes)):
  h_ss_modes_rows.append(h_ss_modes[i].shape[0])
h_ss_modes_rows = np.array(h_ss_modes_rows).astype('int32')
h_ss_modes_rows = h_ss_modes_rows[newaxis, :]
h_ss_modes_data = np.concatenate(h_ss_modes, axis=0)
np.save('data/h_ss_modes_data.npy', h_ss_modes_data)
np.save('data/h_ss_modes_rows.npy', h_ss_modes_rows)

e_ss_modes_goal_rows = []
for i in range(0,len(e_ss_modes_goal)):
  e_ss_modes_goal_rows.append(e_ss_modes_goal[i].shape[0])
e_ss_modes_goal_rows = np.array(e_ss_modes_goal_rows).astype('int32')
e_ss_modes_goal_rows = e_ss_modes_goal_rows[newaxis, :]
e_ss_modes_goal_data = np.concatenate(e_ss_modes_goal, axis=0)
np.save('data/e_ss_modes_goal_data.npy', e_ss_modes_goal_data)
np.save('data/e_ss_modes_goal_rows.npy', e_ss_modes_goal_rows)

h_ss_modes_goal_rows = []
for i in range(0,len(h_ss_modes_goal)):
  h_ss_modes_goal_rows.append(h_ss_modes_goal[i].shape[0])
h_ss_modes_goal_rows = np.array(h_ss_modes_goal_rows).astype('int32')
h_ss_modes_goal_rows = h_ss_modes_goal_rows[newaxis, :]
h_ss_modes_goal_data = np.concatenate(h_ss_modes_goal, axis=0)
np.save('data/h_ss_modes_goal_data.npy', h_ss_modes_goal_data)
np.save('data/h_ss_modes_goal_rows.npy', h_ss_modes_goal_rows)

print('all printed to file.')

ws.wrenchSpaceAnalysis(J_e, J_h, eCone_allFix, hCone_allFix, F_G,
    kContactForce, kFrictionE, kFrictionH, kCharacteristicLength, kNumSlidingPlanes,
    e_cs_modes, e_ss_modes, h_cs_modes, h_ss_modes, G, b_G,
    e_cs_modes_goal, e_ss_modes_goal, h_cs_modes_goal, h_ss_modes_goal)