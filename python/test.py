import matlab.engine
eng = matlab.engine.connect_matlab()
import numpy as np
from numpy import newaxis

import contact_modes as cm
import wrenchStampingLib as ws

# Parameters
kFrictionH = 0.5
kFrictionE = 0.25
kContactForce = 15.0
kObjWeight = 10.0
kCharacteristicLength = 0.15

##
## Geometrical Problem definition
##
# 3       +Y        1
#          |
# -------------------> + X
#          |
# 4        |        2
#
kW = 0.0435 # object width
kH = 0.0435 # object height

# list of contact points and contact normals
p_H_e = np.array(([kW/2, kW/2, -kH],
                  [kW/2, -kW/2, -kH],
                  [-kW/2, kW/2, -kH],
                  [-kW/2, -kW/2, -kH])).T
n_H_e = np.array(([0, 0, 1],
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

CP_H_G = np.array([0, 0, kH/2]);
CP_H_G = CP_H_G[:, newaxis]
z_H = np.array([0, 0, -1]);
z_H = z_H[:, newaxis]

##
## Geometrical Pre-processing
##

kNumSlidingPlanes = 2
jacs = eng.preProcessing(matlab.double([kFrictionE]),
        matlab.double([kFrictionH]),
        matlab.double([kNumSlidingPlanes]),
        matlab.double([kObjWeight]),
        matlab.double(p_H_e.tolist()),
        matlab.double(n_H_e.tolist()),
        matlab.double(p_H_h.tolist()),
        matlab.double(n_H_h.tolist()),
        matlab.double(CP_H_G.tolist()),
        matlab.double(z_H.tolist()),
        nargout=7)
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

J_e = np.vstack((N_e, T_e))
J_h = np.vstack((N_h, T_h))

e_modes, cs_lattice, info = cm.enum_sliding_sticking_3d_proj(N_e, b_e, T_e, t_e)
# divide into cs modes and sliding modes
kNumContactsE = p_H_e.shape[1];
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
# # debug
# e_cs_modes = np.array([[1, 1, 0, 1]]);
# e_ss_modes = [np.array([[0, 0, 0, 0, 1, 1, 0, 0]])];
# h_cs_modes = np.array([[0, 1, 1, 1]]);
# h_ss_modes = [np.array([[1, 1, 0, 0, 0, 0, 0, 0]])];

G = np.array([0., 0., 1., 0., 0., 0., 0, 0, 0, 0, 0, 0]);
G = G[newaxis, :]
b_G = np.array([0.1]);

e_cs_modes_goal = np.array([[0, 0, 1, 1]]).astype('int32');
e_ss_modes_goal = [np.array([[0, 0, 0, 0, 0, 0, 0, 0]]).astype('int32')];
h_cs_modes_goal = np.array([[0, 0, 0, 0]]).astype('int32');
h_ss_modes_goal = [np.array([[0, 0, 0, 0, 0, 0, 0, 0]]).astype('int32')];

# # Test modeCleaning()
# s_modes = em.modeCleaning(h_cs_modes, h_ss_modes, 4)
# print(s_modes)
# num_s = 0
# for i in range(len(s_modes)):
#   num_s += s_modes[i].shape[0]
# print('final s modes:')
# print(num_s)

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