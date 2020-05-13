import matlab.engine
eng = matlab.engine.connect_matlab()
import numpy as np
from numpy import newaxis


import contact_modes as cm
import example as em;
# Parameters
kFrictionH = 0.7
kFrictionE = 0.25
kForceMagnitude = 15

##
## Geometrical Problem definition
##

kW = 0.0435 # object width
kH = 0.0435 # object height

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

R_WH = np.eye(3)
p_WH = np.array(([kW/2, 0, 0]))
p_WH = p_WH[:, newaxis]

kNumContactsE = p_W_e.shape[1];
##
## Geometrical Pre-processing
##

kNumSlidingPlanes = 4
jacs = eng.preProcessing(matlab.double([kFrictionE]),
        matlab.double([kFrictionH]),
        matlab.double([kNumSlidingPlanes]),
        matlab.double(p_W_e.tolist()),
        matlab.double(n_W_e.tolist()),
        matlab.double(p_H_h.tolist()),
        matlab.double(n_H_h.tolist()),
        matlab.double(R_WH.tolist()),
        matlab.double(p_WH.tolist()), nargout=6)
# print('jacs:')
# print(np.array(jacs))

# N_e, T_e, N_h, T_h, eCone_allFix, hCone_allFix
N_e = np.asarray(jacs[0])
T_e = np.asarray(jacs[1])
N_h = np.asarray(jacs[2])
T_h = np.asarray(jacs[3])
eCone_allFix = jacs[4]
hCone_allFix = jacs[5]
b_e = np.zeros((N_e.shape[0], 1))
t_e = np.zeros((T_e.shape[0], 1))
b_h = np.zeros((N_h.shape[0], 1))
t_h = np.zeros((T_h.shape[0], 1))

# get rid of 0 sliding modes
e_modes, cs_lattice, info = cm.enum_sliding_sticking_3d_proj(N_e, b_e, T_e, t_e)
e_cs_modes = np.zeros((len(e_modes), kNumContactsE));
# e_modes_plusminus = []
for i in range(len(e_modes)):
      e_cs_modes[i,:] = e_modes[i][0, 1:kNumSlidingPlanes+1]
      e_modes[i] = e_modes[i][:, kNumSlidingPlanes:]
#     e_modes_sum = np.sum(np.abs(e_modes[i]), 1)
#     valid_row_ids = e_modes_sum >= number_of_none_zeros[i]
#     e_modes_per_cs = np.array(e_modes[i])[valid_row_ids,]
#     e_modes_plusminus.append(e_modes_per_cs)
e_cs_modes = np.array(e_cs_modes)

s_modes = em.modeCleaning(e_cs_modes, e_modes, 4)
# print(s_modes)
num_s = 0
for i in range(len(s_modes)):
  num_s += s_modes[i].shape[0]


print('final s modes:')
print(num_s)