import numpy as np
from numpy import newaxis
import wrenchStampingLib as ws

a = np.load('wrenchspaceanalysis_out.npy', allow_pickle=True);
ws.wrenchSpaceAnalysis_2d(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
    a[10],a[11],a[12],a[13],a[14]);
