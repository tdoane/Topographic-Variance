import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

import pdb

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif'})
#rc('text', usetex=True)

def geom(x0, y0, theta, A, l, w):
    ##Function that is fit to topography
    x = np.arange(-25,25,0.25)
    y = np.arange(-25, 25, 0.25)
    X, Y = np.meshgrid(x,y)
    thetaMat = np.arctan2(X-x0,Y-y0)-theta+np.pi/2
    rMat = np.sqrt((X-x0)**2 + (Y-y0)**2)
    zFit=2*A**2*rMat*np.cos(thetaMat)/((l)**2)*np.exp(-rMat**2*(np.cos(thetaMat)**2/l**2 + np.sin(thetaMat)**2/w**2))
    return zFit


temp = np.load('rayParams.npy', allow_pickle=True).item()
temp2 = np.load('rayCouplets.npy', allow_pickle=True).item()
params= defaultdict(list, temp)
couplets = defaultdict(list, temp2)

dx = 0.25
A = np.abs(params['A'])**2
W = np.abs(params['W'])
L = np.abs(params['L'])
sq = np.abs(params['sd'])
slp = np.array(params['Slp'])
x = params['X']
y = params['Y']
pdb.set_trace()
