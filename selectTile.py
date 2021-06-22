import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt
import pdb
from scipy import signal, optimize
from collections import defaultdict



def geom( coords, x0, y0, theta, A, l, w):
    lD = int(np.sqrt(np.shape(coords)[1]))
    X = coords[0].reshape((lD, lD))
    Y = coords[1].reshape((lD, lD))
    thetaMat = np.arctan2(X-x0,Y-y0)-theta+np.pi/2
    rMat = np.sqrt((X-x0)**2 + (Y-y0)**2)
    zFit=2*A**2*rMat*np.cos(thetaMat)/((l)**2)*np.exp(-rMat**2*(np.cos(thetaMat)**2/l**2 + np.sin(thetaMat)**2/w**2))
    return np.ravel(zFit)
    #return(zFit)

raster = gdal.Open('/home/tyler/Research/pitAndMound/Data/rayProperty/025m.tif')
z = raster.ReadAsArray()
dx = 0.25
grad = np.gradient(z)
slp = np.sqrt(grad[0]**2 + grad[1]**2)/dx
cont='y'

try:
    temp = np.load('rayParams.npy', allow_pickle=True).item()
    temp2 = np.load('rayCouplets.npy', allow_pickle=True).item()
    params= defaultdict(list, temp)
    couplets = defaultdict(list, temp2)
except:
    params=defaultdict(list)
    couplets = defaultdict(list)
    couplets['X']=[]
    couplets['Y']=[]

keys=['X','Y','A','l', 'w', 'sd']
pdb.set_trace()
while cont=='y':
    plt.imshow(slp, vmin=0, vmax = 1.5)
    plt.plot(couplets['X'], couplets['Y'], '+r')
    zoom_ok = False
    while not zoom_ok:
        zoom_ok = plt.waitforbuttonpress()
    pts = plt.ginput(-1, timeout = -1)
    plt.close()

    lD = 50
    tile = z[int(pts[0][1])-lD:int(pts[0][1])+lD, int(pts[0][0])-lD:int(pts[0][0])+lD]
    tile = signal.detrend(tile)

    lY, lX = np.shape(tile)
    x, y = np.arange(-int(lX/2),int(lX/2)), np.arange(-int(lY/2),int(lY/2))
    X, Y = np.meshgrid(x,y)

    plt.pcolor(X, Y, tile)
    plt.grid()
    plt.show()
    lD = int(input('define radius in pixels'))
    plt.close()
    tile = tile[int(lY/2)-lD:int(lY/2)+ lD, int(lX/2)-lD:int(lX/2)+lD]
    tile = signal.detrend(tile)

    lY, lX = np.shape(tile)
    x, y = np.arange(-int(lX/2),int(lX/2)), np.arange(-int(lY/2),int(lY/2))
    X, Y = np.meshgrid(x,y)

    coords = np.array([X.ravel(), Y.ravel()])

    p0 = np.array([0.0, 0.0, np.random.random()*np.pi, 0.5, 2.5, 2.5])
    p1=[]
    try:
        p1, pCov = optimize.curve_fit(geom, coords, tile.ravel(), p0)
        temp = geom(coords, p1[0], p1[1],p1[2],p1[3],p1[4],p1[5])
        temp = temp.reshape(lX, lX)
        sqD = np.sum((temp - tile)**2)/(lX**2)
        p1 = np.append(p1, sqD)
        plt.imshow(tile)
        plt.contour(temp,  cmap = 'hot')
        plt.show()
        keep = input('do you wish to keep?y/n')
        if keep=='y':
            k=0
            for key in keys:
                params[key].append(p1[k])
                k+=1
            print(len(params['X']))
            couplets['X'].append(pts[0][0])
            couplets['Y'].append(pts[0][1])
            np.save('rayParams.npy', params, allow_pickle=True)
            np.save('rayCouplets.npy', couplets, allow_pickle=True)
            print('all good')
        cont = input('do you wish to continue? y/n')
        print(len(couplets['X']))
    except:
         print('no acceptable fit')
         cont = input('do you wish to continue? y/n')

pdb.set_trace()
