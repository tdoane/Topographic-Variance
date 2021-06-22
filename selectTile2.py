##This script evaluates the fit for selected pits and allows for the user to decide to keep or remove couplets from the analysis based on the visual fit and sum of squared differences between the idealized geometry and the natural geometry.

import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt
import pdb
from scipy import signal, optimize
from collections import defaultdict
import matplotlib.colors as mcolors

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))
def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]
def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list.

        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.

        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))

    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp

def geom(coords, x0, y0, theta, A, l, w):
    lD = int(np.sqrt(np.shape(coords)[1]))
    dx = 0.25
    X = coords[0].reshape((lD, lD))
    X=X*0.25
    Y = coords[1].reshape((lD, lD))
    Y=Y*0.25
    thetaMat = np.arctan2(X-x0,Y-y0)-theta+np.pi/2
    rMat = np.sqrt((X-x0)**2 + (Y-y0)**2)
    zFit=2*A**2*rMat*np.cos(thetaMat)/((l)**2)*np.exp(-rMat**2*(np.cos(thetaMat)**2/l**2 + np.sin(thetaMat)**2/w**2))
    return np.ravel(zFit)

def lowPass(z, l0):
    #create a low-pass filter that smooths topography using a Gaussian kernel
    from scipy.signal import detrend
    lY, lX = np.shape(z)
    x, y = np.arange(-lX/2, lX/2), np.arange(-lY/2, lY/2)
    X, Y = np.meshgrid(x, y)
    filt = 1/(2*np.pi*l0**2)*np.exp(-(X**2 + Y**2)/(2*l0**2))
    ftFilt = np.fft.fft2(filt)
    #z = detrend(z, axis=0)
    #z = detrend(z, axis=-1)
    ftZ = np.fft.fft2(z)
    ftZNew = ftZ*ftFilt
    zNew = np.fft.ifft2(ftZNew).real
    zNew = np.fft.fftshift(zNew)
    return(zNew)

myColors = np.load('myColors.npy', allow_pickle = 'True').item()
hex_list = myColors['Arizona Sunset']
colors = plt.get_cmap(get_continuous_cmap(hex_list))

raster = gdal.Open('/home/tyler/Research/pitAndMound/Data/rayProperty/025m.tif')
z = raster.ReadAsArray()
dx = 0.25
l0 = 5*2.5*0.3048/dx
lPass = lowPass(z, l0)
hPass = z - lPass
grad = np.gradient(lPass)
slp = np.sqrt(grad[0]**2 + grad[1]**2)/dx

temp2 = np.load('rayCouplets.npy', allow_pickle=True).item()
couplets = defaultdict(list, temp2)
params = defaultdict(list)

keys = ['X', 'Y', 'Theta', 'A', 'L', 'W', 'SD']

X = couplets['X']
Y = couplets['Y']
dt=15

x = np.arange(-dt, dt)
y = np.arange(-dt, dt)
XX, YY = np.meshgrid(x,y)
XX = XX.ravel()
YY = YY.ravel()
coords = np.array([XX, YY])

k=0
while k <len(Y):
    i = int(Y[k])
    j = int(X[k])
    tile = z[i-dt:i+dt, j-dt:j+dt]
    tile = signal.detrend(tile)
    lY, lX = np.shape(tile)
    p0 = np.array([0.0, 0.0, np.random.random()*np.pi, 0.75, 2.5, 2.5])
    p1=[]
    try:
        p1, pCov = optimize.curve_fit(geom, coords, tile.ravel(), p0)
        temp = geom(coords, p1[0], p1[1],p1[2],p1[3],p1[4],p1[5])
        temp = temp.reshape(lX, lX)
        sqD = np.sum((temp - tile)**2)/(lX**2)
        p1 = np.append(p1, sqD)
        fig, ax = plt.subplots()
        im=ax.pcolor(YY.reshape(2*dt, 2*dt), XX.reshape(2*dt, 2*dt), tile, cmap = colors)
        fig.colorbar(im)
        CS =ax.contour(YY.reshape(2*dt, 2*dt), XX.reshape(2*dt, 2*dt), temp,  colors = 'k')
        ax.clabel(CS, inline=1, fontsize=12)
        ax.text(-7.5, 12. ,'Average Deviation=' + str(np.round(np.sqrt(sqD), 5)), fontsize='large')
        plt.show()
        muSlp = slp[i,j]
        print(p1[3])
        keep = input('do you wish to keep?y/n')
        if keep=='y':
            k1=0
            for key in keys:
                params[key].append(p1[k1])
                k1+=1
            params['Slp'].append(muSlp)
            k+=1
        if keep=='n':
            again = input('try again? y/n')
            if again=='y':
                continue
            elif again =='n':
                couplets['X'] = np.delete(couplets['X'],k)
                couplets['Y'] = np.delete(couplets['Y'],k)
                k+=1
    except:
        print('error somewhere')
        pdb.set_trace()
        continue

#np.save('rayParams.npy', params, allow_pickle=True)
#np.save('rayCouplets.npy', couplets, allow_pickle=True)
pdb.set_trace()
