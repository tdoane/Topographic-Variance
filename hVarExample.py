import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from collections import Counter
from scipy.stats import mode
from scipy.signal import detrend
import os
import pdb
from collections import defaultdict
from osgeo import gdal
import ee
from skimage import io
cwd = os.getcwd()

from matplotlib import rc
#rc('text', usetex = True)
plt.rcParams.update({'font.size':14})
plt.rcParams.update({'font.family': 'serif'})

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

def contributingArea(z):
    #Define contributing area for extracting the network
    lY, lX = np.shape(z)
    grad =np.gradient(z)
    slp = np.sqrt(grad[0]**2 + grad[1]**2)
    fDir = np.arctan2(-grad[0],-grad[1])

    fDir /= np.pi/4
    fDir = np.round(fDir).astype('int')
    fDir[fDir==-1]=7
    fDir[fDir==-2]=6
    fDir[fDir==-3]=5
    fDir[fDir==-4]=4
    fDir[fDir==8] =0

    zVec = np.argsort(z.ravel('F'))
    k = len(zVec)-1
    iNeighbors = [0, 1, 1, 1, 0,-1,-1,-1]
    jNeighbors = [1, 1, 0,-1,-1,-1, 0, 1]

    area = np.ones_like(z)

    while k>-1:
        i= np.mod(zVec[k],lY)
        j = np.floor(zVec[k]/lY).astype('int')

        fDirInt = int(fDir[i,j])
        if (i>1 and i<(lY-2) and j>1 and j<(lX-2)):
            area[i+iNeighbors[fDirInt],j+jNeighbors[fDirInt]]+=area[i,j]
        k-=1
    return (area, fDir)

varRast={}
k=0

#Path to topographic data in this format
#Change this to your specific folder location
#pathname = '/home/tyler/Research/pitAndMound/Data/2019_03_26/025m_no0.tif'
pathname = 'griffyDEM.npy'

#This block of commands attempts to open an existing dataset if it exists. If it doesn't exists, it creates new data. 
try:
    temp=np.load('collectVarGriffy.npy', allow_pickle=True).item()
    CollectVar = defaultdict(list, temp)
    temp=np.load('CollectZGriffy.npy', allow_pickle=True).item()
    CollectZ=defaultdict(list, temp)
    #startTile = np.load('CurrentTile2.npy', allow_pickle=True)
    temp=np.load('polyPtsGriffy.npy', allow_pickle=True).item()
    polyPts = defaultdict(list, temp)
except:
    print('starting new dictionaries')
    CollectVar ={'fDir': [], 'var': [], 'domain': [], 'muSlp': [], 'stdSlp':[]}
    CollectZ={'hSlope':[]}
    polyPts=defaultdict(list)

hNum = len(CollectVar['var'])
dx = 0.25 #define the grid spacing (25 cm in this case)
#raster=gdal.Open(pathname) #open the DEM using gdal
#z = raster.ReadAsArray() #convert the raster in an array that numpy can handle (a Numpy Array)
z = np.load(pathname)
pdb.set_trace()

#ulX, xRes, xSkew, ulY, ySkew, yRes = raster.GetGeoTransform()#Get Extent of DEM in form of upper left, spacing, and offset

lY, lX = np.shape(z) #the dimensions of the domain (height and width)
x = np.arange(0, lX) #array of column numbers
y = np.arange(0, lY) #array of row numbers
XX, YY = np.meshgrid(x, y) #matrix of row and col positions
X, Y = XX.flatten(), YY.flatten() #turn those into 1-dimensional arrays
points = np.vstack((X, Y)).T #create an array of all  x,y pairs for the domain.

varMat = np.ones_like(z)*np.nan #initialize a matrix that records the variance of hillslopes

slpReal = np.sqrt((np.gradient(z)[0]/dx)**2+(np.gradient(z)[1]/dx)**2) #Calculate the magnitude ot the slope
zFilt = lowPass(z, 10) #perform a low pass filter
hPass = (z-zFilt) #subtract the low-pass filter from the topography leaving only roughness elements
grad = np.gradient(zFilt) #calculate the components of slope from the low pass filter
slp = np.sqrt(grad[0]**2+grad[1]**2) #calculate the magnitude of slope from the low pass filter

fDir = np.arctan2(-grad[0], -grad[1]) #Get flow directions using the low pass filter. Flow directions range from 0 (East) to 7 (NE), moving clock-wise
fDir[fDir<0.0]+=2*np.pi
FDIR=np.zeros_like(fDir)
FDIR[fDir>=np.pi/8]=1
FDIR[fDir>=3*np.pi/8]=2
FDIR[fDir>=5*np.pi/8]=3
FDIR[fDir>=7*np.pi/8]=4
FDIR[fDir>=9*np.pi/8]=5
FDIR[fDir>=11*np.pi/8]=6
FDIR[fDir>=13*np.pi/8]=7
FDIR[fDir>=15*np.pi/8]=0
FDIR[fDir<np.pi/8]=0

value = 'y'
while value=='y':
    fig= plt.figure(figsize=(20,20)) #create plot to visualize topography
    plt.imshow(hPass, vmin=-0.5, vmax=0.5) #show high pass topography 
    plt.contour(z, 20, colors = 'k', linewidths=0.5) #show contours 
    zoom_ok = False # this allows for us to zoom in and out before we pick a location
    while not zoom_ok: #once you press a button, then you will begin selecting points to define a polygon
        zoom_ok = plt.waitforbuttonpress()
    pts = plt.ginput(-1, timeout = -1) #collects the points
    p = Path(pts) #creates a collection of points that define a polygon
    grid = p.contains_points(points) #creates a mask of points within the polygon with '1' indicating it is within the polygon and '0' indicating it is outside of the polygon.
    plt.close() #close the figure
    mask = grid.reshape(lY, lX) #reshape the collectino
    inds = np.where(mask==1) #get locations inside polygon
    minI = np.min(inds[0]) #find minimum row
    maxI = np.max(inds[0]) #find maximum row
    minJ= np.min(inds[1])  #find minumum column
    maxJ=np.max(inds[1]) #find maximum column
    maskTemp = mask[minI:maxI,minJ:maxJ] #mask a rectangle of that size
    
    hSlope=hPass[minI:maxI, minJ:maxJ] #select the rectangle
    hSlope=hSlope*maskTemp #multiply the topography by the mask so that areas outside of the polygon are not included
    hSlope[hSlope==0]=np.nan #if things equal 0 assign them 'Not a Number' so that they are not included in the analysis

    hSlope[maskTemp=='False']==np.nan #label those points 'False'
    varTemp = np.nanvar(hSlope) #calculate the variance of the hillslope section
    fTemp = mode(FDIR[mask==1])[0][0] #select the mode of the aspect

    plt.imshow(hSlope, cmap = 'plasma') #show the selected hillslope
    plt.colorbar()
    plt.show()

    keep = input('keep? (y/n)') #if the hillslope looks good, print 'y' and it will record data. If not, type 'n' and it will not collect data.
    if keep=='y':
        CollectVar['fDir'] = np.append(CollectVar['fDir'], fTemp)
        CollectVar['var'] = np.append(CollectVar['var'], varTemp)
        CollectVar['domain']=np.append(CollectVar['domain'], np.sum(maskTemp.ravel()))
        CollectVar['muSlp']=np.append(CollectVar['muSlp'], np.average(slpReal[mask==1]))
        CollectVar['stdSlp']=np.append(CollectVar['stdSlp'], np.std(slpReal[mask==1]))
        CollectZ[str(hNum)]=hSlope.ravel()
        ptsTemp = np.array(pts)
        #polyPts[str(hNum)]=np.hstack((ulX+ ptsTemp[:,0]*dx, ulY - ptsTemp[:,1]*dx))
        polyPts[str(hNum)] = np.hstack(ptsTemp[:,0], ptsTemp[:,1])
        hNum+=1

        varMat[mask==1] = varTemp
        slpReal[mask==1] = varTemp
        hPass[mask==1]=varTemp

    print(hNum)
    value = input('continue (y/n)')

#varRast[os.path.splitext(filename)[0]]=varMat
np.save('polyPtsGriffy.npy', polyPts, allow_pickle=True)
np.save('CollectZGriffy.npy', CollectZ, allow_pickle=True)
np.save('CollectVarGriffy.npy', CollectVar, allow_pickle=True)
