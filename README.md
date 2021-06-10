# Topographic-Variance

This repository contains data and Python scripts for a manuscript titled: "Tree throw, topographic variance, and mapping process contributions from topography" which is currently under review for publication with Geophysical Research Letters. A preprint is available at:

This repository includes two folders that contain data and scripts for extracting topographic variance measurements from individual hillslopes and for parameterizing natural pit-mound couplets. 

Topographic Variance folder contains:
  1. a shapefile of hillslope locations, average slope, topographic variance, and standard deviation of slope
  2. hVa.py: a Python script that takes a DEM tile as an input and creates a high pass filter of topography, displays the topography and a satelite image so that the user may select hillslopes by clicking. Hillslopes measures are saved in a python dictionary. 
  3. CollectVarBrown2.npy: A python dictionary of hillslope characteristics

Parameterization folder contains:
  1. selectTile.py: A python script for selecting tiles that contain pit-mound couplets
  2. SelectTile2.py: A python script for visually inspecting the fit between idealized and natural couplets and allows the user to decide to keep the couplet or remove it from analysis
  3. rayParams.npy: a python dictionary with couplet parameters
  4. rayCouplets: a python dicionary with tiles containing couplets
  5. params.py: a python script for reading the dictionaries. 
