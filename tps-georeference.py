# --- HOWTO

# tps-georeference.py is a command line tool to transform images to georeferenced point-clouds through landmark definition based on thin-plate splines. works with qgis-.points-files

# usage:
# tps-georeference.py [input-path] [image-width] [image-height] [output-path]

# e. g.:
# tps-georeference.py test.points 1920 1080 tps.csv

'''
# Description: Thin-Plane spline warping of the input image (img) to
# output image (imgw). The warping is defined by a set of reference points
# (Zp=[Xp Yp]) on the [img] image and a set of corresponding points (Zs)
# on the (imgw) image. The warping will translate Zp in img to Zs in imgw.
# 
# Input:
#	img - input image
#	outDim - Output canvas ([W H])
#	Zp - landmark in img
#	Zs - landmark in canvas
#	interp.method - interpolation mode('nearest', 'invdist', 'none')
#	interp.radius - Max radius for nearest neighbor interpolation or
#					Radius of inverse weighted interpolation
#	interp.power - power for inverse weighted interpolation
# 
# Output:
#	imgw - warped image with no holes
#	imgwr - warped image with holes
#	map - Map of the canvas with 0 indicating holes and 1 indicating pixel
# 
# Reference:
# F.L. Bookstein, "Principal Warps: Thin-Plate splines and the
# decomposition of deformations", IEEE Transaction on Pattern Analysis and
# Machine Intelligence, Vol. 11, No. 6, June 1989
# 
# Author: Fitzgerald J Archibald
# Date: 07-Apr-09
'''

# --- NEW

# translated to python-3: TM@CSL
# on 10/11-Apr-18

# --- used modules

import sys # for command line args handling

from datetime import datetime

import numpy as np
from numpy import matlib, linalg

import csv

# --- initialize!

 # sys.argv[0] would be name of py-file
input = sys.argv[1] # path
imgW = int(sys.argv[2]) # img dim (width)
imgH = int(sys.argv[3]) # img dim (height)
output = sys.argv[4]

# print(input)
# print(imgW)
# print(imgH)

xraw_, yraw_, uraw_, vraw_, enable = np.loadtxt(input, delimiter=',', skiprows=1, unpack=True)

Xs = np.asmatrix(yraw_).getH() # flip x and y -> geom 2 math
Ys = np.asmatrix(xraw_).getH()

Xp = np.asmatrix(1080+vraw_).getH() # flip u and v, invert v so that [0,0] is bottom left
Yp = np.asmatrix(uraw_).getH()

NPs = len(Xp)

# --- functions

def radialBasis(ri):
	ri_ = ri
	
	for i in range(0, ri.size):
		if ri.item(i) == 0:
			j = i # +1
			u = j%ri.shape[1]
			v = int(np.floor(j/ri.shape[1]))
			
			ri_[v,u] = 1. # np.finfo(ri_.dtype).min # avoid log(0)=inf
	
	ko = np.multiply(
			2* np.power(ri,2),
			np.log(ri_)
	)
	
	return ko

def computeWl(xp, yp, n):
	wP_ = np.matrix(np.zeros([len(xp), 3]))

	rXp = matlib.repmat(xp, 1, n)
	rYp = matlib.repmat(yp, 1, n)
	
	wR = np.sqrt(np.power(rXp-rXp.getH(),2) + np.power(rYp-rYp.getH(),2)) # .getH() is the "complex conjugate transpose"
	
	wK = radialBasis(wR)
	
	for i in range(0,len(xp)):
		wP_[i] = [1, xp[i], yp[i]]
	
	wP = wP_.getH() # .getH() is the "complex conjugate transpose"
	
	wL_01 = np.vstack((wK, wP))
	
	wL_02 = np.vstack((wP_, np.zeros((3,3))))
	
	wL = np.hstack((wL_01, wL_02))
	
	return wL

def tpsMap(wW, imgH, imgW, xp, yp, n):
	X_, Y_ = np.meshgrid(np.linspace(1,imgH,imgH), np.linspace(1,imgW,imgW))
	
	X = X_.flatten()
	Y = Y_.flatten()
	
	NWs = len(X)
	
	rX = matlib.repmat(X,n,1)
	rY = matlib.repmat(Y,n,1)
	
	rxp = matlib.repmat(xp,1,NWs)
	ryp = matlib.repmat(yp,1,NWs)
	
	wR = np.sqrt(np.power(rxp-rX,2) + np.power(ryp-rY,2))
	
	wK = radialBasis(wR)
	
	wP_ = np.matrix(np.zeros([len(X), 3]))
	
	for i in range(0,len(X)):
		wP_[i] = [1,X[i],Y[i]]
	
	wP = wP_.getH()
	
	wL = np.hstack((wK.getH(), wP.getH()))
	
	Xw = wL * wW[:,0] # np.multiply(wL, wW[:,0])
	Yw = wL * wW[:,1] # np.multiply(wL, wW[:,1])
	
	return [Xw, Yw]

# --- main

print(str(datetime.now()) + ': start crunching data...')

# --- 1/3
# --- algebra of thin-plate splines: *

wL = computeWl(Xp, Yp, NPs) # *

# print(wL) # debug: compare w/ MATLAB

wY_ = np.hstack((Xs, Ys))

wY = np.vstack((wY_, np.zeros((3, 2)))) # *

wW = linalg.inv(wL)*wY # *

# print(wW) # debug: compare w/ MATLAB

# --- thin-plate spline mapping (map all points in the plane)

print(str(datetime.now()) + ': map all points in the plane...')

[Xw, Yw] = tpsMap(wW, imgH, imgW, Xp, Yp, NPs)

# --- 2/3
# --- warping

# --- input grid for warping

print(str(datetime.now()) + ': warp...')

[X_01, Y_01] = np.meshgrid(np.linspace(imgH,1,imgH), np.linspace(1,imgW,imgW)) # HxW

X = np.matrix(np.zeros([imgH*imgW, 1]))

X_02 = X_01.flatten()

for i in range(0, imgH*imgW):
	X[i] = (imgH+1) - X_02[i] # set pixel-(0,0) to bottomleft

Y = np.matrix(np.zeros([imgH*imgW, 1]))

Y_02 = Y_01.flatten()

for i in range(0, imgH*imgW):
	Y[i] = Y_02[i]

# --- needed for OLD (lines)
XwYwRC = np.hstack((Yw, Xw, Y, X)) # u = Y, v = X # X, Y [math.] will be X, Y [geom.]

print(str(datetime.now()) + ': calculation succeded.')

# --- 3/3
# --- OLD (lines)

print(str(datetime.now()) + ': write to file...')

with open(output, 'w') as out:
	np.savetxt(out, XwYwRC, comments='', delimiter=',', fmt=('%.9f','%.9f','%d','%d'), header='x,y,u,v')

print(str(datetime.now()) + ': values written to file.')

'''
# --- NEW (matrix)
table = np.zeros((imgH+1,imgW+1), dtype=object) # 1081x1921, so [0]-Indices remain empty

# fill table's unused indices
for i in range(0,imgW+1):
	table[0][i] = np.array([0.0, 0.0])

for i in range(0,imgH+1):
	table[i][0] = np.array([0.0, 0.0])

# fill table with data
for i in range(0, imgH*imgW):
	tempX = Yw.item(i) # X = Y
	tempY = Xw.item(i) # Y = X
	tempU = Y.item(i)
	tempV = X.item(i)
	
	table[int(tempV)][int(tempU)] = np.array([tempX, tempY]) # np.array([tempX, tempY])

print(str(datetime.now()) + ': write to file...')

with open(output, 'w') as out:
	np.savetxt(out, table, comments='', delimiter=',', fmt=('%s'))

print(str(datetime.now()) + ': values written to file.')
'''
