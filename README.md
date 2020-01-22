# HOWTO

tps-georeference.py is a command line tool to transform images to georeferenced point-clouds through landmark definition based on thin-plate splines. works with qgis-.points-files

## usage:
    tps-georeference.py [input-path] [image-width] [image-height] [output-path]

## e. g.:
    tps-georeference.py test.points 1920 1080 tps.csv

# Before we start:
Create a landmark file (pixel coordinates to "world" coordinates). Therefore the QGIS Georeferencer data format *.points is adopted, so you may start with QGIS software and create your file there. For the case of Hamburg there was enough ground control point data available coming from geodata of the federel geo-portal. Then check what size is your source image (width and height) and save this as parameters in this script.

# Description: 
Thin-Plane spline warping of the input image (img) to output image (imgw). The warping is defined by a set of reference points (Zp=[Xp Yp]) on the [img] image and a set of corresponding points (Zs) on the (imgw) image. The warping will translate Zp in img to Zs in imgw.
 
## Input:
    img           - input image
	outDim        - Output canvas ([W H])
	Zp            - landmark in img
	Zs            - landmark in canvas
	interp.method - interpolation mode('nearest', 'invdist', 'none')
	interp.radius - Max radius for nearest neighbor interpolation or
					Radius of inverse weighted interpolation
	interp.power  - power for inverse weighted interpolation
 
## Output:
	imgw          - warped image with no holes
	imgwr         - warped image with holes
	map           - Map of the canvas with 0 indicating holes and 1 indicating pixel
 
## Reference:
F.L. Bookstein, "Principal Warps: Thin-Plate splines and the decomposition of deformations", IEEE Transaction on Pattern Analysis and Machine Intelligence, Vol. 11, No. 6, June 1989
 
## Author/ Acknowledgement: 
Fitzgerald J. Archibald
Date: 07-Apr-09

## Edit:
### From MATLAB to Python3:
CityScienceLab Hamburg
Date: 10/11-Apr-18
