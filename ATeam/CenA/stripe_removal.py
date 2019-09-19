#!/usr/bin/env python
#find and remove horizontal stripes in image
import matplotlib
matplotlib.use('Agg')

import cv2
import numpy as np
from astropy.io import fits

image_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0_polar.fits"
out_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0_polar_filtered.fits"

hdulist = fits.open(image_filename)
img = hdulist[0].data
#img = img[0,0,:,:]

#needs to be in Uint 8 for Canny
#Will need to do some thresholding/clipping and or select a 'quiet' part of the image 
#to be able to get sufficient dynamic range to detect the lines
#img_max = np.max(img)
#print img_max
#img_min = np.min(img)
#print img_min
#img_norm = img/img_max

img_clipped = np.clip(img,a_min=-0.5,a_max=1)

img_clipped_max = np.max(img_clipped)
img_clipped_norm = img_clipped/img_clipped_max

img_uint8 = (img_clipped_norm*255).astype(np.uint8)

#check conversion
test_filename = 'test.fits'
fits.writeto(test_filename,img_uint8,clobber=True)
print "saved %s" % test_filename


#img = img.astype(np.float32)

#img = cv2.imread(image_filename)
#gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
#edges = cv2.Canny(gray, 75, 150)
edges = cv2.Canny(img_uint8, 75, 150)


lines = cv2.HoughLinesP(edges, 1, np.pi/180, 30, maxLineGap=250)
for line in lines:
   x1, y1, x2, y2 = line[0]
   cv2.line(img, (x1, y1), (x2, y2), (0, 0, 128), 1)
cv2.imshow("linesEdges", edges)
cv2.imshow("linesDetected", img)
cv2.waitKey(10)
cv2.destroyAllWindows()

