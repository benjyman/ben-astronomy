#!/usr/bin/env python
#find and remove horizontal stripes in image
import matplotlib
matplotlib.use('Agg')

import cv2
import numpy as np

image_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0-MFS-image-pb.fits"
out_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0-MFS-image-pb_filtered.fits"


img = cv2.imread(image_filename)
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
edges = cv2.Canny(gray, 75, 150)
lines = cv2.HoughLinesP(edges, 1, np.pi/180, 30, maxLineGap=250)
for line in lines:
   x1, y1, x2, y2 = line[0]
   cv2.line(img, (x1, y1), (x2, y2), (0, 0, 128), 1)
cv2.imshow("linesEdges", edges)
cv2.imshow("linesDetected", img)
cv2.waitKey(0)
cv2.destroyAllWindows()

