#!/usr/bin/env python
#convert to polar coords
#see https://stackoverflow.com/questions/51675940/converting-an-image-from-cartesian-to-polar-limb-darkening
import matplotlib
matplotlib.use('Agg')

import cv2
import numpy as np

image_filename = "/md0/ATeam/CenA/stripes/polar/CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0-MFS-image-pb.fits"
out_filename = "CenA_2015_2018_joint_idg_12_obs_145_selfcal_03_robust0-MFS-image-pb_polar.fits"

source = cv2.imread(image_filename, 1)

#--- ensure image is of the type float ---
img = source.astype(np.float32)

#--- the following holds the square root of the sum of squares of the image dimensions ---
#--- this is done so that the entire width/height of the original image is used to express the complete circular range of the resulting polar image ---
value = np.sqrt(((img.shape[0]/2.0)**2.0)+((img.shape[1]/2.0)**2.0))

polar_image = cv2.linearPolar(img,(img.shape[0]/2, img.shape[1]/2), value, cv2.WARP_FILL_OUTLIERS)

polar_image = polar_image.astype(np.uint8)
cv2.imshow("Polar Image", polar_image)

cv2.waitKey(10)
cv2.destroyAllWindows()