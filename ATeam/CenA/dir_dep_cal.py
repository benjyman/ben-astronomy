#Direction-dpendent calibration with 'pre-peel'

#1. Take the best image of CenA
#2. Set the core to zero (-43 +- 7 arcmin)
#3. make a uvmodel with the core subtractedi (wsclean -predict)
#4  make a copy of the ms's
#5. subtract the core-less modell from the visibilities (use subtrmodel -usemodelcol)
#6. Self-cal on the core-only ms
#7. Reimage the 'just core' at higher angular res
#8. Do another self cal
#9. Apply these solutions to the 'unsubtracted ms'
#10 re-image at full image size


