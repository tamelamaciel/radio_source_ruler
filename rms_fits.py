#===============================================================================
# python script to calculate the standard deviation of background pixels
# Masks source pixels to do this.
# Expects 'optical_positions.dat' as input file, and uses FIRST fits files
# stored in FITS_files/
#-------------------------------------------------------------------------------
# Tamela Maciel Jun 2012
#===============================================================================
import sys
import string
import math
import numpy as np
import pyfits
from scipy import ndimage

#################### INPUT and OUTPUT FILES ####################################
sourcefile='optical_positions.dat'
outfile='optical_positions_w_rms.dat'
################################################################################

######### Write header #########################
fileout = open(outfile, 'w')
fileout.write('#Sample CoNFIG Name RAJ2000 - - DEJ2000 - - S_1.4/mJy Alpha Type note_type uncertainty_type Redshift Error Z source_z RAO - - DEO - - Length/cm Scale (arcmin/cm) Length/arcsec S(178MHz)/Jy (for 3CRR only) Deconvolved_Maj_axis/arcsec (for C or U) LAS/arcsec(manual measure for C or U) rmag imag zmag Ksmag first_rms/mJy/bm calc_rms/Jy/beam Array\n')
fileout.close()
################################################
  
    
      
########## VLA exceptions ######################
VLAA_files=['3C230']
     
       

################################################
  
file  = open(sourcefile, 'r')

while 1:
  line=file.readline()
  #skip first header line with 'Name'
  if 'Name' in line:
    continue
  if not line: break
  items=string.split(line)
  name=items[2]
  fitsfile='FITS_files/'+name+'.fits'
  dim=4
  #read image data
  f=pyfits.open(fitsfile)
  img=f[0].data
  f.close()
  #set NaN pixels (empty pixels) to zero
  img[img != img]=0.0
  img[img<0.]=0.0
  if dim==4:
    img=img[0,0,:,:] #gets rid of 3 and 4 dimensions, since FIRST fits files have four axis but only first two have data
  
  ### calculate standard deviation from fits file ###
  T=ndimage.standard_deviation(img)
  sourcelabels,num_sources=ndimage.label(img>T)
  backgroundlabels,num_background=ndimage.label(img<T)
  #we create the mask of the background
  background=img*backgroundlabels
  #mean of background
  meanback=background.mean()
  #standard deviation of background
  stdev=ndimage.standard_deviation(background) #in Jy/beam
  #calc rms (= sqrt(mean^2+stdev^2)
  rms=math.sqrt(meanback**2+stdev**2)

  items.append(str(rms))
  items.append('FIRST')
  items.append('\n')
  #create joined string separated by single space to print to outfile
  line=' '.join(items)
  #write to output file
  fileout = open(outfile, 'a')
  fileout.write(line)
  fileout.close()
  print '----------------------------------------------------'
  print name+' rms calculated'
  print 'rms: '+str(rms)
  print '----------------------------------------------------'
file.close()