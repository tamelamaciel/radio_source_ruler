#===============================================================================
# python script to automatically measure largest angular size from CoNFIG
# FITS files.
# Expects '/optical_positions.dat' as input file, and uses FIRST fits files
# stored in /data/

#-------------------------------------------------------------------------------
# Tamela Maciel Feb 2012
#===============================================================================
import os
import sys
import string
import math
import matplotlib as mpl
from shapely.geometry import Polygon
from shapely.geometry import Point
import fileinput
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import LogNorm
import numpy as np
import pyfits
from scipy import ndimage
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
import pywcsgrid2
import pywcs
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
#using locally downloaded gaussfit code from Adam Ginsberg
import mpfit
import gaussfitter as g
from matplotlib.patches import Ellipse

all=False
startwith=False
plot=False
if sys.argv[-1] == 'auto_size_measure.py':
  print('--------------------------------------------------------')
  print('Auto Size Measure')
  print('--------------------------------------------------------')
  print('Usage:')
  print('\'object name\' to process one galaxy,')
  print('\'object name:\' to start from given source,')
  print('\'all\' to process all sources in upload_file.')
  print('Add \'p\' to create and save plots.')
  print('ie: \'3C334 p\' or \'3C253:\' or \'all\'')
  print('--------------------------------------------------------')
  Input = string.split(raw_input('Input: '))
  print('--------------------------------------------------------')

  if len(Input)>1: 
    plot=True
  if Input[0]=='all':
    all=True
  elif Input[0].endswith(':'):
    startwith=True
    name=Input[0][:-1] #remove colon from name
  else: 
    name=Input[0]  
else:
  if sys.argv[1]=='all': 
    all=True
  elif sys.argv[1].endswith(':'):
    startwith=True
    name=sys.argv[1][:-1] #remove colon from name
  else: 
    name=sys.argv[1]
  if sys.argv[-1]=='p': 
    plot=True
    


######################################################################
########## FUNCTIONS #################################################
###################################################################### 

### 1. Read In Lists Of Data ###
def get_data(line):
  items=string.split(line)
  name=items[2]
  RA=items[4]+' '+items[5]+' '+items[6]
  DEC=items[7]+' '+items[8]+' '+items[9] #will have either minus sign ('-') or positive sign 
  #convert h:m:s to decimal degrees
  degRA=(float(items[4])+float(items[5])/60+float(items[6])/3600)*15.0 #in degrees. cos(dec) factor will be accounted for in sepn function
  if items[7].startswith('-'):
    radio_sign=-1. #plus or minus sign for declination
    #convert h:m:s to decimal degrees
    degDEC=radio_sign*(float(items[7].strip('-'))+float(items[8])/60+float(items[9])/3600) #in degrees
  else:
    radio_sign=1.
    #convert h:m:s to decimal degrees 
    degDEC=radio_sign*(float(items[7].strip('+'))+float(items[8])/60+float(items[9])/3600) #in degrees

  man_morph=items[15]
  redshift=items[10] #string value. Will have some '-' 
  RAO=items[35]+' '+items[36]+' '+items[37] 
  DEO=items[38]+' '+items[39]+' '+items[40] #will have some '-'
  
  optID=True
  try:
    degRAO=(float(items[35])+float(items[36])/60+float(items[37])/3600)*15.0
    if items[38].startswith('-'):
      opt_sign=-1.
      degDEO=opt_sign*(float(items[38].strip('-'))+float(items[39])/60+float(items[40])/3600)
    else:
      opt_sign=1.
      degDEO=opt_sign*(float(items[38].strip('+'))+float(items[39])/60+float(items[40])/3600)

  except ValueError: 
    degRAO=('-')
    degDEO=('-')
    optID=False

  manual_size=True
  try: 
    las=float(items[24]) 
    if man_morph in ['U','C','C*','S','S*']:
      print name+': check! las measure in catalog'
      las=str(las)+'*'
  except ValueError:
    try:
      las=float(items[29]) 
    except ValueError:
      las=items[24] #for sources that have no size information at all
      manual_size=False
      print name+': check! has no size info at all?'

  rms=float(items[91])/1000 #calculated rms from FITS file, using rms_fits.py script.

  return (name,RA,DEC,degRA,degDEC,man_morph,redshift,RAO,DEO,degRAO,degDEO,las,optID,rms)

### 2. FIND PEAKS IN IMAGE ###
def find_peaks(file):
  #read image data
  f=pyfits.open(file)
  img=f[0].data

  #set NaN pixels (empty pixels) to zero
  img[img != img]=0.0
  img[img<0.]=0.0
  if dim==4:
    img=img[0,0,:,:] #gets rid of 3 and 4 dimensions, since FIRST fits files have four axis but only first two have data
    
  T=ndimage.standard_deviation(img)
  sourcelabels,num_sources=ndimage.label(img>T)
  backgroundlabels,num_background=ndimage.label(img<T)
  # define an 8-connected neighbourhood
  neighborhood = generate_binary_structure(2,2)
  fimg=img*sourcelabels
  #apply the local maximum filter; all pixel of maximal value 
  #in their neighbourhood are set to 1
  local_max=maximum_filter(fimg,footprint=neighborhood)==fimg
  #In order to isolate the peaks we must remove the background from the mask.
  #we create the mask of the background
  background=img*backgroundlabels
  #we must erode the background in order to 
  #successfully subtract it form local_max, otherwise a line will 
  #appear along the background border (artifact of the local maximum filter)
  eroded_background=binary_erosion(background,structure=neighborhood,border_value=1)
  
  #we obtain the final mask, containing only peaks, 
  #by removing the background from the local_max mask
  detected_peaks=local_max-eroded_background
  #contains some peaks not in source (background bright features), but code can remove these
  #now need to find positions of these maximum
  #label peaks
  peaklabels,num_peaks=ndimage.measurements.label(detected_peaks)
  
  #get peak positions  
  slices = ndimage.find_objects(peaklabels)
  x, y = [], []
  for dy,dx in slices:
      x_center = (dx.start + dx.stop - 1)/2
      x.append(x_center)
      y_center = (dy.start + dy.stop - 1)/2    
      y.append(y_center)
      
  peak_positions=zip(x,y)
  
  #get peak values, in Jy/beam
  peak_fluxes=[]
  for coord in peak_positions:
      peak_fluxes.append(img[coord[1],coord[0]])
  peaks=zip(peak_positions,peak_fluxes)
  
  #sort by peak_fluxes. Two brightest peaks will be the first two in the list
  peaks=sorted(peaks,key=lambda l: l[1])  
  peaks.reverse()
  return peaks,img,f

### 3. GET FIRST LEVEL CONTOURS ###
def get_contours1(img,rms):
  maxpix = img.max()     # find maximum pixel value
  
  ### contours ###
  
  skip=False
  #if ratio of maxpix to noise is small (less than 40 ish)
  if maxpix/rms<300.:
    print name+' has small ratio of maxpix to noise'
    l=np.logspace(math.log(rms*3,2),math.log(maxpix,2),num=6,base=2.0) #log base 2 contour levels, starting at 5*rms value
  else:
    l=np.logspace(math.log(rms*5,2),math.log(maxpix,2),num=6,base=2.0) #log base 2 contour levels, starting at 5*rms value
  cl=plt.contour(img,levels=l,colors='w',origin='lower',interpolation='nearest')
  paths1=cl.collections[0].get_paths() #get outermost contour path
  contours1=[]
  for path in paths1: 
    contours1.append(Polygon(path.vertices)) #note that contours1 is a list of polygons which need to be converted to coords
  return contours1,skip,cl
    
    
### 4. TEST FOR GOOD CONTOUR ###
def isgood_contour(contour,goodcontour,rms):
  # for x in range(len(img)):
  #   for y in range(len(img[0])):  
  for y in range(len(img)):
    for x in range(len(img[0])):    
      #walk through fits file array until find point inside contour
      #if point in contour >30*rms, return goodcontour
      if contour.contains(Point(x,y)):
        #2*last_contour_level
        if img.max()/rms<40.:
          if img[y][x]>=10*rms:
            goodcontour.append(contour)
            return goodcontour
        else: 
          if img[y][x]>=30*rms:
            #print contour
            goodcontour.append(contour)
            return goodcontour
  #if contour doesn't contain bright point, return same goodcontour list as was fed to the function
  return goodcontour
        
### 5. COMPACT GAUSSIAN FIT FUNCTION ###
def compact_fit():
  params,gaussdata=g.gaussfit(img,returnfitimage=True) #gives height, amp, center_x, center_y, widthx, widthy, rotation
  height,amp,cx,cy,widthx,widthy,rot=params
  return height,amp,cx,cy,widthx,widthy,rot,gaussdata

### 6. SEPARATION FUNCTION ###
#takes ra and dec in degrees, without cosine(dec) factor on ra, but after multiplying by 15, and calculates distance in converted to degrees
def sepn(r1,d1,r2,d2):
  r1=r1*math.pi/180.
  r2=r2*math.pi/180.
  d1=d1*math.pi/180.
  d2=d2*math.pi/180.
  cos_sepn=math.sin(d1)*math.sin(d2)+math.cos(d1)*math.cos(d2)*math.cos(r1-r2)
  return math.acos(cos_sepn)*180./math.pi #in radians

### 7. FUNCTION TO FIND HOTSPOTS ###
def find_hotspots(R,D,optID,peaks,rms):
  if not isinstance(R,float):
    print 'RA DEC not defined'
    return

  wcs=pywcs.WCS(fitshead,naxis=2)
  #convert opt_degRA and opt_degDEC from WC to pixel
  #use optical ID coordinates for center, if available
  #might want to keep track of whether center coords are from optical or radio

  #calc center value in pixels
  center=wcs.wcs_sky2pix([[R,D]],1)[0]
  centerx=center[0]
  centery=center[1]

  #find hotspot a (brightest hotspot). peaks is a sorted tuple with the brightest two peaks first
  #take brightest peak, make sure it's away from the centre. Label hotspot a. It will be in good contour because it's the brightest
  #but check it's larger than 30*rms anyway:
  for peak in peaks:
    if peak[1]>30*rms:
      hs_x,hs_y=peak[0][0]-centerx,peak[0][1]-centery #hotspot coordinates, centered on the optical ID
      peak_dist=math.sqrt(hs_x**2+hs_y**2)
      if peak_dist > 3: # if distance between optical centre and hotspot is greater than 3px, then assume outside core
        HSa_flux=peak[1] #in Jy/beam
        HSa_x,HSa_y=hs_x,hs_y
        HSa_dist=peak_dist #picks brightest peak, >30*rms, and further than 3 px away from centre
        HSa=True
        break
    else:
      HSa_flux='-'
      HSa_dist='-'
      HSa=False
      
  #find hotspot b
  for peak in peaks:
    if peak[1]>30*rms:
      hs_x,hs_y=peak[0][0]-centerx,peak[0][1]-centery #hotspot coordinates, centered on the optical ID
      peak_dist=math.sqrt(hs_x**2+hs_y**2)
      if peak_dist > 2: #  if distance between optical centre and hotspot is greater than 3px, then assume outside core
        try: peak_peak_angle=math.acos((hs_x*HSa_x+hs_y*HSa_y)/(HSa_dist*peak_dist))
        except ValueError: peak_peak_angle=0 #from rounding errors, if vectors are the same, set angle between to 0
        if peak_peak_angle>math.pi/2:
          HSb_flux=peak[1] #in Jy/beam
          HSb_x,HSb_y=hs_x,hs_y
          HSb_dist=peak_dist #picks next brightest peak, that is >30*rms, further away than 3px, and is 90 degrees away from hotspot a
          HSb=True
          break
    else:
      HSb_flux='-'
      HSb_dist='-'
      HSb=False
      
  if HSa==True and HSb==True:
    #convert hotspot positions to world coordinates:
    HSa_px_x=HSa_x+centerx #coordinates of brightest hotspot in pixels
    HSa_px_y=HSa_y+centery
    HSa_deg=wcs.wcs_pix2sky([[HSa_px_x,HSa_px_y]],1)[0]
    HSa_deg_x=HSa_deg[0]  #coordinates of brightest hotspot in degrees
    HSa_deg_y=HSa_deg[1]
    HSa_length_deg=sepn(HSa_deg_x,HSa_deg_y,R,D) #HSa length in degrees (just multiply by 3600 for arcsec. sepn function should account for cos(dec)
      
    HSb_px_x=HSb_x+centerx #coordinates of second brightest hotspot in pixels
    HSb_px_y=HSb_y+centery
    HSb_deg=wcs.wcs_pix2sky([[HSb_px_x,HSb_px_y]],1)[0]
    HSb_deg_x=HSb_deg[0]  #coordinates of brightest hotspot in degrees
    HSb_deg_y=HSb_deg[1]
    HSb_length_deg=sepn(HSb_deg_x,HSb_deg_y,R,D) #HSb length in degrees (just multiply by 3600 for arcsec. sepn function should account for cos(dec)

  elif HSa==True and HSb==False:
    #convert hotspot positions to world coordinates:
    HSa_px_x=HSa_x+centerx #coordinates of brightest hotspot in pixels
    HSa_px_y=HSa_y+centery
    HSa_deg=wcs.wcs_pix2sky([[HSa_px_x,HSa_px_y]],1)[0]
    HSa_deg_x=HSa_deg[0]  #coordinates of brightest hotspot in degrees
    HSa_deg_y=HSa_deg[1]
    HSa_length_deg=sepn(HSa_deg_x,HSa_deg_y,R,D) #HSa length in degrees (just multiply by 3600 for arcsec. sepn function should account for cos(dec)
    
    #set HSb to '-'
    HSb_px_x,HSb_py_y='-','-'
    HSb_deg_x,HSb_deg_y='-','-'
    HSb_length_deg='-'
    
  else:
    HSa_px_x,HSa_px_y='-','-'
    HSa_deg_x,HSa_deg_y='-','-'
    HSa_length_deg='-'
    HSb_px_x,HSb_px_y='-','-'
    HSb_deg_x,HSb_deg_y='-','-'
    HSb_length_deg='-'
    
  return HSa_flux,HSa_px_x,HSa_px_y,HSa_deg_x,HSa_deg_y,HSa_length_deg,HSb_flux,HSb_px_x,HSb_px_y,HSb_deg_x,HSb_deg_y,HSb_length_deg,HSa_dist,HSb_dist,HSa,HSb

### FUNCTION TO CONVERT COMPACT FWHM to ANGULAR SIZE
def gauss_to_angsize(majorax,R,D,optID):
  if not isinstance(R,float):
    print 'RA DEC not defined'
    return
  wcs=pywcs.WCS(f[0].header,naxis=2)
  #convert opt_degRA and opt_degDEC from WC to pixel
  #use optical ID coordinates for center, if available
  center=wcs.wcs_sky2pix([[R,D]],1)[0]
  centerx=center[0] #id from catalog
  centery=center[1]
  
  halflength=majorax/2. #half length of major axis from gaussfit
  
  pix_x1=halflength*math.cos(rot*math.pi/180.)+cx #coord of fwhm of major axis fit, added to cx, the center from gaussfit
  pix_y1=halflength*math.sin(rot*math.pi/180.)+cy
  
  #use x1, y1 and cx, cy (the center found from gaussfit) to compute fwhm/2 in arcsec 
  arm_1=wcs.wcs_pix2sky([[pix_x1,pix_y1]],1)[0]    
  deg_x1=arm_1[0] 
  deg_y1=arm_1[1]  
  gauss_center=wcs.wcs_pix2sky([[cx,cy]],1)[0]
  deg_cx=gauss_center[0] #cx in degrees
  deg_cy=gauss_center[1]  
  ang_arm1=sepn(deg_x1,deg_y1,deg_cx,deg_cy)
  total_jet_length=(ang_arm1*2.)*3600. #total ellipse fwhm in arcsec  
   
  return centerx,centery,pix_x1,pix_y1,ang_arm1,total_jet_length

### 8. FUNCTION TO FIND LENGTH OF JETS ###
def find_dist(goodcontour,R,D,optID):
  if not isinstance(R,float):
    print 'RA DEC not defined'
    return
  wcs=pywcs.WCS(f[0].header,naxis=2)
  #convert opt_degRA and opt_degDEC from WC to pixel
  #use optical ID coordinates for center, if available
  #might want to keep track of whether center coords are from optical or radio
  center=wcs.wcs_sky2pix([[R,D]],1)[0]
  centerx=center[0]
  centery=center[1]

  max_coord=[] #in pixels in reference frame of ra_dec=(0,0)
  for contour in goodcontour:
    x=[coord[0] for coord in list(contour.exterior.coords)] #converts polygon object to coordinate list
    y=[coord[1] for coord in list(contour.exterior.coords)]
    #set center coordinates to (0,0)
    x=[element-centerx for element in x] 
    y=[element-centery for element in y] 
    max_dist1=0.
    for coord in zip(x,y):
      dist1=math.sqrt(coord[0]**2+coord[1]**2)
      if dist1>max_dist1:
        max_dist1=dist1 #largest arm length in pixels
        x1,y1=coord[0],coord[1]
        
        
    ## if only one continuous contour, find arm2 by looking >90 deg away:
    if len(goodcontour)==1:
      #convert pixels/measure distance for first arm:
      pix_x1=x1+centerx #coordinates of largest arm in pixels
      pix_y1=y1+centery
      arm_1=wcs.wcs_pix2sky([[pix_x1,pix_y1]],1)[0]
      deg_x1=arm_1[0] 
      deg_y1=arm_1[1]
      ang_arm1=sepn(deg_x1,deg_y1,R,D) #arm1 length in degrees (just multiply by 3600 for arcsec. sepn function should account for cos(dec)
      #measure second arm
      max_dist2=0.
      for coord in zip(x,y):
        dist2=math.sqrt(coord[0]**2+coord[1]**2)
        #angle between two arms using dot product
        try:
          angle=math.acos((coord[0]*x1+coord[1]*y1)/(max_dist1*dist2))
        except ValueError:
          angle=0 #from rounding errors, if vectors are the same, set angle between to 0
        if angle>math.pi/2:
          if dist2>max_dist2:
            max_dist2=dist2 #second arm length in pixels
            x2,y2=coord[0],coord[1]
            
      pix_x2=x2+centerx #coordinates of second arm in pixels
      pix_y2=y2+centery
      arm_2=wcs.wcs_pix2sky([[pix_x2,pix_y2]],1)[0]
      deg_x2=arm_2[0] 
      deg_y2=arm_2[1]
      ang_arm2=sepn(deg_x2,deg_y2,R,D)
      total_jet_length=(ang_arm1+ang_arm2)*3600 #total jet length in arcsec
      return centerx,centery,max_dist1,pix_x1,pix_y1,ang_arm1,max_dist2,pix_x2,pix_y2,ang_arm2,total_jet_length
          
    #otherwise, continue loop over contours in goodcontours to find max distances and add to max_coord list:      
    else:
      max_coord.append((x1,y1))
  ### calculate max distances from core for each coordinate in max_coord list ###
  max_dists=[math.sqrt(coord[0]**2+coord[1]**2) for coord in max_coord]
  coords_dists=zip(max_dists,max_coord)
  #sort by max distances value:
  coords_dists=sorted(coords_dists,key=lambda arm: arm[0])
  #pick last two (therefore largest two) dists and calc angle. if greater than 90 deg, return two arms
  try:
    angle=math.acos((coords_dists[-1][1][0]*coords_dists[-2][1][0]+coords_dists[-1][1][1]*coords_dists[-2][1][1])/(coords_dists[-1][0]*coords_dists[-2][0]))
  except ValueError:
    angle=0 #from rounding errors, if vectors are the same, set angle to 0.
  if angle>math.pi/2:
    x1,y1=coords_dists[-1][1][0],coords_dists[-1][1][1]
    max_dist1=coords_dists[-1][0]
    x2,y2=coords_dists[-2][1][0],coords_dists[-2][1][1]
    max_dist2=coords_dists[-2][0]
    
    #convert pixels/measure distance for first arm:
    pix_x1=x1+centerx #coordinates of largest arm in pixels
    pix_y1=y1+centery
    arm_1=wcs.wcs_pix2sky([[pix_x1,pix_y1]],1)[0] 
    deg_x1=arm_1[0]
    deg_y1=arm_1[1]
    ang_arm1=sepn(deg_x1,deg_y1,R,D) #arm1 length in degrees

    #compare_dist1=math.sqrt((deg_x1-degRA)**2+(deg_y1-degDEC)**2) 

    pix_x2=x2+centerx #coordinates of second arm in pixels
    pix_y2=y2+centery
    arm_2=wcs.wcs_pix2sky([[pix_x2,pix_y2]],1)[0]
    deg_x2=arm_2[0] 
    deg_y2=arm_2[1]
    ang_arm2=sepn(deg_x2,deg_y2,R,D) #arm2 length in degrees
    #compare_dist2=math.sqrt((deg_x2-degRA)**2+(deg_y2-degDEC)**2) 
    total_jet_length=(ang_arm1+ang_arm2)*3600 #total jet length in arcsec
    #degree_ratio=(total_jet_length/linear_size)*100/3600 #for 100 kpc, in degrees
    #size=wcs.wcs_sky2pix([[degRA+degree_ratio,degDEC]],1)[0]
    #radio_centerx=wcs.wcs_sky2pix([[degRA,degDEC]],1)[0]
    #radio_centerx=radio_centerx[0]
    #size=size[0]-radio_centerx #extract x length and subtract radio center(in pixels)
    return centerx,centery,max_dist1,pix_x1,pix_y1,ang_arm1,max_dist2,pix_x2,pix_y2,ang_arm2,total_jet_length

  #if angle is not greater than pi/2:
  else:
    #find next largest distance that is pi/2 away from max dist. if more than one good contour, will need to find max dist to appropriate contour
    try:
      angle=math.acos((coords_dists[-1][1][0]*coords_dists[-3][1][0]+coords_dists[-1][1][1]*coords_dists[-3][1][1])/(coords_dists[-1][0]*coords_dists[-3][0]))
    except ValueError:
      angle=0 #from rounding errors, if vectors are the same, set angle to 0.
    if angle>math.pi/2:
      x1,y1=coords_dists[-1][1][0],coords_dists[-1][1][1]
      max_dist1=coords_dists[-1][0]
      x2,y2=coords_dists[-3][1][0],coords_dists[-3][1][1]
      max_dist2=coords_dists[-3][0]
      
      #convert pixels/measure distance for first arm:
      pix_x1=x1+centerx #coordinates of largest arm in pixels
      pix_y1=y1+centery
      arm_1=wcs.wcs_pix2sky([[pix_x1,pix_y1]],1)[0]
      deg_x1=arm_1[0] 
      deg_y1=arm_1[1]
      ang_arm1=sepn(deg_x1,deg_y1,R,D) #arm1 length in degrees
      #compare_dist1=math.sqrt((deg_x1-degRA)**2+(deg_y1-degDEC)**2) 

      pix_x2=x2+centerx #coordinates of second arm in pixels
      pix_y2=y2+centery
      arm_2=wcs.wcs_pix2sky([[pix_x2,pix_y2]],1)[0]
      deg_x2=arm_2[0]
      deg_y2=arm_2[1]
      ang_arm2=sepn(deg_x2,deg_y2,R,D) #arm2 length in degrees
      #compare_dist2=math.sqrt((deg_x2-degRA)**2+(deg_y2-degDEC)**2) 
      total_jet_length=(ang_arm1+ang_arm2)*3600 #total jet length in arcsec
      #degree_ratio=(total_jet_length/linear_size)*100/3600 #for 100 kpc, in degrees
      #size=wcs.wcs_sky2pix([[degRA+degree_ratio,degDEC]],1)[0]
      #radio_centerx=wcs.wcs_sky2pix([[degRA,degDEC]],1)[0]
      #radio_centerx=radio_centerx[0]
      #size=size[0]-radio_centerx #extract x length and subtract radio center(in pixels)    
      return centerx,centery,max_dist1,pix_x1,pix_y1,ang_arm1,max_dist2,pix_x2,pix_y2,ang_arm2,total_jet_length
### end of dist function ###  
#######################################################################################################################################

### 9. FUNCTION TO ASSOCIATE HOTSPOTS WITH APPROPRIATE LOBE ###
def associate_hotspots_lobes(centerx,centery,pix_x1,pix_y1,pix_x2,pix_y2,HSa_px_x,HSa_px_y,HSb_px_x,HSb_px_y):
  #convert to coordinates centered on optical/radio ID
  x1=pix_x1-centerx #lobe1 coordinates
  y1=pix_y1-centery
  x2=pix_x2-centerx #lobe2 coordinates
  y2=pix_y2-centery
  if HSa==True:
    HSa_x=HSa_px_x-centerx #hotspota coordinates
    HSa_y=HSa_px_y-centery
    #check angle between arm1 and HSa
    try:
      angle=math.acos((HSa_x*x1+HSa_y*y1)/(max_dist1*HSa_dist))
    except ValueError:
      angle=0 #from rounding errors, if vectors are the same, set angle between to 0
    if angle<math.pi/2:
      #rename HSa info to HS1, corresponding to arm1 (longest arm)
      HS1_px_x=HSa_px_x
      HS1_px_y=HSa_px_y
      HS1_flux=HSa_flux
      HS1_deg_x=HSa_deg_x 
      HS1_deg_y=HSa_deg_y
      HS1_length_deg=HSa_length_deg
      #rename HSb info to HS2, corresponding to the opposite arm2 (shorter arm)
      HS2_px_x=HSb_px_x
      HS2_px_y=HSb_px_y
      HS2_flux=HSb_flux
      HS2_deg_x=HSb_deg_x 
      HS2_deg_y=HSb_deg_y
      HS2_length_deg=HSb_length_deg
    elif angle>math.pi/2:
      #rename HSa info to HS2, corresponding to arm2 (shorter arm)
      HS2_px_x=HSa_px_x
      HS2_px_y=HSa_px_y
      HS2_flux=HSa_flux
      HS2_deg_x=HSa_deg_x
      HS2_deg_y=HSa_deg_y
      HS2_length_deg=HSa_length_deg
      #rename HSb info to HS1, corresponding to arm1 (longest arm)
      HS1_px_x=HSb_px_x
      HS1_px_y=HSb_px_y
      HS1_flux=HSb_flux
      HS1_deg_x=HSb_deg_x
      HS1_deg_y=HSb_deg_y
      HS1_length_deg=HSb_length_deg
  else:
    HS1_px_x=HSa_px_x
    HS1_px_y=HSa_px_y
    HS1_flux=HSa_flux
    HS1_deg_x=HSa_deg_x
    HS1_deg_y=HSa_deg_y
    HS1_length_deg=HSa_length_deg

    HS2_px_x=HSb_px_x
    HS2_px_y=HSb_px_y
    HS2_flux=HSb_flux
    HS2_deg_x=HSb_deg_x
    HS2_deg_y=HSb_deg_y
    HS2_length_deg=HSb_length_deg
    
  return HS2_px_x, HS2_px_y,HS2_flux,HS2_deg_x, HS2_deg_y, HS2_length_deg, HS1_px_x, HS1_px_y, HS1_flux, HS1_deg_x, HS1_deg_y, HS1_length_deg

### 10. FUNCTION TO CALCULATE LUMINOSITY DISTANCE AND LINEAR SIZE ###
def calc_size(z,total_jet_length):
  try: 
    z=float(z)
    #initialize constants
    #assume flat Universe with Ho=71
    ################################################
    H0 = 71                       # Hubble constant
    WM = 0.27                     # Omega(matter)
    WV = 0.73                     # Omega(vacuum) or lambda [1.0 - WM - 0.4165/(H0*H0)]
    WK=0.                         # Omega curvaturve = 1-Omega(total)
    c = 299792.458                # velocity of light in km/sec
    ################################################
    WR = 0.        # Omega(radiation)
    Tyr = 977.8    # coefficent for converting 1/H into Gyr
    DTT = 0.5      # time from z to now in units of 1/H0
    DTT_Gyr = 0.0  # value of DTT in Gyr
    age = 0.5      # age of Universe in units of 1/H0
    age_Gyr = 0.0  # value of age in Gyr
    zage = 0.1     # age of Universe at redshift z in units of 1/H0
    zage_Gyr = 0.0 # value of zage in Gyr
    DCMR = 0.0     # comoving radial distance in units of c/H0
    DCMR_Mpc = 0.0 
    DCMR_Gyr = 0.0
    DA = 0.0       # angular size distance, radians 
    DA_Mpc = 0.0
    DA_Gyr = 0.0
    kpc_DA = 0.0
    DL = 0.0       # luminosity distance
    DL_Mpc = 0.0
    DL_Gyr = 0.0   # DL in units of billions of light years
    V_Gpc = 0.0
    a = 1.0        # 1/(1+z), the scale factor of the Universe
    az = 0.5       # 1/(1+z(object))
  
    h = H0/100.
    WR = 4.165E-5/(h*h)    # includes 3 massless neutrino species, T0 = 2.72528
    az = 1.0/(1+1.0*z)
    age = 0.
    n=1000                 # number of points in integrals
    for i in range(n):
      a = az*(i+0.5)/n
      adot = math.sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
      age = age + 1./adot
  
    zage = az*age/n
    zage_Gyr = (Tyr/H0)*zage #age at redshift z, in Gyr
    DTT = 0.0
    DCMR = 0.0
    ##############################################
    # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
    for i in range(n):
      a = az+(1-az)*(i+0.5)/n
      adot = math.sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
      DTT = DTT + 1./adot
      DCMR = DCMR + 1./(a*adot)
  
    DTT = (1.-az)*DTT/n
    DCMR = (1.-az)*DCMR/n
    age = DTT+zage
    age_Gyr = age*(Tyr/H0)
    DTT_Gyr = (Tyr/H0)*DTT
    DCMR_Gyr = (Tyr/H0)*DCMR
    DCMR_Mpc = (c/H0)*DCMR          #comoving radial distance in Mpc
  
    # tangential comoving distance
    ratio = 1.00
    x = math.sqrt(abs(WK))*DCMR
    y = x*x
    if WK < 0: y = -y
    ratio = 1. + y/6. + y*y/120.
    
    DCMT = ratio*DCMR
    DA = az*DCMT                #angular size distance (unitless)
    DA_Mpc = (c/H0)*DA            #angular size distance in Mpc
    kpc_DA = DA_Mpc/206.264806      #kpc/arcsec scale factor
    DL = DA/(az*az)
    DL_Mpc = (c/H0)*DL            #luminosity distance in Mpc
    DL_m = DL_Mpc*(10**6)*3.08*10**16 #luminosity distance in metres
    
    linear_size = total_jet_length*kpc_DA               #linear size conversion from arcsec to kpc
  except ValueError:
    #no redshift measured
    print 'No redshift available'
    linear_size='-'     
  return linear_size 
#################################################################    
    
### 11. PLOT FUNCTION ###
def plot_figure(f,skip,total_jet_length,linear_size,rms):

  from matplotlib import rcParams
  rcParams['axes.labelsize'] = 11
  rcParams['axes.unicode_minus'] = False
  rcParams['axes.formatter.limits']=-3,3
  rcParams['xtick.labelsize'] = 10
  rcParams['ytick.labelsize'] = 10
  rcParams['xtick.major.size']=3 
  rcParams['ytick.major.size']=3 
  rcParams['legend.fontsize'] = 10
  rcParams['font.family'] = 'serif'
  rcParams['font.serif'] = ['Palatino']
  rcParams['text.usetex'] = True
  rcParams['text.latex.preamble'] = \
      '\usepackage{amsfonts},' \
      '\usepackage{amsmath},' \
      '\usepackage[T1]{fontenc},' \
      '\usepackage{pxfonts},' \
      '\DeclareMathAlphabet{\mathbf}{OML}{pxmi}{bx}{it}'
  
  #clear axes
  plt.cla()
  plt.clf()

  def setup_axes(fig):
    ax=pywcsgrid2.subplot(111,header=f[0].header)
    divider=make_axes_locatable(ax)  
    # cax=divider.new_horizontal('5%',pad=0.05,axes_class=Axes)
    # fig.add_axes(cax)
    return ax #,cax

  if all==False and startwith==False:
    WIDTH = 412. # the width, in points, of a thesis page (textwidth) divided by three images across
    FACTOR = 0.65  # the fraction of the width you'd like the figure to occupy
    fig_width_pt  = WIDTH * FACTOR
    inches_per_pt = 1.0 / 72.27
    fig_width_in  = fig_width_pt * inches_per_pt  # figure width in inches
    fig_height_in = fig_width_in   # figure height in inches. want square image
    fig_dims      = [fig_width_in,fig_width_in]#[5.7,5.7] #dim in inches
  else:
    # images scaled for a directory
    fig_dims      = [5.7,5.7]

  #prepare figure and axes
  fig = plt.figure(figsize=fig_dims)
  ax=setup_axes(fig)
  #ax.set_display_coord_system('fk5')
  
  img=f[0].data
  #set NaN pixels (empty pixels) to zero
  img[img != img]=0.0
  if dim==4:
    img=img[0,0,:,:] 
  
  #Jy-> mJy
  data=img*1000
  rms=rms*1000
  maxpix=data.max()
  minpix=data.min()

  # ### Uncomment to show color image as well
  # #draw image
  # im=ax.imshow(data,cmap=plt.cm.cubehelix_r,vmin=0.01,vmax=maxpix,origin='lower',interpolation='nearest') #could use cmap=plt.cm.cubehelix, but would want to start in the middle (ie not at black or white). original RdBu_r. another good option YlGnBu_r
  # #set color limits
  # im.set_clim(0.01,maxpix)

  #draw contour
  #if ratio of maxpix to noise is small (less than 40 ish):
  if maxpix/rms<300.:
    l=np.logspace(math.log(rms*3,2),math.log(maxpix,2),num=6,base=2.0) #log base 2 contour levels
  else:
    l=np.logspace(math.log(rms*5,2),math.log(maxpix,2),num=6,base=2.0) #log base 2 contour levels
  cl=ax.contour(data,l,colors='#000066',alpha=1,origin='lower',interpolation='nearest')
  for col in cl.collections: 
    col.set_linewidth(0.75)
  # cbar=plt.colorbar(im,cax=cax)
  # #adjust cbar ticks and add levels for contour lines
  # #cbar.set_ticks(l)
  # cbar.add_lines(cl)
  
  #mark centre (use optical ID otherwise radio)
  if optID==True:
    wcs=pywcs.WCS(f[0].header,naxis=2)
    opt_center=wcs.wcs_sky2pix([[R,D]],1)[0]
    opt_centerx=opt_center[0]
    opt_centery=opt_center[1]
    ax.plot(opt_centerx,opt_centery,color='#00FF00',marker='+',markersize=6,markeredgecolor='#00FF00',markeredgewidth=1.2,zorder=6)
  else:
    wcs=pywcs.WCS(f[0].header,naxis=2)
    radio_center=wcs.wcs_sky2pix([[R,D]],1)[0]
    radio_centerx=radio_center[0]
    radio_centery=radio_center[1]
    ax.plot(radio_centerx,radio_centery,color='#DC143C',marker='+',markersize=6,markeredgecolor='#DC143C',markeredgewidth=1.2,zorder=6) #if no optical ID, use centerx (radio assigment)
  f.close()
  #mark arms  
  if skip==False and compact==False:
    ax.plot(pix_x1,pix_y1,color='#FF0000',marker='x',markersize=6,markeredgecolor='#FF0000',markeredgewidth=1.2)
    ax.plot(pix_x2,pix_y2,color='#FFD700',marker='x',markersize=6,markeredgecolor='#FFD700',markeredgewidth=1.2)
  
  #mark hotspots
  if skip==False and compact==False:
    try:
      ax.plot(HS1_px_x,HS1_px_y,color='#FF0000',marker='x',markersize=6,markeredgecolor='#FF0000',markeredgewidth=1.2)  #red
    except ValueError:
      pass
    try:
      ax.plot(HS2_px_x,HS2_px_y,color='#FFD700',marker='x',markersize=6,markeredgecolor='#FFD700',markeredgewidth=1.2)  #yellow
    except ValueError:
      pass
  ##add scale bar for 100 kpc
  try:
    deg_100kpc_ratio=(total_jet_length/linear_size)*100/3600. #100 kpc to degrees
    ax.add_size_bar(deg_100kpc_ratio/0.0005,label='100 kpc',loc=4) #specifies length in pixels using axis coord increment from fits header. 
  except TypeError: pass #no scalebar since redshift='-' and thus linear_size='-'
  
  if compact==True:
    ##add ellipse for gaussian fit
    hm=[gaussdata.max()/2.]
    ax.contour(gaussdata,hm,colors='g',lw=4.,alpha=1,origin='lower',interpolation='cubic',zorder=5)
    #ellipse=Ellipse(xy=(cx,cy), width=fwhmx, height=fwhmy, angle=rot,edgecolor='r', fc='None', lw=1)
    #ax.add_patch(ellipse)
    
    
  ##add beam size
  #(major,minor)=5.4",5.4", angle=0.0
  ax.add_beam_size(0.0015/0.0005,0.0015/0.0005,0.0,loc=3,patch_props={'facecolor':'None','edgecolor':'k'}) #beam size in deg / axis coord increment from fits header. 'hatch':'////'
  
  # #Figure title
  # ax.add_inner_title(name,loc=2)
  
  ## set axes labels
  ## cax.set_ylabel('mJy/Beam',fontdict={'fontname':'Bitstream Vera Sans','fontsize':11},labelpad=0)
  # cax.set_ylabel('mJy/Beam')

  ax.axis['bottom'].label.set_text("Right Ascension (J2000)")
  ax.axis['left'].label.set_text("Declination (J2000)")
  # ax.axis["bottom"].label.set_visible(False)
  # ax.axis["left"].label.set_visible(False)
  #save image
  filename=name+'_size.pdf'
  #ax.set_rasterization_zorder(2.1) #not sure what this does. From IC443 radio continuum code using pywcsgrid2
  #cax.set_rasterization_zorder(2.1)
  plt.tight_layout(pad=0.1)
  plt.savefig('autosize_images/'+filename)
  plt.close('all')
  #plt.show()
  return 


#######################################################################################################################################
################## BEGIN PROGRAM ######################################################################################################
#######################################################################################################################################

### Input file ###
upload_file=open('optical_positions.dat', 'r') 
##################
outfile='auto_sizes.dat'
################################################################################


### Exceptions ###
### List of exceptions, where morphology is determined by NVSS, and resolved out by FIRST, or need bigger image cutout of first ### 
exceptions=['3C230']


lines=upload_file.readlines()
upload_file.close()
print '----------------------------------------------------'
print '----------------------------------------------------'
### need to read in data from definitive_CONFIG for just one source if all=False
if all==False:
  for i, s in enumerate(lines):
    if name in s:
      index=i
      #if startwith=False,get only specific line of data for source
      if startwith==False: 
        line=lines[i]
        lines=[line]
      if startwith==True:
        lines=lines[i:]
      break  
if all==True:
  #write create new auto_sizes.dat file and write header
  ######### Write header #########################
  fileout = open(outfile, 'w')
  fileout.write('#Name RAJ2000 - - DEJ2000 - - Redshift Morph RAO - - DEO - - Man_ang_size/arcsec Auto_ang_size/arcsec Auto_lin_size/kpc HS1/Lobe1_ratio HS2/Lobe2_ratio Man_morph Auto_Morph Check_contours Note_(r=radioID)\n')
  fileout.close()
  ################################################
  #write header for fail list
  failfile=open('fail_list_auto_sizes.dat','w')
  failfile.write('#Fail list from auto_size_measure.py')
  failfile.close()
  ################################################


### iterate for element in lines (which may have len 1, part of total, or total lines of upload_file)      
for row in lines:
  #skip first header line with '#':
  if '#' in row:
    continue 
  name,RA,DEC,degRA,degDEC,man_morph,redshift,RAO,DEO,degRAO,degDEO,las,optID,rms = get_data(row)

  file='FITS_files/'+name+'.fits'
  dim=4
  
  
  #get peak positions and values in image (will be narrowed down to hotspots later on)
  peaks,img,f=find_peaks(file)

  fitshead=f[0].header

  #center coords (either optical or radio)
  if optID==True:
    R=degRAO
    D=degDEO
  if optID==False:
    R=degRA
    D=degDEC

  #find first level contours
  contours1,skip,cl=get_contours1(img,rms)
  check_contours='-'
  if skip==False: 
    goodcontour=[]
    ### find goodcontours list from first level contours ###
    for contour in contours1:
      goodcontour=isgood_contour(contour,goodcontour,rms)
    #check for lots of spurious flux that may need cleaning. rerun contours
    if len(goodcontour)>6:
      print 'Might have lots of spurious flux'
      check_contours='c'
      new_rms=rms*2
      contours1,skip,cl=get_contours1(img,new_rms)
      goodcontour=[]
      for contour in contours1:
        goodcontour=isgood_contour(contour,goodcontour,rms)   
           
    if man_morph in ['C','C*','S','S*']:
      compact=True
      print name+': compact. fitting gaussian' 
      height,amp,cx,cy,widthx,widthy,rot,gaussdata=compact_fit()

      fwhmx=2.*math.sqrt(2.*math.log(2.))*widthx
      fwhmy=2.*math.sqrt(2.*math.log(2.))*widthy
      if fwhmx>fwhmy:
        majorax=fwhmx
      else:
        majorax=fwhmy
      HS1_arm1_ratio='-'
      HS2_arm2_ratio='-' 
      auto_morph='-'
      
      #convert fwhm to angular size (total_jet_length)
      centerx,centery,pix_x1,pix_y1,ang_arm1,total_jet_length=gauss_to_angsize(majorax,R,D,optID)
      ### convert angular size to linear size ###
      linear_size=calc_size(redshift,total_jet_length) #linear size of jet in kpc, using script from Schombert
      
    else:  
      ### find hotspots a and b for all non-compact sources
      compact=False
      try:
        HSa_flux,HSa_px_x,HSa_px_y,HSa_deg_x,HSa_deg_y,HSa_length_deg,HSb_flux,HSb_px_x,HSb_px_y,HSb_deg_x,HSb_deg_y,HSb_length_deg,HSa_dist,HSb_dist,HSa,HSb=find_hotspots(R,D,optID,peaks,rms)  
      except UnboundLocalError:
        failline=' '.join([name,'could not find hotspots','\n'])
        #write to fail file
        failfile = open('fail_list_auto_sizes.dat', 'a')
        failfile.write(failline)
        failfile.close()
        print '----------------------------------------------------'
        print '----------------------------------------------------'
        print name+': could not find hotspots. adding to fail list'
        print '----------------------------------------------------' 
        total_jet_length='-'
        linear_size='-'
        HS1_arm1_ratio='-'
        HS2_arm2_ratio='-' 
        auto_morph='-'
        skip=True
      
    if skip==False and compact==False:
      ### find lengths of two max distances of goodcontour list ###
      try:
        centerx,centery,max_dist1,pix_x1,pix_y1,ang_arm1,max_dist2,pix_x2,pix_y2,ang_arm2,total_jet_length=find_dist(goodcontour,R,D,optID)
        ### associate hotspots a and b with arms 1 and 2: ###
        HS2_px_x,HS2_px_y,HS2_flux,HS2_deg_x,HS2_deg_y,HS2_length_deg,HS1_px_x,HS1_px_y,HS1_flux,HS1_deg_x,HS1_deg_y,HS1_length_deg=associate_hotspots_lobes(centerx,centery,pix_x1,pix_y1,pix_x2,pix_y2,HSa_px_x,HSa_px_y,HSb_px_x,HSb_px_y)
        
        ### calculate hotspot/lobe length ratio ###
        try:
          HS1_arm1_ratio=HS1_length_deg/ang_arm1
          HS2_arm2_ratio=HS2_length_deg/ang_arm2
          ### from here, could estimate FRI or FRII based on this ratio
          if HS1_arm1_ratio>0.5 and HS2_arm2_ratio>0.5:
            auto_morph='II'
          elif HS1_arm1_ratio<0.5 and HS2_arm2_ratio<0.5:
            auto_morph='I'
          else:
            auto_morph='U/Hybrid'
        except TypeError:
          HS1_arm1_ratio,HS2_arm2_ratio='-','-'
          auto_morph='-'
        
        ### convert angular size to linear size ###
        linear_size=calc_size(redshift,total_jet_length) #linear size of jet in kpc, using script from Schombert
        angle='good'
        
      except TypeError:
        #angle between two max arms is less than 90
        angle='toosmall'
        total_jet_length='-'
        linear_size='-'
        HS1_arm1_ratio='-'
        HS2_arm2_ratio='-'  
        
        failline=' '.join([name,'angle < 90 degrees','\n'])
        #write to fail file
        failfile = open('fail_list_auto_sizes.dat', 'a')
        failfile.write(failline)
        failfile.close()
        print '----------------------------------------------------'
        print '----------------------------------------------------'
        print name+': angle between two longest arms is < 90 degrees'
        print '----------------------------------------------------' 
      

  ### otherwise don't calc total jet length ###
  else:
    print name+': not calculating jet length. source is an exception'   
    total_jet_length='-'
    linear_size='-'
    HS1_arm1_ratio='-'
    HS2_arm2_ratio='-'
    auto_morph='-'
    
  #close fits file
  f.close()
  
  ### IF MULTIPLE SOURCES, WRITE TO FILE ###
  #note if startwith=True, existing outfile will be appended
  if all==True or startwith==True:
    #create joined string separated by single space to print to outfile
    if optID==False:
      #append 'r' note
      line=' '.join([name,RA,DEC,redshift,man_morph,RAO,DEO,str(las),str(total_jet_length),str(linear_size),str(HS1_arm1_ratio),str(HS2_arm2_ratio),man_morph,auto_morph,check_contours,'r','\n'])
    elif optID==True:
      line=' '.join([name,RA,DEC,redshift,man_morph,RAO,DEO,str(las),str(total_jet_length),str(linear_size),str(HS1_arm1_ratio),str(HS2_arm2_ratio),man_morph,auto_morph,check_contours,'r','\n'])
    else:
      line=' '.join([name,RA,DEC,redshift,man_morph,RAO,DEO,str(las),str(total_jet_length),str(linear_size),str(HS1_arm1_ratio),str(HS2_arm2_ratio),man_morph,auto_morph,check_contours,'-','\n'])
    #write to output file
    fileout = open(outfile, 'a')
    fileout.write(line)
    fileout.close()
  
  ### PLOT FILE IF PLOT=TRUE
  if plot==True:
    f=pyfits.open(file)
    plot_figure(f,skip,total_jet_length,linear_size,rms)  
    f.close()
    
  #print progress in terminal
  print 'Name: '+name
  print 'Total jet length: '+str(total_jet_length)+' arcsec, '+str(linear_size)+' kpc'
  print 'Hotspot lobe ratio: Lobe 1: '+str(HS1_arm1_ratio)+', Lobe 2: '+str(HS2_arm2_ratio)
  print 'Morphology: Manual: '+man_morph+', Auto: '+auto_morph
  if optID==False: print 'No available optical ID, using radio ID'
  if plot==True: print 'Plotted and saved in /autosize_images/'
  print '----------------------------------------------------'
  print '----------------------------------------------------'
    
### END OF PROGRAM ###
