# -*- coding: utf-8 -*-
"""
Created on Thu Nov 07 16:00:58 2013

@author: Ben Hudson
"""
import os
import arcpy
import matplotlib.pyplot as plt
import scipy.ndimage
import numpy as np

from datetime import date
import scipy
import scipy.stats
from scipy.interpolate import griddata
import mpl_toolkits.basemap.pyproj as pyproj
import matplotlib.cm as cm

from astropy.convolution import convolve, Tophat2DKernel

from osgeo import gdal
drv = gdal.GetDriverByName('GTiff')

arcpy.env.overwriteOutput = True

sr = arcpy.SpatialReference("WGS 1984 UTM Zone 22N")

# ----------------------------------------------------------
# ----------------------------------------------------------
# FUNCTIONS 
# ----------------------------------------------------------
# ----------------------------------------------------------

def imshowwithZ(X):

    import matplotlib.cm as cm
    
    fig, ax = plt.subplots()
    
    
    ax.imshow(X, cmap=cm.jet, interpolation='nearest',vmin=400,vmax=4000)
    
   
    
    numrows, numcols = X.shape
    def format_coord(x, y):
        col = int(x+0.5)
        row = int(y+0.5)
        if col>=0 and col<numcols and row>=0 and row<numrows:
            z = X[row,col]
            return 'x=%1.4f, y=%1.4f, z=%1.4f'%(x, y, z)
        else:
            return 'x=%1.4f, y=%1.4f'%(x, y)
    
    ax.format_coord = format_coord
    #ax.set_cmap('spectral')
    plt.show()

def myround(x,base=10):
    return np.int(base * np.round(np.float(x)/base))
    
def smoothingFunction(x):
    
    #print x
    #aa = scipy.stats.nanmean(scipy.stats.nanmean(x))
    # if doesnt work try this 
    
    # trim the highest x value
    x = x[x < np.max(x)]

    # interquartile mean
    #x = x[np.logical_and(x < np.percentile(x,75),x > np.percentile(x,25))]
    
    bottom = np.size(x) - np.sum(np.isnan(x))
    
    
    if bottom == 0:
        aa = np.NaN
        
    
    else:
        
        aa = np.nansum(x) / bottom
    
    return aa
    #return (x*0.5).sum()

    
def L7toReflectance(landsatImageSubsetDN,band,L_MAX,L_MIN,QCAL_MAX,QCAL_MIN,sun_elevation_angle,earth_sun_distance):
    # http://landsathandbook.gsfc.nasa.gov/data_prod/prog_sect11_3.html
    # http://www.yale.edu/ceo/Documentation/Landsat_DN_to_Reflectance.pdf
    # define some needed parameters 
    ESUN_band= {'ESUN_1':1997.,'ESUN_2':1812.,'ESUN_3':1533.,'ESUN_4':1039.0,'ESUN_5':230.8,'ESUN_7':84.9}
    
    landsatImageSubsetRadiance = ((L_MAX -L_MIN)/(QCAL_MAX - QCAL_MIN)) * (landsatImageSubsetDN.astype(float) - QCAL_MIN) + L_MIN                          
                       
    # NOW REFLECTANCE
    
    landsatImageSubsetReflectance = (np.pi * landsatImageSubsetRadiance * (earth_sun_distance ** 2)) / (ESUN_band['ESUN_'+str(band)] * np.sin((np.pi /180.0) * sun_elevation_angle))
    
    landsatImageSubsetReflectance[landsatImageSubsetReflectance <= 0.0]= np.NaN
    
    # NOW DARK OBJECT SUBTRACT
    
    landsatImageSubsetReflectance= landsatImageSubsetReflectance - np.nanmin(landsatImageSubsetReflectance[np.nonzero(landsatImageSubsetReflectance)])
    
    
    print "band minimum is:", np.nanmin(landsatImageSubsetReflectance[np.nonzero(landsatImageSubsetReflectance)])

    return landsatImageSubsetReflectance

def ALItoReflectance(landsatImageSubsetDN,band,sun_elevation_angle,earth_sun_distance):
    ESUN_band= {'ESUN_PAN':1724.,'ESUN_2':1857.,'ESUN_3':1996.,'ESUN_4':1807.,'ESUN_5':1536,'ESUN_6':1145.,'ESUN_7':955.8,'ESUN_8':452.3,'ESUN_9':235.1,'ESUN_10':82.38}
    
    # TURNS INTO RADIANCE
    scalingFactor = {'2':.045, '3':.043,'4':.028,'5':.018,'6':.011,'7':.0091,'8':.0083,'9':.0028,'10':.00091}
    offset = {'2':-3.4, '3':-4.4,'4':-1.9,'5':-1.3,'6':-0.85,'7':-0.65,'8':-1.3,'9':-0.6,'10':-0.21}                
    
    landsatImageSubsetRadiance = (landsatImageSubsetDN.astype(float) * scalingFactor[str(band)]) + offset[str(band)]                       
           
    # NOW REFLECTANCE
    
    landsatImageSubsetReflectance = (np.pi * landsatImageSubsetRadiance * (earth_sun_distance ** 2)) / (ESUN_band['ESUN_'+str(band)] * np.sin((np.pi * sun_elevation_angle) /180.0))

    landsatImageSubsetReflectance[landsatImageSubsetReflectance < 0.0]= np.NaN

    # NOW DARK OBJECT SUBTRACT
    
    landsatImageSubsetReflectance= landsatImageSubsetReflectance - np.nanmin(landsatImageSubsetReflectance)

    return landsatImageSubsetReflectance   


# OPEN MULTIBAND LANDSAT 7 IMAGERY, (LEVEL 1G)
# TTHEN ANALIZE FOR BRAIDPLAIN AND PLUME CHARACTERISTICS  

# ----------------------------------------------------------
# TURN THE PLOTTING ROUTINE OFF / ON ON=1, OFF=0
# ----------------------------------------------------------

PLOT_AS_YOU_GO = 0

fjord = 'K'

sensorUsed = 'L7'
#'L7'
#'ALI'

import datetime
dateTimeNow = datetime.datetime.now()

dateTimeNowString = dateTimeNow.strftime('%Y%m%d_%H%M%S')
#dateTimeNow.strftime('%Y_%m_%d_%H_%M_%S')

saveAs = "Z:\\Documents\\DATA\\LANDSAT\\MAT_save\\"+fjord+"_Landsat_SSC_"+dateTimeNowString

if fjord  == 'N':
    
    
    rf_mask_path = 'Z:\\Documents\\DATA\\LANDSAT\\Nuuk\\riverCorridor\\NK_rf_mask.tif'

    # nuuk ROI
    utm_w = 520635.0
    utm_e = 581085.0  
    utm_s = 7090815.0
    utm_n = 7124085.0
    
elif fjord == 'K':
   
    rf_mask_path = 'Z:\\Documents\\DATA\\LANDSAT\\Kanger\\riverCorridor\\K_rf_mask_PolygonToRaster1.tif'

    # watson ROI
    # go into ENVI and find utm of upper left and lower right corners 
    utm_w = 492615.0
    utm_e = 546495.0 
    utm_s = 7419255.0
    utm_n = 7448085.0
    
   
    
elif fjord == 'P':
   
    rf_mask_path = 'Z:\\Documents\\DATA\\LANDSAT\\Ilulissat\\riverCorridor\\P_rf_mask.tif'
    
    
    # pakitsup ROI
    utm_w = 507645.0
    utm_e = 537495.0  
    utm_s = 7698555.0
    utm_n = 7714035.0
    
elif fjord == 'Q':
   
    rf_mask_path = 'Z:\\Documents\\DATA\\LANDSAT\\Ilulissat\\riverCorridor\\P_rf_mask.tif'
    
    
    # 
    utm_w = 484635.0
    utm_e = 536475.0  
    utm_s = 7609875.0
    utm_n = 7635975.0    
    
    
elif fjord == 'U':
    # Umiiviit River 
    rf_mask_path = 'Z:\\Documents\\DATA\\LANDSAT\\Kanger\\riverCorridor\\K_rf_mask_PolygonToRaster1.tif'
    
    #W, S   : E,N
    #479127,7367717 : 597277,7425285
    #ISSORTOQ 479182,7439887 : 550785,7474775
    # NEXT BIG NORTH 497544,7473862 : 558418,7503522
    # pakitsup ROI
    utm_w = 479127.0+500.0*30.0
    utm_e = 597277.0-1000.0*30.0 
    utm_s = 7367717.0 
    utm_n = 7425285.0
    
elif fjord == 'I':
    # Isortoq River 
    rf_mask_path = 'Z:\\Documents\\DATA\\LANDSAT\\Kanger\\riverCorridor\\K_rf_mask_PolygonToRaster1.tif'
    
    #W, S   : E,N
    #479127,7367717 : 597277,7425285
    #ISSORTOQ 479182,7439887 : 550785,7474775
    
    # pakitsup ROI
    utm_w = 479182.0
    utm_e = 550785.0 
    utm_s = 7439887.0 
    utm_n = 7474775.0

elif fjord == 'N2':
    #W, S   : E,N
    #N2 = next two north of isortoq
    # 497389,7474172 : 585720,7517170
    utm_w = 497389.0
    utm_e = 585720.0 
    utm_s = 7474172.0 
    utm_n = 7517170.0
    rf_mask_path = 'Z:\\Documents\\DATA\\LANDSAT\\Kanger\\riverCorridor\\K_rf_mask_PolygonToRaster1.tif'
# ----------------------------------------------------------
# define IR water threshold 
# ----------------------------------------------------------

# all values less than this are defined as water 
band5_thresh_define_as_water = .075
# kanger was .03 as of 3-10-2014
# was as low as .03 

# ----------------------------------------------------------
# Tell it what river/fjord corridor mask to use 
# ----------------------------------------------------------


rf_mask = arcpy.RasterToNumPyArray(rf_mask_path)
## river fjord gets spit out at 0, not river/fjord as 255
#
rf_mask = np.where(rf_mask == 0,1,0)

# ----------------------------------------------------------
# DEFINE THE ROI 
# ----------------------------------------------------------


# LOAD MASKS 
#landsatWater_B5_RIVER_BASE_clean = np.load('Z:\Documents\DATA\LANDSAT\Kanger\GOOD_2007_2012\landsatWater_B5_RIVER_BASE_clean.npy')
#LAND_WATER_MASK_Subset_BASE_clean = np.load('Z:\Documents\DATA\LANDSAT\Kanger\GOOD_2007_2012\LAND_WATER_MASK_Subset_BASE_clean.npy')

# 

#regional_x = np.round(((utm_e-utm_w)/30.),0)-1.0
#    
#regional_y = np.round(((utm_n-utm_s)/30.),0)+1.0

# the mask  sets the dimnesions that will be used to grid
regional_y = np.shape(rf_mask)[0]
regional_x = np.shape(rf_mask)[1]

# this is the spacing of the desired grid
xi = np.linspace(utm_w,utm_e,regional_x)
yi = np.linspace(utm_n,utm_s,regional_y)


# ----------------------------------------------------------
# NECESSARY STEPS TO ALLOW FOR LATER CONVERSION FROM DIGITAL NUMBERS TO RADIANCE TO REFLECTANCE 
# ----------------------------------------------------------

# READ Earth-Sun Distance table into PYTHON 

eS_path = 'Z:\\Documents\\DATA\\EARTH_SUN_DISTANCE\\Earth_Sun_DOY_dist.csv'
eS_dict = {}

with open(eS_path) as es:                  
    for line in es:
        #print line        
        
        es_split = line.split(",")
        eS_dict[int(es_split[0].strip())] = float(es_split[1].strip())



if sensorUsed == 'L7':
    if fjord  == 'N':
        
        fileDirectoryMaster = 'Z:\\Documents\\DATA\\LANDSAT\\Nuuk\\GOOD_2011_2013\\'
        fileDirectoryOUTPUT = 'Z:\\Documents\\DATA\\LANDSAT\\Nuuk\\pngOUT\\'
        fileDirectoryOUTPUT_Tiff = 'Z:\\Documents\\DATA\\LANDSAT\Nuuk\\tiffOUT\\'
        
        
    
    elif fjord == 'K':
        fileDirectoryMaster = 'Z:\\Documents\\DATA\\LANDSAT\\Kanger\\GOOD_2007_2012\\FULL_IMAGES\\'
        fileDirectoryOUTPUT = 'Z:\\Documents\\DATA\\LANDSAT\\Kanger\\pngOUT\\'
        fileDirectoryOUTPUT_Tiff = 'Z:\\Documents\\DATA\\LANDSAT\\Kanger\\tiffOUT\\'
        
        all_good_QC = [0,1,1,1,0,0,0,0,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,0,0,1,1,1,1,0,1,1]
        
    elif fjord == 'P':
        fileDirectoryMaster = 'Z:\\Documents\\DATA\\LANDSAT\\Ilulissat\\GOOD_2011_2012\\'
        fileDirectoryOUTPUT = 'Z:\\Documents\\DATA\\LANDSAT\\Ilulissat\\pngOUT\\'
        fileDirectoryOUTPUT_Tiff = 'Z:\\Documents\\DATA\\LANDSAT\\Ilulissat\\tiffOUT\\'
        
    elif fjord == 'Q':
        fileDirectoryMaster = 'Z:\\Documents\\DATA\\LANDSAT\\Qa\\GOOD\\'
        fileDirectoryOUTPUT = 'Z:\\Documents\\DATA\\LANDSAT\\Qa\\pngOUT\\'
        fileDirectoryOUTPUT_Tiff = 'Z:\\Documents\\DATA\\LANDSAT\\Qa\\tiffOUT\\'
    
    elif fjord == 'U':
        fileDirectoryMaster = 'Z:\\Documents\\DATA\\LANDSAT\\Kanger\\GOOD_2007_2012\\FULL_IMAGES\\'
        fileDirectoryOUTPUT = 'Z:\\Documents\\DATA\\LANDSAT\\Umiiviit\\pngOUT\\'
        fileDirectoryOUTPUT_Tiff = 'Z:\\Documents\\DATA\\LANDSAT\\Umiiviit\\tiffOUT\\'
        
    elif fjord == 'I':
        fileDirectoryMaster = 'Z:\\Documents\\DATA\\LANDSAT\\Kanger\\GOOD_2007_2012\\FULL_IMAGES\\'
        fileDirectoryOUTPUT = 'Z:\\Documents\\DATA\\LANDSAT\\Isortoq\\pngOUT\\'
        fileDirectoryOUTPUT_Tiff = 'Z:\\Documents\\DATA\\LANDSAT\\Isortoq\\tiffOUT\\'
        
    elif fjord == 'N2':
        fileDirectoryMaster = 'Z:\\Documents\\DATA\\LANDSAT\\Kanger\\GOOD_2007_2012\\FULL_IMAGES\\'
        fileDirectoryOUTPUT = 'Z:\\Documents\\DATA\\LANDSAT\\NorthTwo\\pngOUT\\'
        fileDirectoryOUTPUT_Tiff = 'Z:\\Documents\\DATA\\LANDSAT\\NorthTwo\\tiffOUT\\'
    
    # ----------------------------------------------------------
    # Count number of folders to be read 
    # ----------------------------------------------------------
    folderCount = 0 
    
    for n in os.listdir(fileDirectoryMaster): 
        if os.path.isdir(os.path.join(fileDirectoryMaster, n)):
            folderCount = folderCount + 1
            
    # ----------------------------------------------------------
    # pre allocate 
    # ----------------------------------------------------------
    
    # manual Qc
    manual_QC = np.empty(folderCount)
    manual_QC[:] = np.NAN 
    
    # Time 
    imageYear = np.zeros(folderCount)
    imageDOY  = np.zeros(folderCount)    
    imageTime  = []
    
    if fjord  == 'N':
        # river         
        meanRiverWidth_SAVE = np.zeros(folderCount)
        riverArea_SAVE = np.zeros(folderCount)
        
        meanRiverB4B3_SAVE = np.zeros(folderCount)    
        meanRiverB4B3_terminus_SAVE = np.zeros(folderCount)
        meanRiverB4B3_lake_SAVE = np.zeros(folderCount)
        
        meanRiverArea_SAVE = np.zeros(folderCount)    
        meanRiverArea_terminus_SAVE = np.zeros(folderCount)
        # doesn make sens to calc lake area as relatively unchanging 
    
    elif fjord == 'K':
        meanRiverB4B3_SAVE = np.zeros(folderCount)  
        meanRiverB4B3_above_gauge_SAVE = np.zeros(folderCount) 
        
        meanRiverB4_above_gauge_SAVE = np.zeros(folderCount) 
        meanRiverB3_above_gauge_SAVE = np.zeros(folderCount) 
        
        meanRiverB4B3_north_SAVE = np.zeros(folderCount)  
        meanRiverB4B3_south_SAVE = np.zeros(folderCount)  
        
        meanRiverB4B3_terminus_russell_SAVE = np.zeros(folderCount)
        meanRiverB4B3_terminus_leverett_SAVE = np.zeros(folderCount)
        meanRiverB3_terminus_leverett_SAVE = np.zeros(folderCount)
        meanRiverB4_terminus_leverett_SAVE = np.zeros(folderCount)
        meanRiverB4B3_terminus_south_SAVE = np.zeros(folderCount)
        
        meanRiverArea_below_gauge_SAVE = np.zeros(folderCount)  
        meanRiverArea_above_gauge_SAVE = np.zeros(folderCount) 
        meanRiverArea_north_SAVE = np.zeros(folderCount)  
        meanRiverArea_south_SAVE = np.zeros(folderCount)  
        
        meanRiverArea_terminus_russell_SAVE = np.zeros(folderCount)
        meanRiverArea_terminus_leverett_SAVE = np.zeros(folderCount)
        meanRiverArea_terminus_south_SAVE = np.zeros(folderCount)
        
        riverWidth = np.zeros(folderCount)
        riverArea = np.zeros(folderCount)
    
    elif fjord == 'P':
        meanRiverB4B3_SAVE = np.zeros(folderCount)    
        #meanRiverB4B3_terminus_SAVE = np.zeros(folderCount)
        meanRiverB4B3_lake1_SAVE = np.zeros(folderCount)
        meanRiverB4B3_lake2_SAVE = np.zeros(folderCount)
        meanRiverB4B3_lakeHydro_SAVE = np.zeros(folderCount)
        
        meanRiverArea_SAVE = np.zeros(folderCount)   
        
    elif fjord == 'Q':
        
        meanRiverB4B3_lake1_SAVE = np.zeros(folderCount)
        meanRiverB4B3_lake2_SAVE = np.zeros(folderCount)
        meanRiverB4B3_lake3_SAVE = np.zeros(folderCount)
       
    elif fjord == 'U':
        meanRiverB4B3_SAVE = np.zeros(folderCount)  
        meanRiverB4B3_upper_SAVE = np.zeros(folderCount)  
        meanRiverB4B3_above_lake_SAVE = np.zeros(folderCount)
        
        meanRiverArea_SAVE = np.zeros(folderCount)  
        meanRiverArea_upper_SAVE = np.zeros(folderCount)  
        meanRiverArea_above_lake_SAVE = np.zeros(folderCount)
    
    elif fjord == 'I':
        meanRiverB4B3_SAVE = np.zeros(folderCount)  
        meanRiverB4B3_terminus_SAVE = np.zeros(folderCount)  
        meanRiverB4B3_trib_SAVE = np.zeros(folderCount)
        
        meanRiverArea_SAVE = np.zeros(folderCount)  
        meanRiverArea_terminus_SAVE = np.zeros(folderCount)  
        meanRiverArea_trib_SAVE = np.zeros(folderCount)
    
    elif fjord == 'N2':
        meanRiverB4B3_SAVE = np.zeros(folderCount)  
        meanRiverB4B3_terminus_SAVE = np.zeros(folderCount)  
        meanRiverB4B3_trib_SAVE = np.zeros(folderCount)
        
        meanRiverArea_SAVE = np.zeros(folderCount)  
        meanRiverArea_terminus_SAVE = np.zeros(folderCount)  
        meanRiverArea_trib_SAVE = np.zeros(folderCount)
       
        
    # PRE ALLOCATE IMAGE CUBE   
    waterCube = np.zeros((np.shape(rf_mask)[0],np.shape(rf_mask)[1],folderCount))
    sedCube = np.zeros((np.shape(rf_mask)[0],np.shape(rf_mask)[1],folderCount))   
        
    ## plume 
    #plumeAreaPT4_SAVE = np.zeros(folderCount)
    #plumeLengthPT4_SAVE = np.zeros(folderCount)
    #plumeWidthPT4_SAVE = np.zeros(folderCount)     
    #plume_average_B4B3_SAVE = np.zeros(folderCount) 
    # ----------------------------------------------------------
    # MASTER FILE/FOLDER LEVEL 
    # ----------------------------------------------------------
    
    ticker = 0 
    
    # FOLDER LEVEL 
    for name in os.listdir(fileDirectoryMaster): 
        if name == '.DS_Store':
            continue  # skip the file
    
        # CHECK IF FOLDER IS DIRECTORY 
        # IF ITS A DIRECTORY YOU ARE GOING TO TREAT IT AS LANDSAT IMAGE- AND RUN ALL ADDITIONAL STEPS
        print "Before first if statement"
        if os.path.isdir(os.path.join(fileDirectoryMaster, name)):
            print name        
            fileDirectory = fileDirectoryMaster+name
    
            # LOOP THOOUGH ALL FILES IN THE DIRECTORy TO FIND THE MTL FILE
            for files in os.listdir(fileDirectory):
                #print files
                # need to get metadata out of text file first 
                if files.endswith('MTL.txt'):    
                    #print "STEP 1"
    
                    MTL_fileName = files
                    
            # ONCE CORRECT PATHS ARE ID'ed, start to use tome                 
            MTL_fileFullPath = fileDirectory+'\\'+MTL_fileName
            # create a empty dictionary that will be filled with input from the text file
            MLT_dictionary ={}
            
            with open(MTL_fileFullPath) as f:                  
                for line in f:
                    if "=" in line:
                        l = line.split("=")
                        MLT_dictionary[l[0].strip()] = l[1].strip()
                        
            # new naming convention as of 2012 re-name the metadata that you will need
            
            if MLT_dictionary.has_key('DATE_ACQUIRED') == 1:
                
                # if it has the new name, change it back to the old name 
                # OLD = NEW 
                MLT_dictionary['ACQUISITION_DATE']  = MLT_dictionary['DATE_ACQUIRED']
                
                MLT_dictionary['PRODUCT_SAMPLES_REF'] = MLT_dictionary['REFLECTIVE_SAMPLES']          
                MLT_dictionary['PRODUCT_LINES_REF']= MLT_dictionary['REFLECTIVE_LINES'] 
                
                MLT_dictionary['PRODUCT_LL_CORNER_MAPX'] = MLT_dictionary['CORNER_LL_PROJECTION_X_PRODUCT']                   
                MLT_dictionary['PRODUCT_LR_CORNER_MAPX'] = MLT_dictionary['CORNER_LR_PROJECTION_X_PRODUCT']     
                MLT_dictionary['PRODUCT_UL_CORNER_MAPY'] = MLT_dictionary['CORNER_UL_PROJECTION_Y_PRODUCT']   
                MLT_dictionary['PRODUCT_LL_CORNER_MAPY'] = MLT_dictionary['CORNER_LL_PROJECTION_Y_PRODUCT']
                
                MLT_dictionary['LMAX_BAND3'] = MLT_dictionary['RADIANCE_MAXIMUM_BAND_3']
                MLT_dictionary['LMIN_BAND3'] = MLT_dictionary['RADIANCE_MINIMUM_BAND_3']
                
                MLT_dictionary['LMAX_BAND4'] = MLT_dictionary['RADIANCE_MAXIMUM_BAND_4']
                MLT_dictionary['LMIN_BAND4'] = MLT_dictionary['RADIANCE_MINIMUM_BAND_4']
                
                MLT_dictionary['LMAX_BAND5'] = MLT_dictionary['RADIANCE_MAXIMUM_BAND_5']
                MLT_dictionary['LMIN_BAND5'] = MLT_dictionary['RADIANCE_MINIMUM_BAND_5']
                
                MLT_dictionary['QCALMAX_BAND3'] = MLT_dictionary['QUANTIZE_CAL_MAX_BAND_3']
                MLT_dictionary['QCALMIN_BAND3'] = MLT_dictionary['QUANTIZE_CAL_MIN_BAND_3']                    
                
                MLT_dictionary['QCALMAX_BAND4'] = MLT_dictionary['QUANTIZE_CAL_MAX_BAND_4']
                MLT_dictionary['QCALMIN_BAND4'] = MLT_dictionary['QUANTIZE_CAL_MIN_BAND_4']                    
                
                MLT_dictionary['QCALMAX_BAND5'] = MLT_dictionary['QUANTIZE_CAL_MAX_BAND_5']
                MLT_dictionary['QCALMIN_BAND5'] = MLT_dictionary['QUANTIZE_CAL_MIN_BAND_5']
                
                MLT_dictionary['SCENE_CENTER_SCAN_TIME'] = MLT_dictionary['SCENE_CENTER_TIME']
                
                # FROM METADATA PULL OUT THE FILENAMES YOU ARE GOING TO NEED TO USE 
                
                MLT_dictionary['BAND3_FILE_NAME'] = MLT_dictionary['FILE_NAME_BAND_3']
                MLT_dictionary['BAND4_FILE_NAME'] = MLT_dictionary['FILE_NAME_BAND_4']
                MLT_dictionary['BAND5_FILE_NAME'] = MLT_dictionary['FILE_NAME_BAND_5']
        
                # for checking UTM ZONE
                MLT_dictionary['ZONE_NUMBER'] = MLT_dictionary['UTM_ZONE']
        
        
            # ----------------------------------------------------------   
            #IF UTM ZONE IS NOT 22n do something about it 
            # ----------------------------------------------------------    
            
            
            # convert everything to utm22 n

            utm_zone_used = int(MLT_dictionary['ZONE_NUMBER'])
            print "The UTM Zone is:", utm_zone_used
            
            inputProjection = pyproj.Proj("+init=EPSG:32622")
            # taking advantage of the fact that EPSG codes are sequential like this for utm 
            outputProjection = pyproj.Proj("+init=EPSG:326"+MLT_dictionary['ZONE_NUMBER'])
            
            
#            MLT_dictionary['PRODUCT_SAMPLES_REF']
#            MLT_dictionary['PRODUCT_LINES_REF']
#            
#            warp_LLx, warp_LLy = pyproj.transform(inputProjection, outputProjection,float(MLT_dictionary['PRODUCT_LL_CORNER_MAPX']),float( MLT_dictionary['PRODUCT_LL_CORNER_MAPY']))       
#            warp_LRx, warp_ULy = pyproj.transform(inputProjection, outputProjection,float(MLT_dictionary['PRODUCT_LR_CORNER_MAPX'] ),float( MLT_dictionary['PRODUCT_UL_CORNER_MAPY'] ))       
#            
#            # redefining in the dictionary to reflect hte new UTM coords 
#            MLT_dictionary['PRODUCT_LL_CORNER_MAPX'] = str(warp_LLx)
#            MLT_dictionary['PRODUCT_LR_CORNER_MAPX'] = str(warp_LRx)
#            MLT_dictionary['PRODUCT_UL_CORNER_MAPY'] = str(warp_ULy)
#            MLT_dictionary['PRODUCT_LL_CORNER_MAPY'] = str(warp_LLy)
            
            print utm_w, utm_n, utm_e, utm_s
            utm_w, utm_n = pyproj.transform(inputProjection, outputProjection, utm_w,utm_n)
        
            utm_e, utm_s = pyproj.transform(inputProjection, outputProjection, utm_e,utm_s)
            print utm_w, utm_n, utm_e, utm_s
            # ----------------------------------------------------------   
            #NOW CALL OUT THE BAND FILE NAMES YOU WILL USE! 
            # ----------------------------------------------------------            
            # Band 3
            forFilePath_B3 = MLT_dictionary['BAND3_FILE_NAME']
            filePath_B3 = fileDirectory+'\\'+forFilePath_B3.split('"')[1]
            
            # band 4
            forFilePath_B4 = MLT_dictionary['BAND4_FILE_NAME']
            filePath_B4 = fileDirectory+'\\'+forFilePath_B4.split('"')[1]
            
            # band 5
            forFilePath_B5 = MLT_dictionary['BAND5_FILE_NAME']
            filePath_B5 = fileDirectory+'\\'+forFilePath_B5.split('"')[1]
            
            # ----------------------------------------------------------
            # CALC DAY OF YEAR                    
            # ----------------------------------------------------------
            
            ac_date_string =  MLT_dictionary['ACQUISITION_DATE']        
            time_split = ac_date_string.split('-')
            
            # date take YEAR MONTH DAY 
            tt = date(int(time_split[0]),int(time_split[1]),int(time_split[2]))
            
            DOY = tt.timetuple().tm_yday
            
            # calc / open day specific info    
            
            # dictionary was made above, feed it the day of the year as an int and it will return the earth sun distance in AU
            earth_sun_distance = eS_dict[DOY]          
    
            sun_elevation_angle = float(MLT_dictionary['SUN_ELEVATION'])    
            
            
           
             
    
            # ----------------------------------------------------------
            # OPEN BAND 5 NEAR IR 
            # PROCESS TO REFLECTANCE 
            # ----------------------------------------------------------
    
            landsatImageB5 = arcpy.RasterToNumPyArray(filePath_B5)
    
            # check dimensions are correct
            imageROW = np.shape(landsatImageB5)[0]
            imageCOL = np.shape(landsatImageB5)[1]
            productCOL = int(MLT_dictionary['PRODUCT_SAMPLES_REF'])
            productROW = int(MLT_dictionary['PRODUCT_LINES_REF'])
            
            if imageCOL == productCOL: 
                if imageROW  == productROW:
                    print "Band 5 processing"
                    x =np.linspace(float(MLT_dictionary['PRODUCT_LL_CORNER_MAPX']),float(MLT_dictionary['PRODUCT_LR_CORNER_MAPX']),num=imageCOL)
                    y =np.linspace(float(MLT_dictionary['PRODUCT_UL_CORNER_MAPY']),float(MLT_dictionary['PRODUCT_LL_CORNER_MAPY']),num=imageROW)
                    
                    # imshow is rows then columns 
                    # ----------------------------------------------------------
                    # FOR ROI
                    # ----------------------------------------------------------
                    xSTART = np.nonzero(x > utm_w)[0][0]
                    xSTOP = np.nonzero(x < utm_e)[0][-1]
                    ySTOP = np.nonzero(y < utm_s)[0][0]
                    ySTART = np.nonzero(y > utm_n)[0][-1]   
                            
                    # ----------------------------------------------------------
                    # BAND 5 
                    # ----------------------------------------------------------
                                                 
                    landsatImageSubsetDN_B5 = landsatImageB5[ySTART:ySTOP,xSTART:xSTOP]                      
                    
                    print ySTART
                    print ySTOP
                    print xSTART
                    print xSTOP
                    
                    #imshowwithZ(landsatImageSubsetDN_B5) 
                    
                    landsatImageSubsetRotatePLUME_B5 = L7toReflectance(landsatImageSubsetDN_B5,5,float(MLT_dictionary['LMAX_BAND5']),float(MLT_dictionary['LMIN_BAND5']),float(MLT_dictionary['QCALMAX_BAND5']),float(MLT_dictionary['QCALMIN_BAND5']),sun_elevation_angle,earth_sun_distance)                
                    
            
            # ----------------------------------------------------------
            # OPEN SEDIMENT PLUME BAND/ FJORD 
            # BAND 3 - RED 
            # ----------------------------------------------------------   
                
            landsatImageB3 = arcpy.RasterToNumPyArray(filePath_B3)
            
            # check dimensions are correct
            imageROW = np.shape(landsatImageB3)[0]
            imageCOL = np.shape(landsatImageB3)[1]
            productCOL = int(MLT_dictionary['PRODUCT_SAMPLES_REF'])
            productROW = int(MLT_dictionary['PRODUCT_LINES_REF'])
            
            if imageCOL == productCOL: 
                if imageROW  == productROW:
                    print "Band 3 processing"
                    x =np.linspace(float(MLT_dictionary['PRODUCT_LL_CORNER_MAPX']),float(MLT_dictionary['PRODUCT_LR_CORNER_MAPX']),num=imageCOL)
                    y =np.linspace(float(MLT_dictionary['PRODUCT_UL_CORNER_MAPY']),float(MLT_dictionary['PRODUCT_LL_CORNER_MAPY']),num=imageROW)
                    
                    # imshow is rows then columns 
                    # ----------------------------------------------------------
                    # FOR ROI
                    # ----------------------------------------------------------
                    xSTART = np.nonzero(x > utm_w)[0][0]
                    xSTOP = np.nonzero(x < utm_e)[0][-1]
                    ySTOP = np.nonzero(y < utm_s)[0][0]
                    ySTART = np.nonzero(y > utm_n)[0][-1]                   
                    
                    # ----------------------------------------------------------
                    # BAND 3 - same subset 
                    # ----------------------------------------------------------
                    landsatImageSubsetDN_B3 = landsatImageB3[ySTART:ySTOP,xSTART:xSTOP]                      
                    
                    landsatImageSubsetRotatePLUME_RED_not_masked = L7toReflectance(landsatImageSubsetDN_B3,3,float(MLT_dictionary['LMAX_BAND3']),float(MLT_dictionary['LMIN_BAND3']),float(MLT_dictionary['QCALMAX_BAND3']),float(MLT_dictionary['QCALMIN_BAND3']),sun_elevation_angle,earth_sun_distance)                
                    landsatImageSubsetRotatePLUME_RED = landsatImageSubsetRotatePLUME_RED_not_masked
                    
                
            
            # ----------------------------------------------------------
            # OPEN SEDIMENT PLUME BAND/ FJORD 
            # BAND 4 - NEAR IR (800 ish  nm)
            # ----------------------------------------------------------
    
                
            landsatImageB4 = arcpy.RasterToNumPyArray(filePath_B4)
            
            # check dimensions are correct
            imageROW = np.shape(landsatImageB4)[0]
            imageCOL = np.shape(landsatImageB4)[1]
            productCOL = int(MLT_dictionary['PRODUCT_SAMPLES_REF'])
            productROW = int(MLT_dictionary['PRODUCT_LINES_REF'])
            
            if imageCOL == productCOL: 
                if imageROW  == productROW:
                    print "Band 4 processing"
                    x =np.linspace(float(MLT_dictionary['PRODUCT_LL_CORNER_MAPX']),float(MLT_dictionary['PRODUCT_LR_CORNER_MAPX']),num=imageCOL)
                    y =np.linspace(float(MLT_dictionary['PRODUCT_UL_CORNER_MAPY']),float(MLT_dictionary['PRODUCT_LL_CORNER_MAPY']),num=imageROW)
                    
                    # imshow is rows then columns 
                    
                    xSTART = np.nonzero(x > utm_w)[0][0]
                    xSTOP = np.nonzero(x < utm_e)[0][-1]
                    ySTOP = np.nonzero(y < utm_s)[0][0]
                    ySTART = np.nonzero(y > utm_n)[0][-1]            
                    
                    landsatImageSubsetDN_B4 = landsatImageB4[ySTART:ySTOP,xSTART:xSTOP]   
                    
                    landsatImageSubsetRotatePLUME_NIR = L7toReflectance(landsatImageSubsetDN_B4,4,float(MLT_dictionary['LMAX_BAND4']),float(MLT_dictionary['LMIN_BAND4']),float(MLT_dictionary['QCALMAX_BAND4']),float(MLT_dictionary['QCALMIN_BAND4']),sun_elevation_angle,earth_sun_distance)                
                    
                    
    
            # END LANDSAT 7 SECTION
            # ----------------------------------------------------------


            
            # GRID INPUT DATA SO THAT IF YOU GET A LANDSAT IMAGE FROM A DIFFERENT ROW/PATH IT DOESN'T BLOW UP THE 
            # FOLLOWING NUMPY ARRAYS HAVE BEEN PRODUCED SO FAR
            # BAND 5        
            # landsatImageSubsetRotatePLUME_B5
            # BAND 4        
            # landsatImageSubsetRotatePLUME_NIR
            # BAND 3         
            # landsatImageSubsetRotatePLUME_RED
            # MAY NOT BE NECESSARY 
            # griddata needs a 2d grid to grid the imagery, that has then been flattened. to do this i am repeating the utm cooridinates xi, the correct number of times in the y direction, and vice versa
            rf_mask_TRICK_SHAPE = np.shape(landsatImageSubsetRotatePLUME_B5)
            print "WARNING: TRICK RF MASK SHAPE IS ON!"
            #if np.shape(rf_mask) !=  np.shape(landsatImageSubsetRotatePLUME_B5):
            if rf_mask_TRICK_SHAPE !=  np.shape(landsatImageSubsetRotatePLUME_B5):
                print "RE GRIDDING"
                xIMAGE = np.arange(x[xSTART],x[xSTOP],30)
                yIMAGE = np.arange(y[ySTART],y[ySTOP],-30)
                           
                
                x_subset_utm = np.tile(xIMAGE,np.shape(yIMAGE)[0])
                y_subset_utm = np.repeat(yIMAGE,np.shape(xIMAGE)[0])
                
               
                landsatImageSubsetRotatePLUME_B5 = griddata((x_subset_utm,y_subset_utm), landsatImageSubsetRotatePLUME_B5.flatten(), (xi[None,:], yi[:,None]), method='nearest')
                landsatImageSubsetRotatePLUME_NIR = griddata((x_subset_utm,y_subset_utm), landsatImageSubsetRotatePLUME_NIR.flatten(), (xi[None,:], yi[:,None]), method='nearest')
                landsatImageSubsetRotatePLUME_RED_not_masked = griddata((x_subset_utm,y_subset_utm), landsatImageSubsetRotatePLUME_RED.flatten(), (xi[None,:], yi[:,None]), method='nearest')
                landsatImageSubsetRotatePLUME_RED = landsatImageSubsetRotatePLUME_RED_not_masked
                
            # ----------------------------------------------------------
            # DO SEDIMENT CONCENTRATION WORK 
            print "STEP 7 "
            # ----------------------------------------------------------
            # MAKE FJORD LAND/WATER MASK 
            
            LAND_WATER_MASK = np.zeros(np.shape(landsatImageSubsetRotatePLUME_B5))
            
            # ALL BAND 5 REFLECATNCE VALUES LESS THAT .7 ARE DEFINED AS WATER (1.0)
            LAND_WATER_MASK[landsatImageSubsetRotatePLUME_B5 <= band5_thresh_define_as_water ] =1.0        
            
            B3_threshold_MASK = np.where(landsatImageSubsetRotatePLUME_RED < .05,0.0,1.0)        
            
            # HERE IS WHERE RELATIVE SEDIMENT CONCENTRATION IS CALCUATED 
            
            all_masks = LAND_WATER_MASK  * B3_threshold_MASK # * rf_mask
            print "RF Mask is off!"
            #all_masks = LAND_WATER_MASK  * B3_threshold_MASK
            
            # This erodes the area of pixels that is classified as WATER - helping to mask out pixels that are mixed between dry land and water
            #all_masks = scipy.ndimage.morphology.binary_erosion(all_masks,structure=np.ones((2,2)))               
            print " Erode Water is turned off!"
            
            combo = (np.float32(landsatImageSubsetRotatePLUME_NIR) / np.float32(landsatImageSubsetRotatePLUME_RED)) * all_masks             
            # use above line when rf mask has been created 
            #combo = (np.float32(landsatImageSubsetRotatePLUME_NIR) / np.float32(landsatImageSubsetRotatePLUME_RED)) * LAND_WATER_MASK * B3_threshold_MASK         
            
            combo[combo > 2 ] = np.NaN
            combo[combo < .00001 ] = np.NaN
 
 
            landsatImageSubsetRotatePLUME_NIR = np.float32(landsatImageSubsetRotatePLUME_NIR) * all_masks
            landsatImageSubsetRotatePLUME_NIR[landsatImageSubsetRotatePLUME_NIR < .00001 ] = np.NaN
            
            landsatImageSubsetRotatePLUME_RED = np.float32(landsatImageSubsetRotatePLUME_RED) * all_masks 
            landsatImageSubsetRotatePLUME_RED[landsatImageSubsetRotatePLUME_RED < .00001 ] = np.NaN
            # CALCULATE SED CON
            
            # This calcuated the pixels sediment concentration i mg/ l - see excel file 
            sed_con = 2227.5 * combo ** 7.37
            
 
 
 
        
            # ----------------------------------------------------------
            # NOW CALCUATE NUMBERS THAT WE CAN USE 
            print "STEP 8"
            # ----------------------------------------------------------
            
            if fjord  == 'N':       
                meanRiverB4B3_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[89:200,618:800]))
                meanRiverB4B3_terminus_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[511:566,896:1006]))
                meanRiverB4B3_lake_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[622:805,985:1302]))
      
                meanRiverArea_SAVE[ticker] = np.nansum(all_masks[89:200,618:800])
                meanRiverArea_terminus_SAVE[ticker] = np.nansum(all_masks[511:566,896:1006]) 
      
      
                               
                
            elif fjord  == 'K':
                
                # turbidity 
                meanRiverB4B3_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[533:665,475:700]))
                
                #imshowwithZ(combo) 
                
                
                meanRiverB4B3_above_gauge_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[494:545,701:722]))
                
                # Extract B4 and B3 also  
                meanRiverB4_above_gauge_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(landsatImageSubsetRotatePLUME_NIR[494:545,701:722]))
                meanRiverB3_above_gauge_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(landsatImageSubsetRotatePLUME_RED[494:545,701:722]))
                
                
                meanRiverB4B3_north_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[284:430,1062:1290])) 
                meanRiverB4B3_south_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[475:720,886:1544]))
        
                meanRiverB4B3_terminus_russell_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[225:262,1290:1366]))
                meanRiverB4B3_terminus_leverett_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[300:360,1372:1437]))
                meanRiverB3_terminus_leverett_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(landsatImageSubsetRotatePLUME_RED[300:360,1372:1437]))
                meanRiverB4_terminus_leverett_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(landsatImageSubsetRotatePLUME_NIR[300:360,1372:1437]))
                # this was combo[300:360,1372:1507] before 3-8-2014
                
                
                #imshowwithZ(combo)
                
                meanRiverB4B3_terminus_south_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[580:720,1545:1630]))
                # was combo[580:720,1545:1755] before 3-8-2014
                
                # river area in the GIVE
                meanRiverArea_below_gauge_SAVE[ticker] = np.nansum(all_masks[533:665,475:700])
                
                # rotate the river mouth subset -34 degrees 
                rotate_below = scipy.ndimage.interpolation.rotate(all_masks[533:765,375:700],-34.0)

                # trim fjord/tidal zone                 
                rbtrim = rotate_below[:,250:]
                #rbtrim = rotate_below[:,250:250+100]
                rbtrim_mask = np.zeros(np.shape(rbtrim))
                rbtrim_maskNODATA = np.zeros(np.shape(rbtrim))
                
                # rotating does an interpolation - so re make mask as binary
                rbtrim_mask[rbtrim > .50] = 1.0
                rbtrim_maskNODATA[rbtrim > .50] = 1.0
                
                #gdal.FillNodata(rbtrim_mask,rbtrim_maskNODATA)                
                
                #plt.imshow(rbtrim_mask)      
                #plt.show()
                ave_width = np.mean(np.nansum(rbtrim_mask,0))
                
                riverWidth[ticker] = ave_width
                
                riverArea[ticker] = np.nansum(np.nansum(rbtrim_mask,0))
                
                meanRiverArea_above_gauge_SAVE[ticker] = np.nansum(all_masks[494:545,701:722])
                meanRiverArea_north_SAVE[ticker] = np.nansum(all_masks[284:430,1062:1290])
                meanRiverArea_south_SAVE[ticker] = np.nansum(all_masks[475:720,886:1544])
        
                meanRiverArea_terminus_russell_SAVE[ticker] = np.nansum(all_masks[225:262,1290:1366])
                meanRiverArea_terminus_leverett_SAVE[ticker] = np.nansum(all_masks[300:360,1372:1437])
                meanRiverArea_terminus_south_SAVE[ticker] = np.nansum(all_masks[580:720,1545:1630])           
                
                
                        
                
            elif fjord  == 'P':
              
                meanRiverB4B3_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[312:415,456:531]))   
                
                meanRiverB4B3_lake1_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[349:357,591:600]))
                meanRiverB4B3_lake2_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[323:334,714:732]))
                meanRiverB4B3_lakeHydro_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[156:184,770:831]))
                
                meanRiverArea_SAVE[ticker] = np.nansum(all_masks[312:415,456:531])              
                
            elif fjord  == 'Q':
              
                
                
                meanRiverB4B3_lake1_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[254:266,493:512]))
                meanRiverB4B3_lake2_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[250:268,571:607]))
                meanRiverB4B3_lake3_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[271:288,718:730]))
            
            elif fjord  == 'U':
                # [Y1:Y2, X1:X2]
                
                meanRiverB4B3_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[257:471,500:1033]))
                meanRiverB4B3_upper_SAVE[ticker] =   scipy.stats.nanmean(scipy.stats.nanmean(combo[687:924,1504:1883]))
                meanRiverB4B3_above_lake_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[512:605,1875:1994])) 
                
                meanRiverArea_SAVE[ticker] =  np.nansum(all_masks[257:471,500:1033])    
                meanRiverArea_upper_SAVE[ticker] =   np.nansum(all_masks[687:924,1504:1883])    
                meanRiverArea_above_lake_SAVE[ticker] =   np.nansum(all_masks[512:605,1875:1994])           
                
            elif fjord  == 'I':
                # [Y1:Y2, X1:X2]
                
            
                meanRiverB4B3_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[232:638,3:285]))
                meanRiverB4B3_terminus_SAVE[ticker] =   scipy.stats.nanmean(scipy.stats.nanmean(combo[301:400,1857:2200]))
                meanRiverB4B3_trib_SAVE[ticker] = scipy.stats.nanmean(scipy.stats.nanmean(combo[651:854,1468:1630])) 
                
                meanRiverArea_SAVE[ticker] =  np.nansum(all_masks[232:638,3:285])    
                meanRiverArea_terminus_SAVE[ticker] =   np.nansum(all_masks[301:400,1857:2200])    
                meanRiverArea_trib_SAVE[ticker] =   np.nansum(all_masks[651:854,1468:1630])          
                
                
                # ADD ALSO 
                # add Big wide short South of Nuuk                
                #497548,7038171 : 522544,7050338
                #Godthabfjord 1 
                #531224,7182022 : 571894,7201820
                #Godthabfjord 2 
                # 545608,7140560 : 588361,7161371
            # ----------------------------------------------------------
            # Write data to a form that will last through different images 
            print "STEP 9 - Write"
            print "ticker=", ticker         
            # ----------------------------------------------------------        
        
            # general 
            #fileName.append(files2)        
            
            # time 
           
            imageYear[ticker] = int(time_split[0])
            imageDOY[ticker]  = DOY
            imageTime.append(MLT_dictionary['SCENE_CENTER_SCAN_TIME'])
            
                #imageTime.append(MLT_dictionary['SCENE_CENTER_SCAN_TIME'])
    
    #        # data CUBES    
    #        if all_good_QC[ticker] == 1:
    #            waterCube[:,:,ticker] = all_masks
    #            sedCube[:,:,ticker] = combo
    #        else:
    #            a = np.ones(np.shape(all_masks))
    #            a[a == 1] = np.NaN
    #            
    #            b = np.ones(np.shape(combo))
    #            b[b == 1] = np.NaN
    #            
    #            waterCube[:,:,ticker] = a
    #            sedCube[:,:,ticker] = b
                
    
            # ----------------------------------------------------------
            # PLOT 
            # ----------------------------------------------------------
            
            if PLOT_AS_YOU_GO == 1:
                # ----------------------------------------------------------
                # FIGURE 1         
                #plt.figure(1,figsize=(8.65,2.0),dpi=300)    
#                plt.figure(1)
#                        
#                plt.title(time_split[0]+'-'+str(DOY))
#                
#                sceneMean = scipy.stats.nanmean(scipy.stats.nanmean(combo))
#                imRiver = plt.imshow(combo) # combo[533:665,400:700]
#                #imshowwithZ(combo/sceneMean)
#                print np.nanmax(sed_con)
#                print np.nanmin(sed_con)
#                imRiver.set_cmap('spectral') # PuBu_r, 'BrBG_r'
#                #imRiver.set_clim(0,6000)
#                plt.colorbar() 
#                
#                plt.figure(2)
                #plt.plot(np.ravel(landsatImageSubsetRotatePLUME_NIR * all_masks),np.ravel(landsatImageSubsetRotatePLUME_RED * all_masks),'.')
                #plt.scatter(np.ravel(landsatImageSubsetRotatePLUME_NIR),np.ravel(landsatImageSubsetRotatePLUME_RED),s=20,c=np.ravel(sed_con))
                
#                plt.figure(3)
#                
                #imshowwithZ((combo/meanRiverB4B3_SAVE[ticker])-1.0)
                
                # dimension zero is north south, dim 1 is east west
                
                # VERY SLOW WAY TO DO IT
#                sed_con_smooth = np.zeros(np.shape(combo)) 
#                for j in xrange(0,np.shape(combo)[0]):
#                    
#                        for k in xrange(0,np.shape(combo)[1]):
#                            
#                            sed_con_smooth[j,k] = scipy.stats.nanmean(scipy.stats.nanmean(combo[j-1:j+2,k-1:k+2]))
                                
                  
                #tophat_2D_kernel = Tophat2DKernel(40)     
                #sed_con_smooth = convolve(sed_con,tophat_2D_kernel) # ,'fill',-999
                
                footprint=np.ones((2,2))   
                sed_con_smooth= scipy.ndimage.filters.generic_filter(sed_con,smoothingFunction,footprint=footprint)    
                #denom = 2227.5 * meanRiverB4B3_terminus_leverett_SAVE[ticker] ** 7.37               
                
                sed_con_smooth_masked = (sed_con_smooth) * all_masks #(sed_con_smooth/denom) * all_masks
                sed_con_smooth_masked[np.isnan(sed_con_smooth_masked)] = 0.0
                
                
                
                
                
                #masked_array = np.ma.array(sed_con_smooth_masked, mask=np.isnan(sed_con_smooth_masked))
                masked_array = np.ma.masked_where(sed_con_smooth_masked < 0.2, sed_con_smooth_masked)
                
                fig, ax = plt.subplots(figsize=(33,8.65*3),dpi=300)
#                \
                ax.imshow(landsatImageSubsetRotatePLUME_RED_not_masked, cmap=cm.gray,interpolation='none')
                #plt.hold(True)
                img = ax.imshow(masked_array, cmap=cm.jet,interpolation='none')                
                img.set_clim(0,10000)                
                #img.set_clim(.8,1.2)
                #plt.imshow(landsatImageSubsetRotatePLUME_RED, interpolation='nearest')
                #plt.hold(True)
                #plt.imshow(masked_array, interpolation='nearest', cmap=cmap)
                
                ax.set_title(str(int(imageYear[ticker]))+"-"+str(int(imageDOY[ticker])))

                # Add colorbar, make sure to specify tick locations to match desired ticklabels
                #cbar = fig.colorbar(ax)
                #plt.colorbar(img, ax=ax)
                #cbar.ax.set_yticklabels(['< -1', '0', '> 1'])# vertically oriented colorbar
                #plt.figure(2)

                #imshowwithZ(masked_array)                
                
                plt.show()
                #plt.savefig(fileDirectoryOUTPUT+str(int(imageYear[ticker]))+"-"+str(int(imageDOY[ticker]))+'N2_Colorbar10k.png')
                #imshowwithZ(sed_con_smooth_masked)
                
                
                
    #            plt.figure(2,figsize=(8.65,2.0),dpi=300)            
    #            plt.title(time_split[0]+'-'+str(DOY))
    #            
    #            imRiver = plt.imshow(landsatWATER_B5_RIVER)
    #            imRiver.set_cmap('BrBG_r') # PuBu_r
    #            imRiver.set_clim(0,1.0)
                            
                
                 
                
                
                       
                plt.close()
          
        

#            LL_corner_subset = arcpy.Point(utm_w,utm_s)
#            pixel_size = 30 # meters
#            
#            outRaster = arcpy.NumPyArrayToRaster(sed_con,LL_corner_subset,pixel_size,pixel_size)
#            outASCII = fileDirectoryOUTPUT_Tiff+name+".txt"
#            outGEOTIFF = fileDirectoryOUTPUT_Tiff+name+".tif"                
#            arcpy.RasterToASCII_conversion(outRaster, outASCII)
#    
#            # using GDAL TO TURN INTO GEOTIFF 
#            ds_in = gdal.Open(outASCII)
#            ds_out = drv.CreateCopy(outGEOTIFF, ds_in)
#            ds_in = None
#            ds_out = None
#            #os.sys('gdalwarp infile.tif outfile.tif -t_srs "+proj=longlat +ellps=WGS84"')
#
#            arcpy.DefineProjection_management(outGEOTIFF,sr)  
            #arcpy.DefineProjection_management()
            # add to the ticker to update 
            ticker+=1
        
        # OPTIONAL, TURN INTO GEOTIFF FOR ARC
        
        # units are in meters
        
#            LL_corner_subset = arcpy.Point(utm_w,utm_s)
#            pixel_size = 30 # meters
#            
#            outRaster = arcpy.NumPyArrayToRaster(landsatImageSubsetRotatePLUME_RED_not_masked,LL_corner_subset,pixel_size,pixel_size)
#            outASCII = fileDirectoryOUTPUT_Tiff+name+"_refl3.txt"
#            outGEOTIFF = fileDirectoryOUTPUT_Tiff+name+"_refl3.tif"                
#            arcpy.RasterToASCII_conversion(outRaster, outASCII)
#            
#            # using GDAL TO TURN INTO GEOTIFF 
#            ds_in = gdal.Open(outASCII)
#            ds_out = drv.CreateCopy(outGEOTIFF, ds_in)
#            ds_in = None
#            ds_out = None
#    #       
#            arcpy.DefineProjection_management(outGEOTIFF,sr)   
        
            
            #Clean up
            del landsatImageB3
            del landsatImageB4
            del landsatImageB5
            del combo
            del LAND_WATER_MASK

## aFTER ALL THE LOOPING 
#imageDOY_round = myround(imageDOY)
#tick2 = 1
#plt.figure(3)
#for iu in np.unique(imageYear):
#    
#    #rows =3.0                
#   #plt.subplot(rows,np.ceil(np.size( np.unique(imageYear))/rows),tick2)
#    # year was UI
#    plt.imshow(np.nansum(waterCube[:,:,imageYear == 2007],axis=2))
#    
#    plt.title(str(2007))
#    plt.colorbar()
#    
#    tick2+=1
##plt.show()

#plt.imshow(np.nansum(waterCube[:,:,imageYear == 2007],axis=2))

#data_out = (np.nansum(waterCube[:,:,imageYear == 2007],axis=2))

#
from scipy import io
    


if fjord  == 'N':     
    
    scipy.io.savemat(saveAs, mdict={'year':imageYear,
                                    'doy':imageDOY,
                                    'meanRiverB4B3_SAVE':meanRiverB4B3_SAVE,
                                    'meanRiverB4B3_terminus_SAVE':meanRiverB4B3_terminus_SAVE,
                                    'meanRiverB4B3_lake_SAVE':meanRiverB4B3_lake_SAVE,
                                    'meanRiverArea_SAVE':meanRiverArea_SAVE,
                                    'meanRiverArea_terminus_SAVE':meanRiverArea_terminus_SAVE})                         
                                    
elif fjord  == 'K':
    
    scipy.io.savemat(saveAs, mdict={'year':imageYear,
                                    'doy':imageDOY,
                                    'meanRiverB4B3_SAVE':meanRiverB4B3_SAVE,
                                    'meanRiverB4B3_above_gauge_SAVE':meanRiverB4B3_above_gauge_SAVE,
                                    'meanRiverB4_above_gauge_SAVE':meanRiverB4_above_gauge_SAVE,
                                    'meanRiverB3_above_gauge_SAVE':meanRiverB3_above_gauge_SAVE,
                                    'meanRiverB4B3_north_SAVE':meanRiverB4B3_north_SAVE,
                                    'meanRiverB4B3_south_SAVE':meanRiverB4B3_south_SAVE,
                                    'meanRiverB4B3_terminus_russell_SAVE':meanRiverB4B3_terminus_russell_SAVE,
                                    'meanRiverB4B3_terminus_leverett_SAVE':meanRiverB4B3_terminus_leverett_SAVE,
                                    'meanRiverB3_terminus_leverett_SAVE':meanRiverB3_terminus_leverett_SAVE,
                                    'meanRiverB4_terminus_leverett_SAVE':meanRiverB4_terminus_leverett_SAVE,
                                    'meanRiverB4B3_terminus_south_SAVE':meanRiverB4B3_terminus_south_SAVE,
                                    'meanRiverArea_below_gauge_SAVE':meanRiverArea_below_gauge_SAVE,
                                    'meanRiverArea_above_gauge_SAVE':meanRiverArea_above_gauge_SAVE,
                                    'meanRiverArea_north_SAVE':meanRiverArea_north_SAVE,
                                    'meanRiverArea_south_SAVE':meanRiverArea_south_SAVE,
                                    'meanRiverArea_terminus_russell_SAVE':meanRiverArea_terminus_russell_SAVE,
                                    'meanRiverArea_terminus_leverett_SAVE':meanRiverArea_terminus_leverett_SAVE,
                                    'meanRiverArea_terminus_south_SAVE':meanRiverArea_terminus_south_SAVE,
                                    'riverWidth':riverWidth,
                                    'riverArea':riverArea})
                                    
                                
elif fjord  == 'P':
    
    scipy.io.savemat(saveAs, mdict={'year':imageYear,
                                    'doy':imageDOY,
                                    'meanRiverB4B3_SAVE':meanRiverB4B3_SAVE,
                                    'meanRiverB4B3_lake1_SAVE':meanRiverB4B3_lake1_SAVE,
                                    'meanRiverB4B3_lake2_SAVE':meanRiverB4B3_lake2_SAVE,
                                    'meanRiverB4B3_lakeHydro_SAVE':meanRiverB4B3_lakeHydro_SAVE,
                                    'meanRiverArea_SAVE':meanRiverArea_SAVE})
                                    
elif fjord  == 'Q':
    
    scipy.io.savemat(saveAs, mdict={'year':imageYear,
                                    'doy':imageDOY,
                                    'meanRiverB4B3_lake1_SAVE':meanRiverB4B3_lake1_SAVE,
                                    'meanRiverB4B3_lake2_SAVE':meanRiverB4B3_lake2_SAVE,
                                    'meanRiverB4B3_lake3_SAVE':meanRiverB4B3_lake2_SAVE })
                                    
elif fjord  == 'U':
    
    scipy.io.savemat(saveAs, mdict={'year':imageYear,
                                    'doy':imageDOY,
                                    'meanRiverB4B3_SAVE':meanRiverB4B3_SAVE,
                                    'meanRiverB4B3_upper_SAVE':meanRiverB4B3_upper_SAVE,
                                    'meanRiverB4B3_above_lake_SAVE':meanRiverB4B3_above_lake_SAVE,
                                    'meanRiverArea_SAVE':meanRiverArea_SAVE,
                                    'meanRiverArea_upper_SAVE':meanRiverArea_upper_SAVE,
                                    'meanRiverArea_above_lake_SAVE':meanRiverArea_above_lake_SAVE})
                                    
elif fjord  == 'I':
    
    scipy.io.savemat(saveAs, mdict={'year':imageYear,
                                    'doy':imageDOY,
                                    'meanRiverB4B3_SAVE':meanRiverB4B3_SAVE,
                                    'meanRiverB4B3_terminus_SAVE':meanRiverB4B3_terminus_SAVE,
                                    'meanRiverB4B3_trib_SAVE':meanRiverB4B3_trib_SAVE,
                                    'meanRiverArea_SAVE':meanRiverArea_SAVE,
                                    'meanRiverArea_terminus_SAVE':meanRiverArea_terminus_SAVE,
                                    'meanRiverArea_trib_SAVE':meanRiverArea_trib_SAVE})
                                                                      

                        
else:
    print "no fjord saved"

