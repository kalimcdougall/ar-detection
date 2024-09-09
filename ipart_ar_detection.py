# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 16:34:00 2022

@author: Kali
"""

import os, sys
import pandas as pd
import numpy as np

from slugify import slugify #Removes invalid characters from output file name 

from ipart.utils import funcs
from ipart import thr
from ipart.utils import plot
from ipart.AR_detector import findARs
from ipart.AR_tracer import readCSVRecord, trackARs, filterTracks

import matplotlib.pyplot as plt
import cartopy.crs as ccrs

#------------------------------------------------------------------------------
# Load in eivt & nivt files

uflux = funcs.readNC('C:/Users/Kali/Desktop/School/Masters/Data/ipart_ar_detection/e-ivt-1995-1.nc', 'p71.162')
vflux = funcs.readNC('C:/Users/Kali/Desktop/School/Masters/Data/ipart_ar_detection/n-ivt-1995-1.nc', 'p72.162')

#------------------------------------------------------------------------------
# Create ivt output file

dir = 'C:/Users/Kali/Desktop/School/Masters/Data/ipart_ar_detection/outputs/1995/'

output_file = dir + 'ivt-1995-1.nc' # Create output file for ivt calculation

#------------------------------------------------------------------------------
# Compute & plot IVT

ivtdata = np.ma.sqrt(uflux.data**2 + vflux.data**2) #Calculate ivt

print(np.percentile(ivtdata, 98)) #Find value of 98th ivt percentile

ivtNV = funcs.NCVAR(ivtdata, 'ivt', uflux.axislist, {'name': 'ivt', 'long_name': 'integrated vapor transport (IVT)',
                                            'standard_name': 'integrated_vapor_transport',
                                            'title': 'integrated vapor transport (IVT)',
                                            'units': getattr(uflux, 'units', '')}) #Create attributes for ivt file
print('\n# Saving output to:\n', output_file)
funcs.saveNC(output_file, ivtNV, 'w') #Write output file

# Plot uflux, vflux, & ivt

figure = plt.figure(figsize=(7,10),dpi=100)
idx = 95  #Select time step from the beginning
time_str = uflux.getTime()[idx]

plot_vars = [uflux.data[idx], vflux.data[idx], ivtdata.data[idx]]
titles = ['U', 'V', 'IVT']

for ii, vii in enumerate(plot_vars):
    axii = figure.add_subplot(3, 1, ii+1, projection=ccrs.PlateCarree())
    iso = plot.Isofill(vii, 10, 1, 2)
    plot.plot2(vii, iso, axii,
            title='%s, %s' %(str(time_str), titles[ii]),
            xarray=uflux.getLongitude(),
               yarray=uflux.getLatitude(),
            legend='local',
            fix_aspect=False)

figure.show()

#------------------------------------------------------------------------------
# Compute & plot THR

# Slice data by region - only if you are using file larger than desired region
LAT1 = -40; LAT2 = -80;
SHIFT_LON = 60 #Shift longitude to center AOI

var = ivtNV.sliceData(LAT1, LAT2, axis=2)
var = ivtNV.shiftLon(SHIFT_LON)

# Optional orographic & high terrain info used to enhance continent-penetration of landfalling ARs
#oroNV = None #Orographic parameter - requires file w surface terrain/elevation info
#HIGH_TERRAIN = 2000 #Surface height (m) above which is defined as high terrain

KERNEL = [16,18,18] #Specifies temporal & spatial scales - time steps, number of grids

ivt, ivt_rec, ivt_ano = thr.THR(var, KERNEL) #Perform THR

# Save THR results to .nc file
ivt_file = dir + 'ivt-1995-1.nc'
fname = os.path.split(ivt_file)[1]
file_out_name = '%s-THR-kernel-t%d-s%d.nc'\
    %(os.path.splitext(fname)[0], KERNEL[0], KERNEL[1])

abpath_out = os.path.join(dir, file_out_name)
print('\n### <testrotating filter>: Saving output to:\n', abpath_out)
funcs.saveNC(abpath_out, ivt, 'w')
funcs.saveNC(abpath_out, ivt_rec, 'a')
funcs.saveNC(abpath_out, ivt_ano, 'a')

# Plot THR results

figure = plt.figure(figsize=(7,10),dpi=100)
idx = 95  #Select time step from the beginning
time_str = ivt.getTime()[idx]

plot_vars = [ivt.data[idx], ivt_rec.data[idx], ivt_ano.data[idx]]
iso = plot.Isofill(plot_vars, 12, 1, 1, min_level=0, qr=0.001)
titles=['IVT (=IVT_rec + IVT_ano)', 'IVT_rec', 'IVT_ano']

for ii, vii in enumerate(plot_vars):
    axii=figure.add_subplot(3, 1, ii+1, projection=ccrs.PlateCarree())
    
    plot.plot2(vii, iso, axii,
            title='%s, %s' %(str(time_str), titles[ii]),
               xarray=ivt.getLongitude(),
               yarray=ivt.getLatitude(),
            legend='local',
            fix_aspect=False)

figure.show()

#------------------------------------------------------------------------------
# Detect ARs 

YEAR = 1995
TIME_START = '1995-01-20 00:00:00' #AR began on Jan 24 & ended Jan 25, iceberg calved Jan 25
TIME_END = '1995-01-30 00:00:00'

thr_file = dir + 'ivt-1995-1-THR-kernel-t16-s18.nc'

PLOT = True #Plot figures of ivt with detected ARs
SHIFT_LON = 60

PARAM_DICT={
    'thres_low': 500, #kg/m*s, define AR candidates as regions >= this anomalous ivt
    'min_area': 50*1e4, #km^2, drop AR candidates smaller than this area
    'max_area': 1000*1e4, #km^2, drop AR candidates larger than this area
    'min_LW': 2, #float, minimal length/width ratio
    'min_lat': -40, #degree, exclude systems whose centroids are lower than this latitude
    'max_lat': -80, #degree, exclude systems whose centroids are higher than this latitude
    'min_length': 2000, #km, ARs shorter than this length treated as relaxed
    'min_length_hard': 800, #km, ARs shorter than this length are discarded
    'rdp_thres': 2, #degree lat/lon, error when simplifying axis using rdp algorithm
    'fill_radius': None, #grids, remove small holes in AR contour
    'single_dome': False, #do peak parition or not, used to separated systems merged together with an outer contour
    'max_ph_ratio': 0.6, #max prominence/height ratio of local peak, only used when single_dome=True
    'edge_eps': 0.4 #minimal proportion of flux component in direction to total flux 
    }

ivt = funcs.readNC(thr_file, 'ivt')
ivt_rec = funcs.readNC(thr_file, 'ivt_rec')
ivt_ano = funcs.readNC(thr_file, 'ivt_ano')

# Shift longitude - for Larsen A
uflux = uflux.shiftLon(SHIFT_LON)
vflux = vflux.shiftLon(SHIFT_LON)
ivt = ivt.shiftLon(SHIFT_LON)
ivt_rec = ivt_rec.shiftLon(SHIFT_LON)
ivt_ano = ivt_ano.shiftLon(SHIFT_LON)

# Slice data by year - only if you are using file larger than desired time frame
uflux = uflux.sliceData(TIME_START, TIME_END, axis=0).squeeze()
vflux = vflux.sliceData(TIME_START, TIME_END, axis=0).squeeze()
ivt = ivt.sliceData(TIME_START, TIME_END, axis=0).squeeze()
ivt_rec = ivt_rec.sliceData(TIME_START, TIME_END, axis=0). squeeze()
ivt_ano = ivt_ano.sliceData(TIME_START, TIME_END, axis=0).squeeze()

# Get coordinates
latax = uflux.getLatitude()
lonax = uflux.getLongitude()
timeax = ivt.getTime()
timeax = ['%d-%02d-%02d %02d:00' %(timett.year, timett.month, timett.day, timett.hour) for timett in timeax]

# Call detection function
time_idx, labels, angles, crossfluxes, result_df = findARs(ivt.data, ivt_rec.data,
            ivt_ano.data, uflux.data, vflux.data, latax, lonax, times=timeax, **PARAM_DICT)

# Print number of ARs detected
print('Number of ARs found during %s - %s = %d.' %(TIME_START, TIME_END, len(result_df)))
result_df.head(4)

# Plot detected ARs at single timestep

plot_idx = time_idx[14] #Define timestep

plot_time = timeax[plot_idx]
slab = ivt.data[plot_idx]
slab_rec = ivt_rec.data[plot_idx]
slab_ano = ivt_ano.data[plot_idx]
ardf = result_df[result_df.time==plot_time]

plot_vars = [slab, slab_rec, slab_ano]
titles = ['IVT', 'THR_recon', 'THR_ano']
iso = plot.Isofill(plot_vars, 12, 1, 1, min_level=0, max_level=800)

figure = plt.figure(figsize=(12,10), dpi=100)

for jj in range(len(plot_vars)):
    ax = figure.add_subplot(3, 1, jj+1, projection=ccrs.PlateCarree())
    pobj = plot.plot2(plot_vars[jj], iso, ax,
            xarray=lonax, yarray=latax,                    
            title='%s %s' %(plot_time, titles[jj]),
            fix_aspect=False)
    plot.plotAR(ardf, ax, lonax)
    
# =============================================================================
# # Plot detected ARs at all timesteps
# 
# label_timeax = funcs.num2dateWrapper(labels.getTime(),
#                 labels.axislist[0].units)
#     
# for (ii, timett) in zip(time_idx, label_timeax):
# 
#             timett_str = str(timett)
# 
#             slab = ivt.data[ii]
#             slab_rec = ivt_rec.data[ii]
#             slab_ano = ivt_ano.data[ii]
#             ardf = result_df[result_df.time==timett]
# 
#             plot_vars = [slab, slab_rec, slab_ano]
#             titles = ['IVT', 'THR_recon', 'THR_ano']
#             iso = plot.Isofill(plot_vars, 12, 1, 1, min_level=0, max_level=800)
# 
#             figure = plt.figure(figsize=(12,10), dpi=100)
# 
#             for jj in range(len(plot_vars)):
#                 ax = figure.add_subplot(3, 1, jj+1,
#                             projection=ccrs.PlateCarree())
#                 pobj = plot.plot2(plot_vars[jj], iso, ax,
#                         xarray=lonax, yarray=latax,
#                         title='%s %s' %(timett_str, titles[jj]),
#                         fix_aspect=False)
# 
#                 plot.plotAR(ardf, ax, lonax)        
# =============================================================================
    
# Save AR records to .csv file
abpath_out = os.path.join(dir, 'ar_records_1995.csv')
print('\n# Saving output to :\n', abpath_out)

if sys.version_info.major==2:
    np.set_printoptions(threshold=np.inf)
elif sys.version_info.major==3:
    np.set_printoptions(threshold=sys.maxsize)
result_df.to_csv(abpath_out, index=False)

# Save AR labels to .nc file
abpath_out = os.path.join(dir, 'ar_s_18_%d_labels.nc' %YEAR)
print('\n# Saving output to:\n', abpath_out)
funcs.saveNC(abpath_out, labels, 'w')
funcs.saveNC(abpath_out, angles, 'a')
funcs.saveNC(abpath_out, crossfluxes, 'a')

#------------------------------------------------------------------------------
# Track ARs

print('\n# Read in file:\n', 'ar_records_1995.csv')
ardf = readCSVRecord(dir + 'ar_records_1995.csv')

#PLOT = True #Plot track movements
SCHEMATIC = True #Plot linkage schematic
TIME_GAP_ALLOW = 6 #Int, hours, gap allowed to link 2 records - should be time resolution of data
TRACK_SCHEME = 'simple' #Tracking scheme, 'simple': tracks are simple paths, 'full': network scheme, tracks connected by joint points
MAX_DIST_ALLOW = 1200 #Int, max distance in km to define neighbourhood relationship
MIN_DURATION = 24 #Int, min duration in hrs to keep a track
MIN_NONRELAX = 1 #Int, min number of non-relaxed records in a track in order to keep a track

# Call tracking function
track_list = trackARs(ardf, TIME_GAP_ALLOW, MAX_DIST_ALLOW,
                      track_scheme=TRACK_SCHEME, isplot=SCHEMATIC, plot_dir=dir)
print('Number of AR tracks=', len(track_list))

# Optional filtering of detected AR tracks
track_list = filterTracks(track_list, MIN_DURATION, MIN_NONRELAX)

# Plot single AR sequence

latax = np.arange(-80,-30) #Latitude domain to plot
lonax = np.arange(200,320) #Longitude domain to plot

plot_ar = track_list[3] #Select AR index to track

figure = plt.figure(figsize=(12,6), dpi=100)
ax = figure.add_subplot(111, projection=ccrs.PlateCarree())
plot.plotARTrack(plot_ar, latax, lonax, ax, full=True)
plt.title('AR[3] 1995-01-23/25')

# Save output

for ii in range(len(track_list)):
    tii = track_list[ii]
    trackidii = '%d%d' %(tii.data.loc[0,'time'].year, ii+1)
    tii.data.loc[:,'trackid']=trackidii
    tii.trackid = trackidii

    if ii==0:
        trackdf = tii.data
    else:
        trackdf = pd.concat([trackdf,tii.data],ignore_index=True)


    # Plot all AR sequences
    figure = plt.figure(figsize=(12,6), dpi=100)
    ax = figure.add_subplot(111, projection=ccrs.PlateCarree())
    plot.plotARTrack(tii, latax, lonax, ax, full=True)

    # Save plots
    plot_save_name = 'ar_track_%s' %trackidii
    plot_save_name = os.path.join(dir, plot_save_name) #May need slugify()
    print('\n# : Save figure to', plot_save_name)
    figure.savefig(plot_save_name+'.png', dpi=100, bbox_inches='tight')

    plt.close(figure)


# Save AR tracks to .csv file
abpath_out = os.path.join(dir,'ar_tracks_1995.csv')
print('\n# Saving output to:\n',abpath_out)
if sys.version_info.major==2:
    np.set_printoptions(threshold=np.inf)
elif sys.version_info.major==3:
    np.set_printoptions(threshold=sys.maxsize)
trackdf.to_csv(abpath_out,index=False)