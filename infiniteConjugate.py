# -*- coding: utf-8 -*-
"""
INFINITE CONJUGATE IMAGING
Ray tracing
Created on Thu Nov 1 2018

@author: Carmel Howe and Peter Quicke

https://en.wikipedia.org/wiki/Ray_transfer_matrix_analysis
https://docs.sympy.org/latest/modules/physics/optics/gaussopt.html

Code notes:
Assumes infinite lens diameters so max angle is set using numerical aperture of
objective lens.
Sensor Parameters can be easily changed along with the displacement in z and x (height) of your
point source

"""

import numpy as np
from sympy.physics.optics import RayTransferMatrix, ThinLens, FreeSpace, GeometricRay
import matplotlib.pyplot as plt
import os, sys

#figure font stuff
font = {'family': 'sans',
        'weight': 'normal',
        'size': 14,
        }

plt.rc('font',**font)

#define displacement in z and location of image
init_z = 0#1.5e-6 
init_h = 0#2.5e-6 

#system parameters
fOb=7.2e-3 # focal length objective
fT=180e-3 # focal length tube
NA = 1
n = 1.333 # refractive index of medium
theta = np.arcsin(NA/n) #Theta max for obj lens  
objRad = np.tan(theta)*fOb # objective radius

NAtube = fOb/fT
thetaTube = np.arcsin(NAtube/1) #Theta max for obj lens  
tubeRad = np.tan(thetaTube)*fT # objective radius

## Sensor ##
pixel=6.5e-6 # pixel size 
numPixels = 2048
sensorSize = pixel*numPixels 

rayDensity = 1000
angle=np.linspace(-theta,theta,rayDensity) # max angles
color_idx = np.linspace(0, 1, len(angle)) # for colormap

pixelBins = np.arange(-sensorSize/25,sensorSize/25,pixel) # 25 is a random number

#lens distances from origin in mm
distTube = ((fOb)+(fOb+fT))*1000
distSensor = ((fOb)+(fOb+fT)+fT)*1000

def origin(init_h,init_theta,init_z):
    ray = GeometricRay(init_h,init_theta)
    prop_mat = FreeSpace(fOb+init_z)
    out_ray = prop_mat*ray
    return out_ray
    
def propagate_to_tube(init_h,init_theta,init_z):
    return FreeSpace(fT+fOb)*ThinLens(fOb)*origin(init_h,init_theta,init_z)

def propagate_to_sensor(init_h,init_theta,init_z):
    return FreeSpace(fT)*ThinLens(fT)*propagate_to_tube(init_h,init_theta,init_z)


############################################################
    ############ FIGURES ############
############################################################

############ new directory for figures ############
saveVar = input("Do you want to save the figures? (Y/N):    ")

if saveVar == 'Y':
    newpath = ((r'H:\Python_Scripts\LightFieldRayTracing\Figures\RayDensity%s_z%s_h%s') %(len(angle),init_z,init_h)) 
    if not os.path.exists(newpath): os.makedirs(newpath)

############ plot origin ############
plt.plot([-init_z*2,init_z*2],[0,0],'dimgrey')
plt.plot([init_z],[init_h],'.r')
plt.plot(0,0,'.k')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#if saveVar == 'Y':
   # plt.savefig("Origin_z%s_h%s.png" %(init_z,init_h),format='png', dpi=1000)
   
plt.show()

############ origin to objective ############
#fig = plt.gcf()
#ax = plt.subplot(111)
#for a, i in zip(angle, color_idx):
#    ax.plot([init_z*1000,fOb*1000],[init_h,origin(init_h,a,init_z).height],
#            color=plt.cm.rainbow(i))
#    
#plt.plot([fOb*1000,fOb*1000],[-objRad,objRad],'darkgrey')    
#
##axis/size formatting
#fig.set_size_inches(12,4)
#ax.spines['left'].set_visible(False)
#ax.spines['right'].set_visible(False)
#ax.spines['top'].set_visible(False)
#ax.spines['bottom'].set_linewidth(1)
#ax.axes.get_yaxis().set_ticks([])
##ax.axes.get_xaxis().set_ticks([])
#plt.xlim((init_z,(fOb/20)*1000))
#plt.ylim(-0.0004,0.0004)
#plt.xlabel('Distance (mm)', fontdict = font)
#plt.tight_layout()
#if saveVar == 'Y':
    #plt.savefig("originToObj_z%s_h%s.png" %(init_z,init_h),format='png', dpi=1000)

#plt.show()


############# Plot Tube to Sensor #################
#fig = plt.gcf()
#ax = plt.subplot(111)
#
#for a, i in zip(angle, color_idx):
#    ax.plot([distTube,distMLA],[propagate_to_tube(init_h,a,init_z).height,
#            propagate_to_sensor(init_h,a,init_z).height],color=plt.cm.rainbow(i))
#
#ax.plot([distMLA,distMLA],[-MLArad,MLArad],'darkgrey')
#    
##axis/size formatting
#fig.set_size_inches(12,4)
#ax.spines['left'].set_visible(False)
#ax.spines['right'].set_visible(False)
#ax.spines['top'].set_visible(False)
#ax.spines['bottom'].set_linewidth(1)
#ax.axes.get_yaxis().set_ticks([])
##ax.axes.get_xaxis().set_ticks([])
#plt.xlabel('Distance (mm)', fontdict = font)
#plt.tight_layout()
#plt.xlim(((distMLA)-1,(distMLA)+0.25))
#plt.ylim((-MLArad,MLArad))  
#if saveVar == 'Y':
    #plt.savefig('tubeToMLA.png', format='png', dpi=1000)

#plt.show()
#

############ Plot Tube to sensor ############
#fig = plt.gcf()
#ax = plt.subplot(111)
##plt.plot([0,1],[0,0],'dimgrey')   # optical axis
#
#for a, i in zip(angle, color_idx):
#    ax.plot([distTube,distSensor],[propagate_to_tube(init_h,a,init_z).height,
#            propagate_to_sensor(init_h,a,init_z).height],color=plt.cm.rainbow(i))
#    
#ax.plot([distTube,distTube],[-tubeRad,tubeRad],'darkgrey') # Tube
#ax.plot([distSensor,distSensor],[-(sensorSize/2),(sensorSize/2)])
#
##axis/size formatting
#fig.set_size_inches(12,4)
#ax.spines['left'].set_visible(False)
#ax.spines['right'].set_visible(False)
#ax.spines['top'].set_visible(False)
#ax.spines['bottom'].set_linewidth(1)
#ax.axes.get_yaxis().set_ticks([])
##ax.axes.get_xaxis().set_ticks([])
#plt.xlim((distTube,distSensor))
#plt.xlabel('Distance (mm)', fontdict = font)
#plt.tight_layout()
#if saveVar == 'Y':
    #plt.savefig('MLAToSensor.png', format='png', dpi=1000)

#plt.show()


########## Plot origin to sensor ############
fig = plt.gcf()
ax = plt.subplot(111)

for a, i in zip(angle, color_idx):
    ax.plot([init_z*1000,fOb*1000],[init_h,origin(init_h,a,init_z).height],color=plt.cm.rainbow(i))    
    ax.plot([fOb*1000,distTube],[origin(init_h,a,init_z).height,propagate_to_tube(init_h,a,init_z).height],color=plt.cm.rainbow(i))       
    ax.plot([distTube,distSensor],[propagate_to_tube(init_h,a,init_z).height,propagate_to_sensor(init_h,a,init_z).height],color=plt.cm.rainbow(i))
    
ax.plot([fOb*1000,fOb*1000],[-objRad,objRad],'darkgrey') #objective
#ax.plot([(fOb*1000)+fT,(fOb*1000)+fT],[-objRad,objRad],'--') #BFP
ax.plot([distTube,distTube],[-tubeRad,tubeRad],'darkgrey') # Tube
ax.plot([distSensor,distSensor],[-(sensorSize/2),(sensorSize/2)])
    
fig.set_size_inches(12,4)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_linewidth(1)
ax.axes.get_yaxis().set_ticks([])
#ax.axes.get_xaxis().set_ticks([])
plt.xlabel('Distance (mm)', fontdict = font)
plt.tight_layout()
if saveVar == 'Y':
    plt.savefig("ic_OriginToSensor_z%s_h%s.png" %(init_z,init_h),format='png', dpi=1000)
    
plt.show()


############ Histogram of rays on each pixel of sensor ############
sensorLoc = np.array([propagate_to_sensor(init_h,a,init_z).height for a in angle],dtype=float)

np.histogram(sensorLoc,bins=pixelBins)

fig = plt.gcf()
ax = plt.subplot(111)

ax.hist(sensorLoc,bins=pixelBins)

#axis/size formatting
fig.set_size_inches(12,4)
ax.spines['left'].set_visible(1)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_linewidth(1)
plt.xlabel('Distance Across Sensor', fontdict = font)
plt.ylabel('Count', fontdict = font)
plt.xlim(-2e-4,2e-4)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.tight_layout()
if saveVar == 'Y':
    plt.savefig("ic_Histogram_z%s_h%s.png" %(init_z,init_h),format='png', dpi=1000)
    
plt.show()

    

