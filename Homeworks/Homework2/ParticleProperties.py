#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 14:33:08 2025

@author: aaronhernandez
"""
import numpy as np #needeed modules
import astropy.units as u
from ReadFile import Read

def ParticleInfo(filename,particle_typ,particle_number):
   '''
    This function will read a text file and provide information of the 
    particles position,veloicty and mass based on the users inputs.
    '''
   time,particle_num,data = Read(filename) #set valiables 
    #to the returned values of the Read function
    
   index = np.where(data['type'] == particle_typ) #index will separate the 
    #various particle types in the data
    
   x = data['x'][index] #catesian coordinates pulled from data array
   y = data['y'][index]
   z = data['z'][index]
    
   dist = np.sqrt(x[particle_number - 1]**2 + y[particle_number - 1]**2 + z[particle_number - 1]**2)*u.kpc #Calculates mag of distance in kpc
   dist = np.round(dist,3) #rounds the particles distence to 3 decimal places

   Vx = data['vx'][index] #components of velocity 
   Vy = data['vy'][index]
   Vz = data['vz'][index]
   velo = np.sqrt(Vx[particle_number - 1]**2 + Vy[particle_number - 1]**2+ Vz[particle_number - 1]**2)*u.km/u.s #calculates mag of speed in km/s
   velo = np.round(velo,3) #round the particle velocity to 3 decimal places
   
   
   M = data['m'][index]*(10**10)*u.M_sun #mass in 10^10 solar masses
   Mass = M[particle_number - 1] #retreves the particle's mass
   
   dist_ly = dist.to(u.lyr) #Comverts the distence from kpc to Lyr
   dist_ly = np.round(dist_ly,3) #rounds the distence in lyr to 3 decimals
   print('The 3D Distence is: ', dist) #print out statements for the desired values
   print('The 3D Velocity is: ', velo)
   print('The Mass is ', Mass)
   print('The 3D Distence in Lightyears is: ', dist_ly)
   return dist,velo,Mass #return the distance,velocity, and mass of the
                         # desired particle
'''
   print('The 3D Distence is: ', dist)
   print('The 3D Velocity is: ', velo)
   print('The Mass is ', Mass)
   print('The 3D Distence in Lightyears is: ', dist_ly)
'''
#This will give the distence,velocity,and mass of the 100th particle in the 
# disk region 
ParticleInfo('MW_000_test.txt',2,99)
