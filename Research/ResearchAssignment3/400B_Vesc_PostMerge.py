#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 11:33:23 2025

@author: aaronhernandez
"""

'''
This research projects main goal is to determine the escape velocity of the 
merged M31 and Mw galaxy as a function of radius (R). 
The first plot I create will be the point mass assumption for the merged galaxy
at snap number 650
'''



# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import the GalaxyMass to determine the component masses of the merged
# M31 and MW galaxies
from GalaxyMass_spydr import Component_mass
#Import the Read function to read the files of the MW and M31
from ReadFile_spydr import Read

#Convert G to prefered units
G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)

'''
This Part of the code is for The point mass Approximation
'''
def merged_mass():
    '''
    Function to sum up the mass from the MW and M31 and return the mass
    in units of solar masses
    
    Returns
    -------
    Total_Mass : [float] Total mass of the merged MW and M31 galaxies

    '''
    #Calculate the total mass of the MW
    MW_halo = Component_mass("MW_650.txt",1)
    MW_disk = Component_mass("MW_650.txt",2)
    MW_bulge = Component_mass("MW_650.txt",3)
    
    MW_Mass = MW_bulge + MW_disk + MW_halo
    
    #Calculate the total mass of M31
    M31_halo = Component_mass("M31_650.txt",1)
    M31_disk = Component_mass("M31_650.txt",2)
    M31_bulge = Component_mass("M31_650.txt",3)
    
    M31_Mass = M31_bulge + M31_disk + M31_halo
    
    #Sum to get the total combined mass
    Total_Mass = (MW_Mass + M31_Mass)*1e12
    # returns the total mass of M31 and MW
    return Total_Mass

def point_mass_Vesc():
   '''
    Function to calculate the point mass escape velocity of the merged 
    M31 and MW galaxies

    Returns
    -------
    R : [array] Radii
    Vesc_pm : [array] Escape velocities at every radius from R

    '''
   #Define a radius array 
   R = np.linspace(0.0001,50,300)
   #calc the velocity ar each radius value from the radius array
   Vesc_pm = np.sqrt(2*G*merged_mass()/R)
   #returns the radius array and escape velocity array
   return R,Vesc_pm

'''
This Part of the code is for The Hernquist Approximation
'''
def Hernquist(h_a):
    '''
    Function to calculate the Hernquist mass and then return a radii array and 
    the hernquist escape velocity array
    
            MAY NEED TO CALCULATE NEW h_a FOR MERGED GALAXIES
    
    Parameters
    ----------
    h_a : [int] Hernquist scale radius

    Returns
    -------
    R_hern : [array] Radius out 
    Vesc_hern : [array] Escape Velcoity at every radius in R_hern

    '''
    #Define a new radius array, it is the same as R from befroe but 
    #called something else for simplicity
    R_hern = np.linspace(0.0001,50,300)
    
    # Calc Galaxy mass compnents except for the halo
    MW_disk = Component_mass("MW_650.txt",2)
    MW_bulge = Component_mass("MW_650.txt",3)
    
    M31_disk = Component_mass("M31_650.txt",2)
    M31_bulge = Component_mass("M31_650.txt",3)
    # Sum the components
    components_xHalo = (MW_disk + MW_bulge + M31_bulge + M31_disk)*1e12
    # Calc the Halo masses
    MW_halo = Component_mass("MW_650.txt",1)
    M31_halo = Component_mass("M31_650.txt",1)
    #Sum the halo masses of M31 and the MW
    M_halo = MW_halo + M31_halo
    # Lab 4 code
    a = M_halo*1e12
    b = R_hern**2/(h_a+R_hern)**2
    # Sum the components with the Hernquist halo
    M_hern =  a*b + components_xHalo#Herquist profile
   
    #Vesc with Hernquist mass taken into account
    Vesc_hern = np.sqrt(2*G*M_hern/R_hern)
    #returns the radius array and escape velocity array
    return R_hern,Vesc_hern




'''
Plots
'''
#Point Mass
plt.figure()
r,V_pm = point_mass_Vesc()
plt.plot(r,V_pm,'r')
plt.xlabel('Radius [Kpc]')
plt.ylabel('Escape Velocity [Km/S]')
plt.title('Escape Velocity assuming a Point Mass Approximation')
plt.show()

#Hernquist 
plt.figure()
R_Hernquist,V_Hernquist = Hernquist(60)
plt.plot(R_Hernquist,V_Hernquist,'k')
plt.xlabel('Radius [Kpc]')
plt.ylabel('Escape Velocity with Hernquist Adjustment [Km/S]')
plt.title('Escape Velocity with Hernquist Halo Correction')
plt.show()

# Log log plot of the two escape velocities
plt.figure()
plt.loglog(V_pm,V_Hernquist,'b')
plt.xlabel('Point Mass V_esc [Km/s]')
plt.ylabel('Hernquist V_esc [Km/S]')
plt.title('Point Mass VS Hernquist')
plt.show()

    
