#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 00:12:10 2025

@author: aaronhernandez
"""

'''
This research projects main goal is to determine the escape velocity of the 
merged M31 and Mw galaxy as a function of radius (R). 
The first plot I create will be the point mass assumption for the merged galaxy
at snap number 650

OutLine: 
    1) Point mass aprox function
    2) Hernquist function
    3) M33 speed determination 
    4) Function to feternine if M33 is boaunded or not
    

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
#Import the Mass Profile class to find a new Hernquist scale radius
from MassProfile import MassProfile
#Import the Center of Mass class to get COM velocity for M33
from CenterOfMass2 import CenterOfMass

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
def Hernquist(r, a, mhalo):
    ''' This method returns the mass enclosed within a radius based on
    the analytic Hernquist density profile.

    PARAMETERS
    ----------
    r : `float` 
        radius to compute mass within in kpc
    a : `float`
        Hernquist profile scale radius in kpc
    mhalo : `astropy.Quantity`
        total halo mass in Msun

    RETURNS
    -------
    m_enc : `astropy.Quantity'
        mass enclosed by r in Msun
    '''

    # adding units
    r = r*u.kpc
    a = a*u.kpc
    
    # compute numerator and denominator separately
    A = mhalo * r**2
    B = (a + r)**2
    
    return A/B



def ScaleRadius(a):
    ''' 
    This function will creat a plot of the halo mass distrubution 
    of the merged M31 and MW galaxies and will allow me to find the best
    scale radius (a) that fits the distribution

    PARAMETERS
    ----------
    a : [float] Hernquist profile scale radius in kpc
    
    RETURNS
    -------
    A plot of the mass distribution and the Hernquist approx
    '''
    
    # read in galaxy information
    mwProf = MassProfile('MW', 650) # mass profile
    m31Prof = MassProfile('M31', 650)
   
     
    MW_haloM = Component_mass('MW_650.txt', 1) * 1e12 * u.Msun 
    M31_haloM = Component_mass('M31_650.txt', 1) * 1e12 * u.Msun 
        # halo mass in Msun
    M_halo_tot = M31_haloM + MW_haloM  
        
    # radius array in kpc
    r_arr = np.linspace(0.1, 50, 100)
    r_arrMW = np.linspace(0.1, 50, 100)
    r_arrM31 = np.linspace(0.1, 50, 100)
    
    # calculate mass profiles
    MW_halo = mwProf.massEnclosed(1, r_arrMW)*1e12
    M31_halo = m31Prof.massEnclosed(1, r_arrM31)*1e12
    #Concatenate mass profile arrays and radius arrays
    Halo_combo = np.concatenate((MW_halo,M31_halo), axis=0)
    Rs_array = np.concatenate((r_arrM31,r_arrMW), axis=0)

    # make plot
    fig, ax = plt.subplots()
    # lines
    ax.plot(r_arr, Hernquist(r_arr, a, M_halo_tot), 
            c='r', label='Analytic Halo, a={} kpc'.format(a))
    ax.plot(Rs_array, Halo_combo, c='b', linestyle=':', label='Halo')
    
    
ScaleRadius(140)    
# 140 looks good but somthing isnt right with the plot 
# Need to go to office hours to sort this out



def Hernquist_approx(h_a):
    '''
    Function to calculate the Hernquist mass and then return a radii array and 
    the hernquist escape velocity array
    
            MAY NEED TO CALCULATE NEW h_a FOR MERGED GALAXIES
          $$$ Redo HW 5, with the Dark Matter particles of both galaxies 
          concatenated into one array and plot that with the hernquist profile
          both against R
    
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
    MW_disk = Component_mass("MW_650.txt",2) #Milkey Way
    MW_bulge = Component_mass("MW_650.txt",3)
    
    M31_disk = Component_mass("M31_650.txt",2) #M31
    M31_bulge = Component_mass("M31_650.txt",3)
    # Sum the Disk and bulge components seperatly
    components_bulge = (MW_bulge + M31_bulge )*1e12
    components_disk = (MW_disk + M31_disk)*1e12
    # Calc the Halo masses
    MW_halo = Component_mass("MW_650.txt",1)
    M31_halo = Component_mass("M31_650.txt",1)
    #Sum the halo masses of M31 and the MW
    M_halo = MW_halo + M31_halo
    
    # Lab 4 code
    a = M_halo*1e12
    b = R_hern**2/(h_a+R_hern)**2
    #Herquist profile
    M_hern =  a*b
   
    #Calulate the poentials of the form |P| = |GM/R|
    Halo_pot = G*M_hern/(R_hern + h_a)
    Bulge_pot = G*components_bulge/R_hern
    Disk_pot = G*components_disk/R_hern
    
    #Sum the poentials
    potential = Halo_pot + Bulge_pot + Disk_pot
    #Vesc with Hernquist mass taken into account
    Vesc_hern = np.sqrt(2*potential)
    #returns the radius array and escape velocity array
    return R_hern,Vesc_hern


def M33COM():
    '''
    This function will determine the position and velocity of the 
    center of mass of M33 using the Dark matter particles
    
    Parameters
    ----------
    N/A

    Returns
    -------
    M33posMAG : [float] Magnitude of the position vector
    M33velMAG : [float] Magnitude of the velocity vector
    
    '''
    COM = CenterOfMass('M33_650.txt',1) #Initalizw class
    
    M33pos = COM.COM_P(0.1,2) #calc COM position
    
    x = M33pos[0] #for visual sake I made each component of the position
    y = M33pos[1] # vector a seperate variable so the code to get the 
    z = M33pos[2] # magnitude and the variables inserted into the functions
    # would be easier to read
    
    M33vel = COM.COM_V(x,y,z) #calc COM velocity 
    
    vx = M33vel[0] # Same proccess as mentioned above
    vy = M33vel[1]
    vz = M33vel[2]
    
    M33velMAG = np.sqrt(vx**2 + vy**2 + vz**2) #Velocity vector magnitude
    M33posMAG = np.sqrt(x**2 + y**2 + z**2) #Position vector magnitude
    
    return M33posMAG,M33velMAG 





# Treat M33 as a point mass, what is the gernal velocity compared
# to the remnent


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

    