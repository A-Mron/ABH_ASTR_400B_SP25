#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 18:53:37 2025

@author: aaronhernandez
"""
# import modules
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.table as tbl

#import neede functions
from ReadFile_spydr import Read
from CenterOfMass_spydr import CenterOfMass
from astropy.constants import G
'''
----------PLEASE READ-----------
I ran into an error I could not fix and that means 
I cant plot the curves. But I am uploading my code so you can see 
what I did do and hopefully I get credit for the functions I did write
---------------------------------
'''

class MassProfile:
    
    def __init__(self, gal_name,snap_num):
        '''
        Class to deturnmine the mass distributions of each galaxy at
        a given snap number 
        Parameters
        ----------
        gal_name : Name of the galaxy you want:
            - MW
            - M31
            - M33
        snap_num : The Snap number of the galaxy
            - Example: 0,1,2,3,4
        '''
        # add a string of the filenumber to the value “000”
        ilbl = '000' + str(snap_num)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        self.filename="%s_"%(gal_name) + ilbl + '.txt'
        
        # read data in the given file using Read
        self.time, self.total, self.data = Read(self.filename)                                                                                             
    
        G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        #create an array to store indexes of particles of desired Ptype                                
        #self.index = np.where(self.data['type'] == ptype)
        
        self.x = self.data['x']*u.kpc
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc
        self.m = self.data['m']
        self.gname = self.filename 
    
    def MassEnclosed(self,ptype,r):
        '''
        Parameters
        ----------
        ptype : [integer] particle type we want
                1 - Halo
                2 - Disk
                3 - Bulge
        r : array of radii values

        Returns
        -------
        masses: array of masses that are within a set of radii

        '''
        #create an array to store indexes of particles of desired Ptype 
        index = np.where(self.data['type'] == ptype)
        #call the COM class so we can access the functions inside
        com = CenterOfMass(self.filename,2)
        #array of zeros that the masses will be added to
        masses = np.zeros(len(r))
        galcom = com.COM_P(0.1)*u.kpc
        
        #shift to the COM refrence frame 
        xpos = self.x[index] - galcom[0]
        ypos = self.y[index] - galcom[1]
        zpos = self.z[index] - galcom[2]
        
        xpos=xpos*u.kpc
        ypos=ypos*u.kpc
        zpos=zpos*u.kpc
        
        particle_mass = self.m[index]
    
        
        #for look will iterate over all the r values in the array
        for i in range(len(r)):
            r_measure = r[i]*u.kpc  #the measured radius
            
            #calculate the COM R
            r_part = np.sqrt(xpos**2+ypos**2+zpos**2)*u.kpc
            
            index2 = np.where(r_measure<r_part)
            part_enclosed = particle_mass[index2]
            
            masses += np.sum(part_enclosed)
        masses = masses*(1e10*u.solMass)
    
        return masses
    

    def MassEnclosedTotal(self,r):
        '''
        Parameters
        ----------
        r : raray of radiai values
            
        Returns
        -------
        mass_total : array of the total enclosed mass of the 
        halo + disk + bulge of the desired galaxy
        '''
        #if statment will ensure M33 calculation does not have the bulge 
        #accounted for 
        if self.gname == "M33":
            #calc the mass in the halo within every r value
            halomass = self.MassEnclosed(1,r)
            #calc the mass in the disk within every r value
            diskmass = self.MassEnclosed(2,r)
            mass_total = halomass + diskmass
        else:
            #calc the mass in the halo within every r value
            halomass = self.MassEnclosed(1,r)
            #calc the mass in the disk within every r value
            diskmass = self.MassEnclosed(2,r)
            #calc the mass in the bulge within every r value
            bulgemass = self.MassEnclosed(3,r)
            #adds the three arrays together to get the total mass
            mass_total = halomass + diskmass + bulgemass
        
        return mass_total
    
    def HernquistMass(self,radius,a,Mhalo):
        '''
        Parameters
        ----------
        radius : [float] the set radius you desire 
        a : [float] scale factor 
        Mhalo : [float] mass of the halo

        Returns
        -------
        m_halo : halo mass in units of Solar Masses
        '''
        
        m_hern = (Mhalo*radius**2)/((a+radius)**2)*u.solMass
        
        #rho = (m_halo*a/(2*np.pi*radius))*(1/((a+radius)**3))
        
        return m_hern
    
    def CircularVelocity(self,ptype,r):
        '''
        Parameters
        ----------
        ptype : [integer] particle type we want 
                1 - Halo
                2 - Disk
                3 - Bulge
        r : array of r values

        Returns
        -------
        V_circ : array of cirvular velocities
        '''
        #Import and Convert G to the units we desire
        #from astropy.constants import G
        #G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        
        M = self.MassEnclosed(ptype, r) #get array of masses
        
        r = r*u.kpc #addd usnits to r array
        
        vel_array = np.zeros(len(r)) #initialize velocity array
        #for look will iterate over all the r values in the array
        for j in range(len(r)):
            #circular velocity equation
            V_circ = round(np.sqrt(self.G*M[j]/r[j]),2)
            #add to - array
            vel_array += V_circ
       
        return vel_array #circ vel array
    
    def CircularVelocityTotal(self,r):
        '''
        Parameters
        ----------
        r : array of r values

        Returns
        -------
        V_total : array of the final velocities 
        '''
        #Import and Convert G to the units we desire
        #from astropy.constants import G
        #G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        
        #initalize 0 array
        vel_circ_total = np.zeros(len(r))
        #add the units to r for calculations
        r = r*u.kpc
        #calc the total mass enclosed by the radii
        M = self.MassEnclosedTotal(r)
        
        for k in range(len(r)):
           Vels = np.sqrt(self.G*M[k]/r[k])
           vel_circ_total += Vels
       
        return vel_circ_total
    
    
    def HernquistVCirc(self,radius,a,Mhalo):
        '''
        Parameters
        ----------
        radius : distance in kpc 
        a : scale factor
        Mhalo : TYPE
            DESCRIPTION.

        Returns
        -------
        V_circ_Hern : Hernquist Circular speed in [km/s]
        '''
        #import G again and convert units again
        #from astropy.constants import G
        #G.to(u.kpc*u.km**2/u.s**2/u.Msun)
       
        M_H = self.HernquistMass(radius, a, Mhalo)
        
        #Calcs Hernquist velocty
        V_circ_Hern = np.sqrt(self.G*M_H/radius)
        
        return V_circ_Hern
    
    
'''
PLOTS
'''

if __name__ == '__main__' : 
    
    #initaluze our array of distances from 0.1kpc to 30.2 kpc
    #units are added in the finctions so I dont need it here
    r = np.arange(0.1,30.2,1.5)
    
    #initalize the three mass profiles for the three galaxies 
    MW = MassProfile('MW',0)
    M31 = MassProfile('M31',0)
    M33 = MassProfile('M33',0)
    
    '''
    Mass Profiles
    '''
    #MW components
    MW_halo = MW.MassEnclosed(1,r)
    MW_disk = MW.MassEnclosed(2,r)
    MW_bulge = MW.MassEnclosed(3,r)
    MW_total = MW.MassEnclosedTotal(r)
    
    #M31 Components
    M31_halo = M31.MassEnclosed(1,r)
    M31_disk = M31.MassEnclosed(2,r)
    M31_bulge = M31.MassEnclosed(3,r)
    M31_total = M31.MassEnclosedTotal(r)
    
    #M33 components
    M33_halo = M33.MassEnclosed(1,r)
    M33_disk = M33.MassEnclosed(2,r)
    M33_bulge = M33.MassEnclosed(3,r)
    M33_total = M33.MassEnclosedTotal(r)
    
    
   
    fig2,ax2 = plt.subplots()
    plt.plot(r,MW_halo,color='r',label='MW Halo')
    plt.plot(r,MW_disk,color='g',label='MW Disk')
    plt.plot(r,MW_bulge,color='k',label='MW Bulge')
    plt.plot(r,MW_total,color='b',label='MW Total')
    #plot the Hernquist mass
    plt.plot(r,MW.HernquistMass(30,25,MW.MassEnclosed(1,30)),color='pink'
             ,label='Hernquist')
    
    plt.plot(r,M31_halo,color='r',label='M31 Halo')
    plt.plot(r,M31_disk,color='g',label='M31 Disk')
    plt.plot(r,M31_bulge,color='k',label='M31 Bulge')
    plt.plot(r,M31_total,color='b',label='M31 Total')
    #plot the Hernquist mass
    plt.plot(r,M31.HernquistMass(30,25,M31.MassEnclosed(1,30)),color='pink'
             ,label='Hernquist')
    
    plt.plot(r,M33_halo,color='r',label='M33 Halo')
    plt.plot(r,M33_disk,color='g',label='M33 Disk')
    plt.plot(r,M33_total,color='k',label='M33 Total')
    #plot the Hernquist mass
    plt.plot(r,M33.HernquistMass(30,25,M33.MassEnclosed(1,30)),color='pink'
             ,label='Hernquist')
   
    plt.xlabel('Radius Out [kpc]')
    plt.ylabel('Mass [10^12 SolMass]')
    plt.grid()
    
    '''
    Rotation Curve
    '''
    # MW velocity compoents
    MW_halo_vel = MW.CircularVelocity(1,30)
    MW_disk_vel = MW.CircularVelocity(2,30)
    MW_bulge_vel = MW.CircularVelocity(3,30)
    # MW Hernquist velocity
    MW_total_vel = MW.HernquistVCirc(30, 25,MW.MassEnclosed(1,30))
    
    # M31 velocity compoents
    M31_halo_vel = M31.CircularVelocity(1,30)
    M31_disk_vel = M31.CircularVelocity(2,30)
    M31_bulge_vel = M31.CircularVelocity(3,30)
    # M31 Hernquist velocity
    M31_total_vel = M31.HernquistVCirc(30, 25,M31.MassEnclosed(1,30))
   
    # M33 velocity compoents
    M33_halo_vel = M33.CircularVelocity(1,30)
    M33_disk_vel = M33.CircularVelocity(2,30)
    # M33 Hernquist velocity
    M33_total_vel = M33.HernquistVCirc(30, 25,M33.MassEnclosed(1,30))
    
    fig3,ax3 = plt.subplots()
    # Milky Way Rotation curve
    plt.plot(r,MW_halo_vel,color='r',label='MW Halo')
    plt.plot(r,MW_disk_vel,color='g',label='MW Disk')
    plt.plot(r,MW_bulge_vel,color='k',label='MW Bulge')
    plt.plot(r,MW_total_vel,color='b',label='MW Total')
    #plot the Hernquist mass
    plt.plot(r,MW.HernquistVCirc(30,25,MW.MassEnclosed(1,30)),color='pink'
             ,label='Hernquist Velocity')
    
    # M31 Rotation curve
    plt.plot(r,M31_halo_vel,color='r',label='M31 Halo')
    plt.plot(r,M31_disk_vel,color='g',label='M31 Disk')
    plt.plot(r,M31_bulge_vel,color='k',label='M31 Bulge')
    plt.plot(r,M31_total_vel,color='b',label='M31 Total')
    #plot the Hernquist mass
    plt.plot(r,M31.HernquistVCirc(30,25,M31.MassEnclosed(1,30)),color='pink'
             ,label='Hernquist Velocity')
    
    # M33 Rotation curve 
    plt.plot(r,M33_halo_vel,color='r',label='M33 Halo')
    plt.plot(r,M33_disk_vel,color='g',label='M33 Disk')
    plt.plot(r,M33_total_vel,color='b',label='M33 Total')
    #plot the Hernquist mass
    plt.plot(r,M33.HernquistVCirc(30,25,MW.MassEnclosed(1,30)),color='pink'
             ,label='Hernquist Velocity')

'''
Due to the error I could not get the plots out
but my first guess at a was a=25
'''

    
    

    
    
        
    
            