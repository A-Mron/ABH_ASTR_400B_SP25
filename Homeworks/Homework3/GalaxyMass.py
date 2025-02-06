#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 13:39:58 2025

@author: aaronhernandez
"""
 #needeed modules
import numpy as np 
import astropy.units as u
from ReadFile import Read
#pandas is needed for the cration of the data table
import pandas as pd

def Component_mass(filename,particle_typ):
    '''
    This code will return the total mass of any desired galaxy component
    
    Parameters
    ----------
    filename : File the data comes from
    particle_typ : 1 - Halo
                   2 - Disk
                   3 - Bulge

    Returns
    -------
    M_tot : Total mass of galaxy component

    '''
    time,particle_num,data = Read(filename) #set valiables 
     #to the returned values of the Read function
     
    index = np.where(data['type'] == particle_typ) #index will separate the 
     #various particle types in the data
    
    #Set total mass initially equal to 0
    M_tot = 0

    # Sum over all the masses of the given index into the M_tot variable
    # multiplied by 10^10 to abound for the units the masses are in 
    M_tot = np.sum(data['m'][index])*(10**10)
   
    
    print(M_tot)
    return M_tot

#Calculate the masses of the three components of the MW galaxy
# rounded to 3 decimals
MW_halo = np.round(Component_mass('MW_000_test.txt',1)/(10**12),3)
MW_disk = np.round(Component_mass('MW_000_test.txt',2)/(10**12),3)
MW_bulge = np.round(Component_mass('MW_000_test.txt',3)/(10**12),3)

#Calculate the masses of the three components of the M31 galaxy
# rounded to 3 decimals
M31_halo = np.round(Component_mass('M31_000_test.txt',1)/(10**12),3)
M31_disk= np.round(Component_mass('M31_000_test.txt',2)/(10**12),3)
M31_bulge = np.round(Component_mass('M31_000_test.txt',3)/(10**12),3)

#Calculate the masses of the two components of the M33 galaxy
#M33 does not hace a halo rounded to 3 decimals
M33_halo = np.round(Component_mass('M33_000_test.txt',1)/(10**12),3)
M33_disk = np.round(Component_mass('M33_000_test.txt',2)/(10**12),3)

#Calculate the total masses of the three galxies MW,M31,and M33
# rounded to 3 decimals
MW_total = np.round(MW_bulge + MW_disk + MW_halo,3)
M31_total = np.round(M31_bulge + M31_disk + M31_halo,3)
M33_total = np.round(M33_disk + M33_halo,3)

#Calculate the total mass of the three galxies
Group_total = np.round(MW_total + M31_total + M33_total,3)

#Calculate the baryon fraction (f_bar = total stellar mass/total mass)
# for each of the three galxies MW,M31,and M33 rounded to 3 decimals
fbar_MW = np.round((MW_bulge + M31_disk)/MW_total,3)
fbar_M31 = np.round((M31_bulge+ M31_disk)/M31_total,3)
fbar_M33 = np.round(M33_disk/M33_total,3)

#Creating a table (in form of dictionary
# to store all the information that is calcularted above 
# 'column name':[data1,data2,data3]
data_table = {'Galaxy Name':['MW','M31','M33']
              ,'Halo Mass [10^12 SolMasses]':[MW_halo,M31_halo,M33_halo]
              ,'Disk Mass [10^12 SolMasses]':[MW_disk,M31_disk,M33_disk]
              ,'Bulge Mass [10^12 SolMasses]':[MW_bulge,M31_bulge,'NA']
              ,'Total [10^12 SolMasses]':[MW_total,M31_total,M33_total]
              ,'f_bar':[fbar_MW,fbar_M31,fbar_M33]
              ,'Local Group Mass [10^12 SolMasses]':[Group_total,'NA','NA']}



#Calls and prints the table
table = pd.DataFrame(data_table)
print(table)

#Saves the table to an html that can be opened on an internet browser 
table.to_html('HW3_Mass_Table.html')

