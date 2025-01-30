#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 14:06:15 2025

@author: aaronhernandez
"""
'''
This program will take an imputed file and read it. The output will be an 
array of the data stored in the submitted file
'''

import numpy as np #importing the needed modules
import astropy.units as u

def Read(filename):
    '''
    Parameters
    ----------
    filename : Data File

    Returns 
    -------
    None: The data in 'filename' stored in a data array
    '''
    file = open(filename,'r') # open the file and read it
    
    line1 = file.readline() #read first data file line
    label,value = line1.split() #splits the imputs into two variables
    time = float(value)*u.Myr #sets time 
    line2 = file.readline()  #reads the second line
    label2,numb = line2.split() #splits the second line
    particle_num = float(numb) #maks numb its own vaiable
    file.close() #closes the file
    
    data = np.genfromtxt(filename,dtype=None,names = True, skip_header= 3)
    return time,particle_num,data # This will return the time, the particle  
    # number and the data in an array 
    
#Read('MW_000_test.txt')
