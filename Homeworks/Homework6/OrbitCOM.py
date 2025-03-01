

# Homework 6 Template
# G. Besla & R. Li




# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile_spydr import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass2 import CenterOfMass




def OrbitCOM(galaxy,start_snap=0,end_snap=800,n=5):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
          
    outputs: 
    """
    # compose the filename for output
    if galaxy == 'MW':
        fileout = "Orbit.MW.txt"
    elif galaxy == 'M31':
        fileout = "Orbit.M31.txt"
    else:
        fileout = "Orbit.M33.txt"
    # remove all but the last 3 digits
   
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    delta = 0.1
    voldec = 2
    M33voldec = 4
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    snap_ids = np.arange(start_snap,end_snap,n,dtype=int)
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros([len(snap_ids),7])
    
    # a for loop 
    for i,snap_id in enumerate(snap_ids): # loop over files
        
        # compose the data filename (be careful about the folder)
        dname = '../'+galaxy+'/'
        ilble = '000' + str(snap_id)
        ilble = ilble[-3:]
        filename = dname +"%s_"%(galaxy) + ilble + '.txt'
        # Initialize an instance of CenterOfMass class, using disk particles
        COM = CenterOfMass(filename,2)
        
        # Store the COM pos and vel. Remember that now COM_P required VolDec
        if galaxy == 'M33':
            comP = COM.COM_P(delta, M33voldec)
            comV = COM.COM_V(comP[0],comP[1],comP[2])
        else:
            comP = COM.COM_P(delta, voldec)
            comV = COM.COM_V(comP[0],comP[1],comP[2])
    
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        time = COM.time.value/1000
    
        orbit[i] = np.array([time,*tuple(comP.value),*tuple(comV.value)])
        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))




# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 
#Milky way
MW_orbit = OrbitCOM('MW',0,800) #n is set equal to 5 in the parameters

#M31
M31_orbit = OrbitCOM('M31',0,800)

#M33
M33_orbit = OrbitCOM('M33',0,800)




# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt
MWCOM = np.genfromtxt('Orbit_MW.txt',dtype=None,names=True)

M31COM = np.genfromtxt('Orbit_M31.txt',dtype=None,names=True)

M33COM = np.genfromtxt('Orbit_M33.txt',dtype=None,names=True)


# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  



'''
Code not finished, the stuff below shows my thought process but 
I can wrap my head around this. Ill be in office hours next week
with questions
'''
# Determine the magnitude of the relative position and velocities 

# of MW and M31
diffposmw31 = np.sqrt((MW_orbit[1]-M31_orbit[1])**2+(MW_orbit[2]-M31_orbit[2])**2
                      +(MW_orbit[3]-M31_orbit[3])**2)

diffvelmw31 = np.sqrt((MW_orbit[4]-M31_orbit[4])**2+(MW_orbit[5]-M31_orbit[5])**2
                      +(MW_orbit[6]-M31_orbit[6])**2)
# of M33 and M31
diffpos3331 = np.sqrt((M31_orbit[1]-M33_orbit[1])**2+(M31_orbit[2]-M33_orbit[2])**2
                      +(M31_orbit[3]-M33_orbit[3])**2)

diffvel3331 = np.sqrt((M31_orbit[4]-M33_orbit[4])**2+(M31_orbit[5]-M33_orbit[5])**2
                      +(M31_orbit[6]-M33_orbit[6])**2)


# Plot the Orbit of the galaxies 
#################################
fig,ax = plt.subplots(figsize=(5,5))
plt.xlabel('Time [Gys]')
plt.ylabel('Seperation [Kpc]')

plt.plot(MW_orbit[0],diffposmw31,lebel = 'MW-M31')

plt.plot(MW_orbit[0],diffpos3331,lebel = 'M31-M33')

plt.show()

# Plot the orbital velocities of the galaxies 
#################################
fig2,ax2 = plt.subplots(figsize=(5,5))
plt.xlabel('Time [Gys]')
plt.ylabel('Velocity [Kpc/Gys]')

plt.plot(MW_orbit[0],diffvelmw31,lebel = 'MW-M31')

plt.plot(MW_orbit[0],diffvel3331,lebel = 'MW-M31')

plt.show()