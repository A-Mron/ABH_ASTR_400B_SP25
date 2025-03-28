
# # Homework 7 Template
# 
# Rixin Li & G . Besla
# 
# Make edits where instructed - look for "****", which indicates where you need to 
# add code. 



# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
from IPython.display import Latex

# **** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass2 import CenterOfMass

# **** import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass_spydr import Component_mass

# # M33AnalyticOrbit



class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self,fileout): # **** add inputs
        """Class to predict the future trajectory of M33 in the center of 
        mass frame of M31"""

        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        ### **** store the output file name
        self.filename = "%s.txt"%(fileout)
        
        ### get the current pos/vel of M33 
        # **** create an instance of the  CenterOfMass class for M33 
        M33_COM = CenterOfMass("M33_000.txt", 2)
        # **** store the position VECTOR of the M33 COM (.value to get rid of units)
        M33_COM_P = M33_COM.COM_P(0.1,4).value
        # **** store the velocity VECTOR of the M33 COM (.value to get rid of units)
        M33_COM_V = M33_COM.COM_V(M33_COM_P[0], M33_COM_P[1], M33_COM_P[2]).value
        
        ### get the current pos/vel of M31 
        # **** create an instance of the  CenterOfMass class for M31 
        M31_COM = CenterOfMass("M31_000.txt", 2)
        # **** store the position VECTOR of the M31 COM (.value to get rid of units)
        M31_COM_P = M31_COM.COM_P(0.1,2).value
     
        # **** store the velocity VECTOR of the M31 COM (.value to get rid of units)
        M31_COM_V = M31_COM.COM_V(M31_COM_P[0], M31_COM_P[1], M31_COM_P[2]).value
        
        ### store the DIFFERENCE between the vectors posM33 - posM31
        #POS_diff = [M33_COM_P[0]- M31_COM_P[0], M33_COM_P[1]-M31_COM_P[1], M33_COM_P[2]-M31_COM_P[2]]
        # **** create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        self.r0 = M33_COM_P - M31_COM_P
        
        self.v0 = M33_COM_V - M31_COM_V
        
        ### get the mass of each component in M31 
        ### disk
        # **** self.rdisk = scale length (no units)
        self.rdisk = 5.0
        # **** self.Mdisk set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mdisk = Component_mass("M31_000.txt", 2)*1e12
        ### bulge
        # **** self.rbulge = set scale length (no units)
        self.rbulge = 1.0
        # **** self.Mbulge  set with ComponentMass function. Remember to *1e12 to get the right units Use the right ptype
        self.Mbulge = Component_mass("M31_000.txt", 3)*1e12
        # Halo
        # **** self.rhalo = set scale length from HW5 (no units)
        self.rhalo = 62.0
        # **** self.Mhalo set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mhalo = Component_mass("M31_000.txt", 1)*1e12
    
    
    def HernquistAccel(self,M,r_a,r): # it is easiest if you take as an input the position VECTOR 
        """ Returns the Hernquist Acceleration of the Halo and Bulge
        Input:
            M: [float], Total mass of Halo and Bulge
            r_a: [float], Hernquist scale length
            r: [array], position vector
        Output:
            Hern: [array], Hernaquist acceleration vector
        """
        
        ### **** Store the magnitude of the position vector
        rmag = np.sqrt(r[0]**2+r[1]**2+r[2]**2)
        
        ### *** Store the Acceleration
        Hern = -(self.G*M/(rmag *(r_a + rmag)**2))*r
            #follow the formula in the HW instructions
        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        # use  -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        
        return Hern
    
    
    
    def MiyamotoNagaiAccel(self,M_disk,r_d,r):# it is easiest if you take as an input a position VECTOR  r 
        """Returns the Miyamoto-Nagai Acceleration of the Disk
        Input:
            M_disk: [float], Mass of component
            r_d [float], 
            r: [array], position vector
        Output:
            Miya_Nag: [array], Miyamoto-Nagai acceleration vector
        """
        # Set R and B values
        R = np.sqrt(r[0]**2+r[1]**2)
        B = r_d + np.sqrt(r[2]**2 + (r_d/5)**2)
       
        ### Acceleration **** follow the formula in the HW instructions
        # AGAIN note that we want a VECTOR to be returned  (see Hernquist instructions)
        # this can be tricky given that the z component is different than in the x or y directions. 
        # we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        # differences in the z direction:
        #  multiply the whle thing by :   np.array([1,1,ZSTUFF]) 
        # where ZSTUFF are the terms associated with the z direction
        Miya_Nag = -self.G*M_disk*r/((R**2 + B**2)**1.5)*np.array([1,1,
                                    B/np.sqrt(r[2]**2 + (r_d/5)**2)])
        
       
        return Miya_Nag
        # the np.array allows for a different value for the z component of the acceleration
     
    
    def M31Accel(self,r): # input should include the position vector, r
        """Returns the TOTAL acceleration of M31 using the functions
        HernquistAccel and MiyamotoNagaiAccel 
        Input:
            r: [array], position vector
        Output:
            M31_Acc: [array], Total M31 acceleration vector
        """

        ### Call the previous functions for the halo, bulge and disk
        # **** these functions will take as inputs variable we defined in the initialization of the class like 
        # self.rdisk etc. 
        M31_Acc = (self.HernquistAccel(self.Mbulge, self.rbulge, r) 
        + self.HernquistAccel(self.Mhalo, self.rhalo, r) + self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, r)) 
            # return the SUM of the output of the acceleration functions - this will return a VECTOR 
        return M31_Acc
    
    
    
    def LeapFrog(self,dt,r,v): # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        """Returns the new position and velocity vectors based on a 
        LeapFrog integration method
        Input:
            dt: [float], time interval for the integration
            r: [array], initial position vector
            v: [array], initial velocity vector
        Output:
            rnew: [array], new position vector
            vnew: [array], new velocity vector
        """
        
        # predict the position at the next half timestep
        rhalf = r + v*(dt/2)
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        vnew = v + self.M31Accel(rhalf)*dt
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = rhalf + vnew*(dt/2)
        
        return rnew,vnew#return the new position and velcoity vectors
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        """Integrates the equations of motion to compute the future orbit 
        of M33 for 10 Gyr into the future
        Input:
            t0: [float], stating time
            dt: [float], time interval
            tmax: [float], final time
        Output:
            Orbit: [array], 
        """

        # initialize the time to the input starting time
        t = t0
        
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        orbit = np.zeros((int(tmax/dt)+2,7))
        
        # initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        # this above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        r = self.r0
        v = self.v0
        
        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while (t<=tmax):  # as long as t has not exceeded the maximal time 
        
            # **** advance the time by one timestep, dt
            t += dt
            # **** store the new time in the first column of the ith row
            #orbit[i,0] = t 
        
            # ***** advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like:     a,b,c = function(input)
            r,v= self.LeapFrog(dt,r, v)
            
            orbit[i] = t,*tuple(r),*tuple(v)
    
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
            #orbit[i] = t, *tuple(r), *tuple(v)
            
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
           
            # **** update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i += 1
        
        # write the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                       header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                       .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
        # there is no return function

def relative_mag(a): 
    """
    Function that computes the magnitude of a vector a

    PARAMETERS
    ----------
    a : [array], vector 1

    RETURNS
    -------
    mag : [array], Mag(a)
    """
    
    # compute the difference vector
    x = a[0]  
    y = a[1] 
    z = a[2] 

    # return its magnitude
    return np.sqrt(x**2 + y**2 + z**2)

M33 = M33AnalyticOrbit('Orbit_Prediction')
M33_predict = M33.OrbitIntegration(0, 0.1, 9.9)

M31_orbit = np.genfromtxt('Orbit_Prediction.txt') 

Rel_Position = np.array((M31_orbit['x'],M31_orbit['y'],M31_orbit['z']))

Velocity = np.array((M31_orbit['vx'],M31_orbit['vy'],M31_orbit['vz']))


#Plotting
'''
#Possition plot
plt.figure()
plt.xlabel('Time [Gyr]')
plt.ylabel('Total Possition')
plt.title('')
plt.plot(M31_orbit['t'],Rel_Position)
plt.show()

#Velocity Plot
plt.figure()
plt.xlabel('Time [Gyr]')
plt.ylabel('Total Velocity')
plt.title('')
plt.plot(M31_orbit['t'],Velocity)
plt.show()
'''

'''
Due to some problems with my data file I could not make the plots but I have 
seen the plot I am suppost to get so I will answer the questions to the best
of my ability.

2) The plots are very different, the plot form HW6 shows an occilitory pattern 
    for the distence between the galaxies while this proediction is must more 
    quadratic and they have very little simularities. 

3) The physics that come from the MW may be whats missing in the sumilation. 
    the forces from the MW would undoutibly chamge the reletive positions if 
    M31 and M33 so it mush be at least one of the missing links. 

4) To include the effects from the MW you could try and include its effect on 
    the COM position of M31 and that could in turn effect the reletive position
    between M31 and M33 that is used in the calculations. 
'''

