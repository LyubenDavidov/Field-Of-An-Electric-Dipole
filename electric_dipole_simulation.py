# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 20:44:53 2022

@author: ljubo
"""
import numpy as np
import math
import matplotlib.pyplot as plt



import plotly.io as pio
pio.renderers.default='browser'






# Define constants
lambda0 = 700*10**(-9); # Wavelength of light in vacuum in [m] SET YOUR WAVELENGTH
c = 299792458;          # Speed of light in vacuum in [m/s]
f = c/lambda0;          # Frequency of light in [Hz]
w = 2*np.pi*f;          # Frequency of light in [rad/s]

mu_0 = 4*np.pi*10**(-7);# Magnetic permeability of free space [N/A^2]
etha_0 = 1/(mu_0*c**2); # Electric permittivity of free space [F/m]
etha_r = 1;             # Relative electric permittivity
k = w*np.sqrt(etha_0*etha_r*mu_0);



# Miscellaneous functions

def magnitude(vector): 
    return math.sqrt(sum(pow(element, 2) for element in vector));

#test = np.arange(0,3,1).reshape(3,1);
#print(test);
#print(magnitude(test));




# Define vectors 

ax_length = 50;
points = 26; #The interval between -ax_length and +ax_length is divided in points


rx = np.linspace(-ax_length, ax_length, points);
#rx = np.delete(rx, 0);

ry = np.linspace(-ax_length, ax_length, points);
#ry = np.delete(ry, 0);

rz = np.linspace(-ax_length, ax_length, points);
#rz = np.delete(rz, 0);


p = np.array([0,0,1]); #Electric dipole moment (One Debye = 3.33*10^-30 [C/m] => typical for a molecule)



# Define position and field magnitude storage
E_local = np.array([[0,0,0,0]]);


# Define a formula for the electric field magnitude

def EL_FIELD(r):
    r_mag = magnitude(r);
    r_norm = r/r_mag;
    A = (k**2)*np.cross(r_norm, np.cross(p,r_norm))+(3/(r_mag**2))*(np.dot(r_norm, p))*r_norm-(1/r_mag**2)*p;
    B = (k/r_mag)*p-(3*k/r_mag)*np.dot(r_norm, p)*r_norm;
    E_mag = (1/(4*np.pi*etha_0*r_mag))*np.sqrt(A[0]**2+B[0]**2+A[1]**2+B[1]**2+A[2]**2+B[2]**2);
    return E_mag;



# Calculate E field values for all points

for i in rx:
    for j in ry:
        for k in rz:
            r_vec = np.array([i,j,k]);
            E_local = np.append(E_local, [[r_vec[0], r_vec[1], r_vec[2], EL_FIELD(r_vec)]], axis = 0);


X = E_local[:,0];
Y = E_local[:,1];
Z = E_local[:,2];
VALUE = E_local[:,3];



# Plot

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

img = ax.scatter3D(X, Y, Z, c=VALUE, marker='.', cmap="jet")
ax.set_xlabel('$x$', fontsize=20)
ax.set_ylabel('$y$', fontsize=20)
ax.set_zlabel('$z$', fontsize=20)
plt.title("|E(r)| for a dipole moment = ({},{},{})".format(p[0],p[1],p[2], 's'))
fig.colorbar(img)
plt.show()







