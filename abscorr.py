#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 11:04:10 2024

@author: Jose A. Rodriguez-Rivera
jose.rodriguez@nist.gov
joarodri@gmail.com
University of Maryland
NIST Center for Neutron Research
"""

import numpy as np
import math


def boxsample(lx,ly,lz):
    
    """
    Define a box samle with lx,ly and lz dimentions in mm
    IMPORTANT:
    Boundaries are defined as
    F= Ax^2 + By^2 + Cz^2+ Dxy + Eyz + Gzx +Hx +Iy +Jz = fu(x,y,z)
    The sample must be located at (0,0,0). There is not need to be centered
    The sample must be defined as a set of convex boundaries."   
    """
    
    #We define 6 planes boundaries for the box
    F=np.zeros((6,10),dtype=float)
    #Plane 1
    #plane x=lx/2  ->  H=1  F=lx/2  K=1
    F[0,6]=1.0      #H=1.0
    F[0,9]=lx/2  		#F=lx/2
    #Plane 2
    #plane x=-lx /2 ->  H=1  F=-lx  K=0
    F[1,6]=1.0       #H=1.0
    F[1,9]=-1.0*lx/2      #F=-lx/2
    #Plane 3
    #plane y=ly  ->  I=1  F=ly/2  K=1
    F[2,7]=1.0       #I=1.0
    F[2,9]=ly/2      #F=ly/2
    #Plane 4
    #plane y=-ly/2  ->  I=1  F=-ly/2  K=0
    F[3,7]=1.0       #I=1.0
    F[3,9]=-1.0*ly/2 #F=-ly/2
    #Plane 5
    #plane z=lz/2 ->  J=1, F=lz/2
    F[4,8]=1.0      #J=1.0
    F[4,9]=lz/2     #F=lz/2
    #Plane 6
    #plane z=-lz/2 ->  J=1, F=-lz/2
    F[5,8]=1.0        #J=1.0
    F[5,9]=-1.0*lz/2  #F=-lz/2
    
    VEI=[0,0,0.001]                    #VEI= define volume element inside
    box_bottom=[-lx/2,-ly/2,-lz/2]     #box lower corner
    box_dimensions=[lx,ly,lz]          #ibox volume dimensions
    box_deltas=[20,20,20]              #box number volumeo elements in each direction
    
    sample={'F':F,'VEI':VEI,'box_bottom':box_bottom,'box_dimensions':box_dimensions,'box_deltas':box_deltas}
    return sample

# This procedure returns [t1,t2] distances using the s vector direction
# at the dd(x,y,z) volume element positions for the Fu defined boundary
# function F= Ax^2 + By^2 + Cz^2+ Dxy + Eyz + Gzx +Hx +Iy +Jz.
def tdis(fu,s,dd):
    t1=15.0e20
    t2=15.0e20
    # Defining the boundary function
    A=fu[0]
    B=fu[1]
    C=fu[2]
    D=fu[3]
    E=fu[4]
    G=fu[5]
    H=fu[6]
    I=fu[7]
    J=fu[8]
    F=fu[9]
    # Defining the vector direction S
    sx=s[0]
    sy=s[1]
    sz=s[2]
    # Defining the volume elemet dd
    xp=dd[0]
    yq=dd[1]
    zr=dd[2]
    # Calculate U, V and W
    U = A*sx**2 + B*sy**2 +C*sz**2 + D*sx*sy + E*sy*sz + G*sz*sx
    V = 2*(A*xp*sx + B*yq*sy + C*zr*sz) + D*(xp*sy+yq*sx) + E*(yq*sz +zr*sy) + G*(zr*sx +xp*sz) + H*sx + I*sy + J*sz
    W = A*xp**2 + B*yq**2 +C*zr**2 + D*xp*yq + E*yq*zr + G*zr*xp +H*xp + I*yq + J*zr -F
    if U != 0.0:
        zz= V**2-4*U*W
        if zz >= 0.0:
            t1= (-1*V + math.sqrt(zz))/(2*U)
            t2= (-1*V - math.sqrt(zz))/(2*U)
        else:
            t1=1.5e20
            t2=1.5e20
    else:
        if V == 0:
            t1=1.5E20
        else:
            t1= -1*W/V
        t2= 1.5e20
    treturn=[t1,t2]
    return treturn

# This procedure returns [to1,to2] distances using the s vector direction
# at the dd(x,y,z) volume element positions for the Fu defined boundary
# function F= Ax^2 + By^2 + Cz^2+ Dxy + Eyz + Gzx +Hx +Iy +Jz.
def todis(fu,s,dd):
    t1=15.0e20
    t2=15.0e20
    # Defining the boundary function
    A=fu[0]
    B=fu[1]
    C=fu[2]
    D=fu[3]
    E=fu[4]
    G=fu[5]
    H=fu[6]
    I=fu[7]
    J=fu[8]
    F=fu[9]
    # Defining the vector direction S
    sx=s[0]
    sy=s[1]
    sz=s[2]
    # Defining the volume elemet dd
    xp=dd[0]
    yq=dd[1]
    zr=dd[2]
    # Calculate U, V and W
    U = A*sx**2 + B*sy**2 +C*sz**2 + D*sx*sy + E*sy*sz + G*sz*sx
    V = -1.0*(2*(A*xp*sx + B*yq*sy + C*zr*sz) + D*(xp*sy+yq*sx) + E*(yq*sz +zr*sy) + G*(zr*sx +xp*sz) + H*sx + I*sy + J*sz)
    W = A*xp**2 + B*yq**2 +C*zr**2 + D*xp*yq + E*yq*zr + G*zr*xp +H*xp + I*yq + J*zr -F
    if U != 0.0:
        zz= V**2-4*U*W
        if zz >= 0.0:
            t1= (-1*V + math.sqrt(zz))/(2*U)
            t2= (-1*V - math.sqrt(zz))/(2*U)
        else:
            t1=1.5e20
            t2=1.5e20
    else:
        if V == 0:
            t1=1.5E20
        else:
            t1= -1*W/V
        t2= 1.5e20
    treturn=[t1,t2]
    return treturn

#get t and t0 distances from sample volumeF
# F  Boundary functions
# dd  volume position
# s0 incident beam
# s scattered beam
# returns [tn0,tp0,tn,tp] where:
# tn,tn0 absolute smalest negative distance for the incident and scattered beam
# tp,tp0 absolute smalest positive distance for the incident and scattered beam

def get_t_to_distances(F,dd,s0,s):
    tn=-1.5e20
    tn0=-1.5e20
    tp=1.5E20
    tp0=1.5E20
    for boundary in F:
        [t10,t20]=todis(boundary,s0,dd)
        [t1,t2]=tdis(boundary,s,dd)
        if t1 > 0 and t1 < tp:
            tp=t1
        if t1 < 0 and t1 > tn:
            tn=t1
        if t2 > 0 and t2 < tp:
            tp=t2
        if t2 < 0 and t2 > tn:
            tn=t2
        if t10 > 0 and t10 < tp0:
            tp0=t10
        if t10 < 0 and t10 > tn0:
            tn0=t10
        if t20 > 0 and t20 < tp0:
            tp0=t20
        if t20 < 0 and t20 > tn0:
            tn0=t20
    return [tn0,tp0,tn,tp]
        



def get_t_t0_distances_VEI_out(smpl,ddl,s0l,sl):
    F=smpl["F"]
    outsid0=0
    outsid=0
    tt=-1.5e20
    tt0=-1.5e20
    dd=np.array(ddl)
    s0=np.array(s0l)
    s=np.array(sl)
    [tn0,tp0,tn,tp] = get_t_to_distances(F,dd,-s0,-s)
    if tp0>0 and tp0<1.5E19:
        ddtest=dd+1.01*tp0*s0
        [insid,t0,t] = is_inside(smpl,s0,s0,ddtest)
        if insid==1:
            tt0=t0+t
            outsid0=1     
    if tp>0 and tp<1.5E19:
        ddtest=dd-1.01*tp*s;
        [insid,t0,t] = is_inside(smpl,s,s,ddtest)
        if insid==1:
            tt=t0+t
            outsid=1
    return [tt0,tt,outsid0,outsid]
      
# This function check if the volume element is inside the sample.
# The routine check if there are boundaries between the dd volume position
# and the origin.
# The sample must be defined as a set of convex boundaries.
# If the sample is defined by several bounding functions, the 
# origin [0,0,0] will located at the volume of interest.
# eg: if we define a an sphere x^2+y^2+z^2=5 and a plane z=2.0
# there will be 2 parts of the sphere. And the origin [0,0,0] will be located
# at the sample volume of interest.
# F  Boundary functions
# dd  volume position
# s0 0incident beam vector
# ss scattered beam vector
# returns [inside,t0,t]
# inside=1 if the volume position is inside the samle, 0 otherwise.
#t distance from the sample volume to the beam exit boundary function.
#t0 distance from the sample volume to the bean entrance boundary function.
  
def is_inside(sample,s00,ss,dd1):
    F=sample["F"]
    VEI=np.array(sample["VEI"])
    s0=np.array(s00)
    s=np.array(ss)
    dd=np.array(dd1)
    VEIs0=dd-VEI;
    
    disdd=np.sqrt(VEIs0[0]**2+VEIs0[1]**2+VEIs0[2]**2)
    
    if disdd == 0.0:
        tp1=1.0
    else:    
        s0p=VEIs0/disdd
        sd0=VEI
        insid=False
        t=-1.0
        t0=-1.0
        [tn00,tp00,tn1,tp1]=get_t_to_distances(F,sd0,s0p,s0p) 
    if disdd <= tp1:
        insid=True
        [tn0,tp0,tn,tp]=get_t_to_distances(F,dd,s0,s)
        t=tp
        t0=tp0
    return [insid,t0,t]



#	This function calculates the absorption given an s0,s,F,grids and d.
#	s0: Incoming beam
#	s: exit beam
#	F: bounding functions
#	gridx,gridy,gridz: numberof grid divitions en x, y and z
#	dx,dy,d: Defines thebox enclosing the sample. length (+-d)
#   mu_Ei, muEf: Linear transmission coefficients for Ei and Ef)

def integ(s0,s,samples,mu_Ei,mu_Ef):
    absorp=0
    mm=0.0
    if type(samples)==dict:
        samples=[samples]
    for sample in samples:
        gridx=sample['box_deltas'][0]
        gridy=sample['box_deltas'][1]
        gridz=sample['box_deltas'][2]
        dx=sample['box_dimensions'][0]/2
        dy=sample['box_dimensions'][1]/2
        dz=sample['box_dimensions'][2]/2
        x_bb=sample['box_bottom'][0];
        y_bb=sample['box_bottom'][1];
        z_bb=sample['box_bottom'][2];
        for p in range(gridx):
            xp=x_bb+(dx/gridx)+p*(2*dx/gridx)
            for q in range(gridy):
                yq=y_bb+(dy/gridy)+q*(2*dy/gridy)
                for r in range(gridz):
                    zr=z_bb+(dz/gridz)+r*(2*dz/gridz)
                    element=[xp,yq,zr]
                    [ins,t0,t]=is_inside(sample,s0,s,element)
                    if ins:
                        mm=mm+1
                        t0_out=0
                        t_out=0
                        for smpl_out in samples:
                            if smpl_out != sample:                               
                                [tt0,tt,outsid0,outsid] = get_t_t0_distances_VEI_out(smpl_out,element,s0,s);
                                if outsid0 == 1:
                                    t0_out = t0_out + tt0;
                                if outsid == 1:
                                    t_out = t_out + tt;
                        absorp=absorp + math.exp(-1.0*(mu_Ei*t0+mu_Ef*t))
    return absorp/mm

  

# rotate the s vector around x-y plane alpha degrees
def srotxy(s,alpha):
    alpharad=math.radians(alpha)
    [sx,sy,sz]=s
    ssx=sx*math.cos(alpharad)-sy*math.sin(alpharad)
    ssy=sx*math.sin(alpharad)+sy*math.cos(alpharad)
    return [ssx, ssy, sz]


# rotate the s vector around y-z plane alpha degrees
def srotyz(s,alpha):
    alpharad=math.radians(alpha)
    [sx,sy,sz]=s
    ssy=sy*math.cos(alpharad)-sz*math.sin(alpharad)
    ssz=sy*math.sin(alpharad)+sz*math.cos(alpharad)
    return [sx, ssy, ssz]


def srotzx(s,alpha):
    """rotate the s vector around z-x plane alpha degrees"""
    alpharad=math.radians(alpha)
    [sx,sy,sz]=s
    ssx=sx*math.cos(alpharad)+sz*math.sin(alpharad)
    ssz=-1*sx*math.sin(alpharad)+sz*math.cos(alpharad)
    return [ssx, sy, ssz]


 
def generate_Sample_theta_2theta(range_theta,range_twotheta,**kwargs):
    """
    Parameters
    ----------
    range_theta : float
        Sample theta range list  [Min, Max, Steps]
    range_twotheta : float
        Sample twotheta range  ist  [Min, Max, Steps]
    **kwargs : File name
        String file name to save the file with generated theta-twotheta angles.

    Returns
    -------
    Return a list of theta-twotheta angles.
    if file='file' the theta-twotheta values will be saved in a file
    """   
  
    if "file" in kwargs:
        file2write=kwargs["file"]
        ff=open(file2write, 'w')
    fil=[]
    for theta in range(range_theta[0],range_theta[1],range_theta[2]):
        for twotheta in range(range_twotheta[0],range_twotheta[1],range_twotheta[2]):
            fil.append([theta,twotheta])
            if "file" in kwargs:
                np.savetxt(ff,[[theta,twotheta]])
    if "file" in kwargs:
        ff.close()        
    return(fil) 





        
    
    
    