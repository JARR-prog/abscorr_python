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
# we calculate the absorption of the box given the theta-twotheta position.
# The incident beam s0 is along the y axis s0=[0,1,0] s=[0,1,0] for theta=2theta=0


import numpy as np
import abscorr as abcr
import multiprocessing as multip
import matplotlib.pyplot as plt
import time



def boxsmpltransm(theta,twotheta,**kwargs):
    """

    Parameters
    ----------
    theta : float
        theta value
    twotheta : fluat
        twotheta value.

    Returns
    -------
    absorption : fluat
        Calculated transmission given a theta and twotheta value for a sample define as a box.

    """

    s0=[0.0,1.0,0.0]    #   Incident beam
    s=[0.0,1.0,0.0]     #   Scattered beam  
    mu_ei=2.0           #   Linear absorption factor for the incident beam
    mu_ef=2.0           #   Linear absorption factor for the scattered beam
    samplebox=abcr.boxsample(1,1.55,3.1)
    ss0=abcr.srotxy(s0,-1*theta)
    ss=abcr.srotxy(abcr.srotxy(s,-1*theta),twotheta)
    transm=abcr.integ(ss0, ss, samplebox,mu_ei,mu_ef)
    if 'file_out' in kwargs:
       file2write=kwargs["file_out"]
       ff=open(file2write, 'a')
       np.savetxt(ff,[[theta,twotheta,transm]])
       ff.close()
       
    if 'return_list' in kwargs:
        return_list=kwargs["return_list"] 
        return_list.append([theta,twotheta,transm])       
    return transm



def trans_sample_sequential(angles,**kwargs):
    """
    Parameters
    ----------

    angles:  
                angles :    Two columns list [theta, twotheta]
                            List with two columns list of floats
                            [theta, twotheta] values to compute transmission.
                            
                            String
                            if  angles is a "string" then angles is a file name 
                            with the angles theta_twotheta.
                            The file  must have two columns list of floats  [theta, twotheta]
                                                                             
                         
    **kwargs :  file_out :  file to save the calculated transmission with 3 columns:
                            [theta, twotheta, transmission]

    Returns
    -------
                A list with with three columns [theta, twotheta, transmission]
    """ 
    
    if type(angles)==str:
        file2read=angles
        ff=open(file2read, 'r')
        angles=np.loadtxt(ff)
        ff.close()
        
    if 'file_out' in kwargs:
       file2write=kwargs["file_out"]
       ff=open(file2write, 'w')
       ff.close()
    
    transmi=[boxsmpltransm(angle[0],angle[1],**kwargs) for angle in angles]
        
    
    return transmi
    
    



      
              


              
# calculate the absorption
# in a multitasking loop
def trans_sample_multitasking(angles,**kwargs):
    
    startime=time.time()
    
    if type(angles)==str:
        file2read=angles
        ff=open(file2read, 'r')
        angles=np.loadtxt(ff)
        ff.close()
        
    if 'file_out' in kwargs:
       file2write=kwargs["file_out"]
       ff=open(file2write, 'w')
       ff.close()
       
   
    manager = multip.Manager()
    return_list = manager.list()

    jobs = []

    for tt in angles:
        kwargs['return_list']=return_list
        p = multip.Process(target=boxsmpltransm, args=(tt[0], tt[1]),kwargs=kwargs)
        jobs.append(p)
        p.start()
    for proc in jobs:
        proc.join()
    print('done.', time.time()-startime)
    return(return_list)


