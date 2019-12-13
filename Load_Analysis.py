# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 14:25:20 2019

@author: waldx091
"""

#--------------------- IMPORT MODULES ------------------------------
# Imort FLORIS
import floris
import numpy as np
from math import pi
import floris.tools as tool
import floris.simulation
from floris.utilities import Vec3

# Import CCBlade
from os import path
from ccblade import  CCAirfoil, CCBlade


    
def RootMoments(y_cut_value, blade_radius, yaw_angles, turb_idx, r, chord,  theta, mu, pitch, BN, af, json_file):
    
    """Calculate Blade Root Moments as a funtion of Azimuth angle (every two degrees).
    Parameters
    ----------
    y_cut_value : float (meters)
        value at where wind speed profile slice should be taken
    blade_radius : float (meters)
        length of turbine blade
    yaw_angles : array_like (deg)
        yaw angles for each turbine respectively
    turb_idx : integer
        the index of the turbine being analyzed
    r : array_like (m)
        locations defining the blade along z-axis of :ref:`blade coordinate system <azimuth_blade_coord>`
        (values should be increasing).
    chord : array_like (m)
        corresponding chord length at each section
    theta : array_like (deg)
        corresponding :ref:`twist angle <blade_airfoil_coord>` at each section---
        positive twist decreases angle of attack.
    mu : float
        dynamic viscosity
    pitch : fload (deg)
        blade pitch angle
    BN : integer
        Number of blades analyzed on turbine
    af : database
        Holds cl and cd airfoil informatin for CCBlade
    json_file : ""
        name of json file used for floris simulation. Example: "example_input1.json"
    Returns
    -------
    Np_int : array_like
        The total normal force on blade at each azimuth angle
    Np : array_like
        Normal force on blade as a function of azimuth(in direction of wind)
    Root_Moments : array_like
        Root Moments on blade as a function of azimuth
    notes
    -------
    this function calculates the blade root moments as a function of azimuth angle while
    taking into account the turbulence intensity on the rotor plane"""
        
    ### --------------------- Define Inputs ------------------------------ ###
    # Initialize FLORIS
    json_file = json_file
    F = tool.floris_utilities.FlorisInterface(json_file)
    
    r = r
    chord = chord
    theta = theta
    mu = mu
    tidx = turb_idx
    pitch = pitch
    BN = BN
    af = af
    
    blade_radius = np.ceiling(blade_radius/5) #convert length in m to length in grid boxes
    
    Rhub = 0.0
    Rtip = F.floris.farm.turbines[tidx].rotor_radius
    B = 1
    rho = F.floris.farm.turbines[tidx].air_density
    tilt = F.floris.farm.turbines[tidx].tilt_angle
    yaw = F.floris.farm.turbines[tidx].yaw_angle
    hubHt = F.floris.farm.turbines[tidx].hub_height
    tsr = F.floris.farm.turbines[tidx].tsr #?              
    
    yaw_angles = yaw_angles
    F.calculate_wake(yaw_angles=yaw_angles)
    
    # Vertical cut of flow field
    y_cut_value = y_cut_value                                                                  # where the u_in velocities will be measured
    grid_spacing = 5 #(meters)
    v_plane = tool.cut_plane.VertPlane(F.get_flow_data(grid_spacing=grid_spacing),y_cut_value)
    
    ### --------------------- EXTRACT WIND SPEEDS ------------------------------ ###
    
    # Extracting incoming velocities from vertical cut
    u = v_plane.u_in
    v = v_plane.v_in   # spanwise wind speed will be zero unless using 'curl' model and if a turbine has yaw
    
    # converting to comma separated list
    u = list(u)
    v = list(v)
    
    # define the resulution of the cut plane grid
    res = np.flip(v_plane.resolution)
    
    # combine each horizontal split into a 2d matrix
    u_2d = np.array(u).reshape(res)
    vy_2d = np.array(v).reshape(res)
    
    # define the location of the tubine hub in the 2d array
    center_y = int(round(F.floris.farm.turbines[tidx].hub_height/5))-1
    center_x = 50    # need to figure out way to map layout_x from floris to center_x index value
    
    # take only the velocity in contact with the turbine
    for i in range(len(u_2d)):
        for j in range(len(u_2d[i])):
            if np.sqrt((i-(center_y))**2+(j-(center_x))**2) > blade_radius + 2:
                u_2d[i][j] = 0
    
    for k in range(len(vy_2d)):
        for l in range(len(vy_2d[k])):
            if np.sqrt((k-(center_y))**2+(l-(center_x))**2) > blade_radius + 2:
                vy_2d[k][l] = 0
                
    # calculating grid plane Turbulence intensity
    
#    ---------------turbulence intensity from floris-------------
#    velocity_model = floris.simulation.WakeVelocity([0])
#    turb_coord = [0,0,90]
#    turbine_coord = floris.utilities.Vec3(turb_coord)
#    w_coord = [0,400,90]
#    wake_coord = floris.utilities.Vec3(w_coord)
#    turbine_wake = F.floris.farm.turbines[1]
#    TI = F.floris.farm.turbines[0].calculate_turbulence_intensity(flow_field_ti, velocity_model, turbine_coord, wake_coord, turbine_wake)
    
    TI = 0.15 #will be calulated using floris funciton in future(0.15 is a random value)
    
    ### --------------------- Create CCBlade rotor ------------------------------ ###
    
    rotor = CCBlade(r, chord, theta, af, Rhub, Rtip, B, rho, mu, tilt, yaw, hubHt)
    
    ### --------  OBTAIN WIND SPEEDS AT EACH BLADE SECTOIN AT EACH AZIMUTH ---------------- ###
    BR = len(blade_radius)
    
    azimuth = np.zeros(shape=(180,1))
    x = np.zeros(shape=(180,13))
    y = np.zeros(shape=(180,13))
    x_coord = np.zeros(shape=(180,BR))
    y_coord = np.zeros(shape=(180,BR))
    blade_sections = np.zeros(shape=(BR,1))
    Uinf = np.zeros(shape=(180,BR))
    vy = np.zeros(shape=(180,BR))
    x0 = np.zeros(shape=(180,BR))
    y0 = np.zeros(shape=(180,BR))
    x1 = np.zeros(shape=(180,BR))
    y1 = np.zeros(shape=(180,BR))
    Q11 = np.zeros(shape=(180,BR))
    Q21 = np.zeros(shape=(180,BR))
    Q12 = np.zeros(shape=(180,BR))
    Q22 = np.zeros(shape=(180,BR))
    T11 = np.zeros(shape=(180,BR))
    T21 = np.zeros(shape=(180,BR))
    T12 = np.zeros(shape=(180,BR))
    T22 = np.zeros(shape=(180,BR))
    P = np.zeros(shape=(180,BR))
#    S = np.zeros(shape=(180,BR))
#    TI_values = np.zeros(shape=(180,BR))
    
    for ii in range(len(np.linspace(0,360,180))):
        azimuth[ii] = ii*2                           # analyzing the loads at an azimuth angle every 2 degrees
        for jj in range(len(np.linspace(1,BR,BR))):
            blade_sections[jj] = jj+1                # checking wind speed at every 5m section of blade
            x[ii,jj] = blade_sections[jj] * (np.cos(np.deg2rad(azimuth[ii])))
            y[ii,jj] = blade_sections[jj] * (np.sin(np.deg2rad(azimuth[ii])))
            x_coord[ii,jj] = center_x + x[ii,jj]
            y_coord[ii,jj] = center_y + y[ii,jj]
            vy[ii,:] = vy_2d[y_coord[ii,:].astype(int),x_coord[ii,:].astype(int)]
        
    # ---------------------bilinear interpolation at each point------------------------
    
            x0[ii,jj] = np.floor(x_coord[ii,jj])
            x0 = x0.astype(int)
            y0[ii,jj] = np.floor(y_coord[ii,jj])
            y0 = y0.astype(int)
            x1[ii,jj] = np.floor(x_coord[ii,jj])+1
            x1 = x1.astype(int)
            y1[ii,jj] = np.floor(y_coord[ii,jj])+1
            y1 = y1.astype(int)
            Q11[ii,jj] = u_2d[y1[ii,jj],x1[ii,jj]]
            Q21[ii,jj] = u_2d[y1[ii,jj],x0[ii,jj]]
            Q12[ii,jj] = u_2d[y0[ii,jj],x1[ii,jj]]
            Q22[ii,jj] = u_2d[y0[ii,jj],x0[ii,jj]]
            
            T11[ii,jj] = TI[y1[ii,jj],x1[ii,jj]]
            T21[ii,jj] = TI[y1[ii,jj],x0[ii,jj]]
            T12[ii,jj] = TI[y0[ii,jj],x1[ii,jj]]
            T22[ii,jj] = TI[y0[ii,jj],x0[ii,jj]]
            P[ii,jj] = (((x0[ii,jj]-x_coord[ii,jj])*(y0[ii,jj]-y_coord[ii,jj]))/((x0[ii,jj]-x1[ii,jj])*(y0[ii,jj]-y1[ii,jj])))*Q11[ii,jj] + (((x_coord[ii,jj]-x1[ii,jj])*(y0[ii,jj]-y_coord[ii,jj]))/((x0[ii,jj]-x1[ii,jj])*(y0[ii,jj]-y1[ii,jj])))*Q21[ii,jj] + (((x0[ii,jj]-x_coord[ii,jj])*(y_coord[ii,jj]-y1[ii,jj]))/((x0[ii,jj]-x1[ii,jj])*(y0[ii,jj]-y1[ii,jj])))*Q12[ii,jj] + (((x_coord[ii,jj]-x1[ii,jj])*(y_coord[ii,jj]-y1[ii,jj]))/((x0[ii,jj]-x1[ii,jj])*(y0[ii,jj]-y1[ii,jj])))*Q22[ii,jj] 
            Uinf[ii,:] = P[ii,:]
#            S[ii,jj] = (((x0[ii,jj]-x_coord[ii,jj])*(y0[ii,jj]-y_coord[ii,jj]))/((x0[ii,jj]-x1[ii,jj])*(y0[ii,jj]-y1[ii,jj])))*T11[ii,jj] + (((x_coord[ii,jj]-x1[ii,jj])*(y0[ii,jj]-y_coord[ii,jj]))/((x0[ii,jj]-x1[ii,jj])*(y0[ii,jj]-y1[ii,jj])))*T21[ii,jj] + (((x0[ii,jj]-x_coord[ii,jj])*(y_coord[ii,jj]-y1[ii,jj]))/((x0[ii,jj]-x1[ii,jj])*(y0[ii,jj]-y1[ii,jj])))*T12[ii,jj] + (((x_coord[ii,jj]-x1[ii,jj])*(y_coord[ii,jj]-y1[ii,jj]))/((x0[ii,jj]-x1[ii,jj])*(y0[ii,jj]-y1[ii,jj])))*T22[ii,jj] 
#            TI_values[ii,:] = S[ii,:]
            
    Omega = (np.mean(Uinf))*tsr/Rtip * 30.0/pi
        
    ### --------------------- Introducing Turbulence intensity ------------------------------ ###
    
    # calculating standard deviation of Turbulence Intensity values
#    a_infs = np.where(np.isnan(TI_values))
#    TI_values[a_infs] = 0.0
#    TI_std = np.std(TI_values)
#    print(TI_std)
    
    # new method for implementing TI:
    Uinf_new = np.zeros(shape=(180,BR))
    TI_rand = np.zeros(shape=(180,BR))
    Tp = np.zeros(shape=(180,BR))
    Np = np.zeros(shape=(180,BR))
    Root_Moment = np.zeros(shape=(180,1))
    Np_int = np.zeros(shape=(180,1))
    rn = np.zeros(shape=(BR,1))
    for ii in range(len(np.linspace(0,360,180))):
        for jj in range(len(np.linspace(1,BR,BR))):
            TI_rand[ii,jj] = (np.random.random_sample(1)-.5)*TI
            Uinf_new[ii,jj] = Uinf[ii,jj] + (Uinf[ii,jj]*TI_rand[ii,jj])
            Np[ii],Tp[ii] = rotor.distributedAeroLoads(Uinf_new[ii,:], Omega, pitch, 0, vy[ii,:])
            rn[jj] = r[jj] -r[jj-1]
            rn[0] = r[0]
            Np_int[ii,:] = np.sum(np.multiply(Np[ii,:], rn[:]))
                
    Root_Moment = Np_int * (Rtip / 2)
    
    # Multiplying by number of blades on turbine
    Np_int = Np_int * BN
    Root_Moment = Root_Moment * BN
    return Np_int, Root_Moment, Np



def FatigueAnalysis(Root_Moment, sigma_0, m, cycles, BN):
    """Calculate damage equivalent loads from fatigue loading on a blade.
    Parameters
    ----------
    Root_Moment : array_like (N*m)
        Blade rood moment calculated for each azimuth for one rotation
    RC : integer
        The number of random turbulence intensity cases
    sigma_0 : float
        Initial S-N curve stress value
    m : float
        S-N curve slope
    cycles : float
        specifies number of cycles the turbine will experience for each situation
    BN : integer
        specifies the number of blades on the turbine
    Returns
    -------
    D_tot : float
        Damage equivalent load for specific situation
    notes
    -------
    this function calculates the damage equivalent loads for a specified
    number of cycles using Miners rule. Random cases due to turbulence intensity
    are taken into account. The function assumes that each random case occurs
    for a fraction of the total cycles.
    """
    
    Root_Moment = Root_Moment / BN
    
#    cycle_frac = (1 / RC) * cycles
    
#    RM_DIFF = np.zeros(shape=(RC,1))
#    for q in range(len(np.linspace(1,RC, RC))):
#        RM_DIFF[q] = np.max(Root_Moment[q]) - np.min(Root_Moment[q])
    
    RM_DIFF = np.max(Root_Moment) - np.min(Root_Moment)
    
    sigma_a = RM_DIFF
    
    S_N_cycles = (sigma_a - sigma_0) / m
    
    D = (1 / S_N_cycles)
    
    D_tot = D * cycles
    
    return D_tot
        
        