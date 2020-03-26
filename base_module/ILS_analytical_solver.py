# -*- coding: utf-8 -*- {}
'''
IFS analytical solver
# function used to calculate the analytical solution of soil temperature
# distribution in an Borehole Heat Exchangers (BHEs) array proposed by Bayer2014 
# Peter Bayer et. al(2014), Strategic optimization of borehole heat exchanger field 
# for seasonal geothermal heating and cooling, Applied Engineering 136: 445-453. 
# http://dx.doi.org/10.1016/j.apenergy.2014.09.029
input: computed source term on each BHE point from Beier analytical solution on each timestep_tot
output: the soil temperature at the distance of borehole wall boundary. 

Author: Shuang Chen
'''
import numpy as np
from scipy import special as sp
import math

import base_module.geometry as geometry

#%% basic input
#soil
T0 = 273.15 + 15 #K soil initial temperature
rho_s = 1120 #kg/m3
c_s = 1780 #J/kg*K
k_s = 2.4 #W/m*K
alpha = k_s/(rho_s*c_s) #m^2/s

#time
time_tot = 10*24*60*60 #s
delta_t = 86400 #s
timestep_tot = int(time_tot/delta_t)

#BHE 
BHE_num = geometry.bhe_num
BHE_wall_points_num = geometry.BHE_wall_points_num
BHE_wall_points_num_all = np.size(geometry.bhe_wall_pos_x)
#power
#create 3 dim dataframe to store the st_all for all BHE_wall_points
#the data type in dataframe is:
#axis 0: BHE_wall_points,
#axis 1: BHE num, 
#axis 2: timestep_tot  
st_all_global = np.zeros((BHE_wall_points_num_all, BHE_num, timestep_tot))


#bhe location
bhe_pos_x = geometry.bhe_pos_x
bhe_pos_y = geometry.bhe_pos_y

#import reference points (borehole wall points)
bhe_wall_pos_x = geometry.bhe_wall_pos_x
bhe_wall_pos_y = geometry.bhe_wall_pos_y

#%% functions
#global sourceterm array data container
def st_dataframe(step,BHE_id,st):
    #sourceterm dataframe starts from step = 1.
    cur_step = step - 1
    #first step no need
    if cur_step == 0:
        st_all_global[:,:,cur_step] = st
    for i in range(BHE_wall_points_num_all):
        st_all_global[i,BHE_id,cur_step] = st

#global coefficient data container(sourceterms to each reference point) of whole timesteps
def ILS_solver_global_coeff():
    coeff_all = np.zeros([BHE_wall_points_num_all,BHE_num,timestep_tot])
    
    for currstep in range(0,timestep_tot):
        #data container
        dist_bhe_to_ref_po= np.zeros([BHE_wall_points_num_all,BHE_num])
        localcoeff= np.zeros([BHE_wall_points_num_all,BHE_num])
        
        for i in range(0,BHE_num):
            #coefficient of current timestep
            for j in range(0,BHE_wall_points_num_all):
                dist_bhe_to_ref_po[j,i] = (bhe_pos_x[i] - bhe_wall_pos_x[j] )**2     \
                                        + (bhe_pos_y[i] - bhe_wall_pos_y[j] )**2
                exp1 = dist_bhe_to_ref_po[j,i]/(4*alpha*delta_t*(currstep+1))
                n1 = sp.exp1(exp1)
                localcoeff[j,i] = 1/(4*math.pi*k_s)*n1  
            #coefficient of current timestep after 
            if currstep > 0 :
                for j in range(0,BHE_wall_points_num_all):
                    dist_bhe_to_ref_po[j,i] = (bhe_pos_x[i] - bhe_wall_pos_x[j] )**2     \
                                            + (bhe_pos_y[i] - bhe_wall_pos_y[j] )**2
                    exp1 = dist_bhe_to_ref_po[j,i]/(4*alpha*delta_t*currstep)
                    n1 = sp.exp1(exp1)
                    localcoeff[j,i] = localcoeff[j,i] - 1/(4*math.pi*k_s)*n1  
        
        #reverse the coefficient order
        coeff_all[:,:,1:]=coeff_all[:,:,:timestep_tot-1]
        #store each timestep's localcoefficient into global time coefficient dataframe
        coeff_all[:,:,0]= localcoeff 
        
    return coeff_all
    
    
def ILS_solver(timestep):
    #soil wall temperature data container on the current timestep      
    T_domain=np.zeros(BHE_wall_points_num_all)
    # get the current timestep temperature field by multiplying
    #the global coefficinet dataframe with the sourceterm matrix 
    T_domain[:] =  np.sum(np.sum(coeff_all[:,:,(timestep_tot - timestep):]
                                *st_all_global[:,:,:timestep],axis=1),axis=1) + T0                              
    
    #get each BHE' average wall soil temperature by summarizing the
    #all 4 reference points temperature of each BHE. Then combine all BHEs
    #average wall soil temperatures into one array and send back.
    bhes_avg_wallsoil_T_array = np.zeros(BHE_num)
    #summarizing each 4 points temperature as one BHE's average wall soil temperature
    for i in range(BHE_num):
        bhes_avg_wallsoil_T_array[i] = np.sum(T_domain[i * BHE_wall_points_num:
                                                      (i + 1) * BHE_wall_points_num
                                                ])/BHE_wall_points_num
    return bhes_avg_wallsoil_T_array

#%%initialise global coefficient data container. In the whole main calculation
#procedure the global coefficinet array do not need to change.
coeff_all = ILS_solver_global_coeff()