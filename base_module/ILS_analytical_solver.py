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
import pandas as pd

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
BHE_num = 3
BHE_wall_points_num = 4
#power
#create 3 dim dataframe to store the st_all for all BHE_wall_points
#the data type in dataframe is:
#axis 0: BHE_wall_points,
#axis 1: BHE num, 
#axis 2: timestep_tot  
st_all_global = np.zeros((BHE_wall_points_num, BHE_num, timestep_tot))


#bhe location
bhe_pos_x = geometry.bhe_pos_x
bhe_pos_y = geometry.bhe_pos_y

#import reference points (borehole wall points)
localVars = locals()
for i in range(BHE_num):
    localVars['bhe_'+ str(i) + '_wall_pos_x' ] = geometry.localVars['bhe_'+ str(i) + '_wall_pos_x' ]
    localVars['bhe_'+ str(i) + '_wall_pos_y' ] = geometry.localVars['bhe_'+ str(i) + '_wall_pos_y' ]


#%% functions
def st_dataframe(step,BHE_id,st):
    #sourceterm dataframe starts from step = 1.
    cur_step = step - 1
    #first step no need
    if cur_step == 0:
        st_all_global[:,:,cur_step] = st
    for i in range(BHE_wall_points_num):
        st_all_global[i,BHE_id,cur_step] = st

def ILS_solver(timestep, bhe_id):
    T_domain=np.zeros([BHE_wall_points_num,timestep])
    coeff_all = np.zeros([BHE_wall_points_num,BHE_num,timestep])
    
    for currstep in range(0,timestep):
        #data container
        dist_bhe_to_ref_po= np.zeros([BHE_wall_points_num,BHE_num])
        localcoeff= np.zeros([BHE_wall_points_num,BHE_num])
        
        for i in range(0,BHE_num):
            #coefficient of current timestep
            for j in range(0,BHE_wall_points_num):
                dist_bhe_to_ref_po[j,i] = (bhe_pos_x[i] - localVars['bhe_'+ str(bhe_id) + '_wall_pos_x' ][j] )**2     \
                                        + (bhe_pos_y[i] - localVars['bhe_'+ str(bhe_id) + '_wall_pos_y' ][j] )**2
                exp1 = dist_bhe_to_ref_po[j,i]/(4*alpha*delta_t*(currstep+1))
                n1 = sp.exp1(exp1)
                localcoeff[j,i] = 1/(4*math.pi*k_s)*n1  
            #coefficient of current timestep after 
            if currstep > 0 :
                for j in range(0,BHE_wall_points_num):
                    dist_bhe_to_ref_po[j,i] = (bhe_pos_x[i] - localVars['bhe_'+ str(bhe_id) + '_wall_pos_x' ][j] )**2     \
                                            + (bhe_pos_y[i] - localVars['bhe_'+ str(bhe_id) + '_wall_pos_y' ][j] )**2
                    exp1 = dist_bhe_to_ref_po[j,i]/(4*alpha*delta_t*currstep)
                    n1 = sp.exp1(exp1)
                    localcoeff[j,i] = localcoeff[j,i] - 1/(4*math.pi*k_s)*n1  
        
        #reverse the coefficient order
        coeff_all[:,:,1:]=coeff_all[:,:,:timestep-1]
        #store each timestep's localcoefficient into global time coefficient dataframe
        coeff_all[:,:,0]= localcoeff 
        
    # get the final temperature field by multiplying
    #the global coefficinet dataframe with the sourceterm matrix 
    for currstep in range(timestep):
        T_domain[:,currstep] =  np.sum(np.sum(coeff_all[:,:,timestep-1-currstep:]
                                *st_all_global[:,:,:currstep+1],axis=1),axis=1) + T0                              
    
    #get the selected BHE' average wall soil temperature by summarizing the
    #all 4 reference points temperature 
    bhe_avg_soil_T = np.sum(T_domain[:,timestep-1],axis=0)/BHE_wall_points_num
    return bhe_avg_soil_T

    