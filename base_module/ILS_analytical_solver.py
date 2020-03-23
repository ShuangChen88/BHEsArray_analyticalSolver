# -*- coding: utf-8 -*- {}
'''
IFS analytical solver
# function used to calculate the analytical solution of soil temperature
# distribution in an Borehole Heat Exchangers (BHEs) array proposed by Bayer2014 
# Peter Bayer et. al(2014), Strategic optimization of borehole heat exchanger field 
# for seasonal geothermal heating and cooling, Applied Engineering 136: 445-453. 
# http://dx.doi.org/10.1016/j.apenergy.2014.09.029
input: computed source term on each BHE point from Beier analytical solution on each timestep
output: the soil temperature at the distance of borehole wall boundary. 

Author: Shuang Chen
'''
import numpy as np
from scipy import special as sp
import math
import pandas as pd

import base_module.geometry

#%% basic input
#soil
T0 = 15 #soil initial temperature
rho_s = 1120 #kg/m3
c_s = 1851.4 #J/kg*K
k_s = 1.8 #W/m*K
alpha = k_s/(rho_s*c_s) #m^2/s

#time
timestep = 12
time_trans = 3600 #s
#BHE 
BHE_num = 9
BHE_wall_points_num = 4
#power
'''
TODO: here the interface to Bayer analytical solution

'''
#power
st1 = np.array([10,10,10,10,10,10,-20,-20,-20,-20,-20,-20]).reshape(1,-1)
st2 = np.array([30,30,30,30,30,30,-50,-50,-50,-50,-50,-50]).reshape(1,-1)
st3 = np.array([10,10,10,10,10,10,-20,-20,-20,-20,-20,-20]).reshape(1,-1)
st4 = np.array([10,10,10,10,10,10,-20,-20,-20,-20,-20,-20]).reshape(1,-1)
st5 = np.array([30,30,30,30,30,30,-50,-50,-50,-50,-50,-50]).reshape(1,-1)
st6 = np.array([10,10,10,10,10,10,-20,-20,-20,-20,-20,-20]).reshape(1,-1)
st7 = np.array([10,10,10,10,10,10,-20,-20,-20,-20,-20,-20]).reshape(1,-1)
st8 = np.array([30,30,30,30,30,30,-50,-50,-50,-50,-50,-50]).reshape(1,-1)
st9 = np.array([10,10,10,10,10,10,-20,-20,-20,-20,-20,-20]).reshape(1,-1)

#add all the individul sourceterm into a globle sourceterm array.
st_all = np.concatenate((st1,st2,st3,st4,st5,st6,st7,st8,st9), axis =0)
#create 3 dim dataframe to store the st_all for all BHE_wall_points
#the data type in dataframe is:
#axis 0: BHE_wall_points,
#axis 1: BHE num, 
#axis 2: timestep  
st_all_global = np.zeros((BHE_wall_points_num, BHE_num, timestep))
for i in range(BHE_wall_points_num):
    st_all_global[i,:,:] = st_all


#bhe location
bhe_pos_x = geometry.bhe_pos_x
bhe_pos_y = geometry.bhe_pos_y

#import reference points (borehole wall points)
localVars = locals()
for i in range(BHE_num):
    localVars['bhe_'+ str(i) + '_wall_pos_x' ] = geometry.localVars['bhe_'+ str(i) + '_wall_pos_x' ]
    localVars['bhe_'+ str(i) + '_wall_pos_y' ] = geometry.localVars['bhe_'+ str(i) + '_wall_pos_y' ]


#temp:
bhe_id = 1    
#%% function
def ILS_solver(timestep, bhe_id):
    T_domain=np.zeros([BHE_wall_points_num,timestep + 1])
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
                exp1 = dist_bhe_to_ref_po[j,i]/(4*alpha*time_trans*(currstep+1))
                n1 = sp.exp1(exp1)
                localcoeff[j,i] = 1/(4*math.pi*k_s)*n1  
            #coefficient of current timestep after 
            if currstep > 0 :
                for j in range(0,BHE_wall_points_num):
                    dist_bhe_to_ref_po[j,i] = (bhe_pos_x[i] - localVars['bhe_'+ str(bhe_id) + '_wall_pos_x' ][j] )**2     \
                                            + (bhe_pos_y[i] - localVars['bhe_'+ str(bhe_id) + '_wall_pos_y' ][j] )**2
                    exp1 = dist_bhe_to_ref_po[j,i]/(4*alpha*time_trans*currstep)
                    n1 = sp.exp1(exp1)
                    localcoeff[j,i] = localcoeff[j,i] - 1/(4*math.pi*k_s)*n1  
        
        #reverse the coefficient order
        coeff_all[:,:,1:]=coeff_all[:,:,:timestep-1]
        #store each timestep's localcoefficient into global time coefficient dataframe
        coeff_all[:,:,0]= localcoeff 
        
    # get the final temperature field by multiplying
    #the global coefficinet dataframe with the sourceterm matrix 
    for currstep in range(timestep):
        T_domain[:,currstep+1] =  np.sum(np.sum(coeff_all[:,:,timestep-1-currstep:]*st_all_global[:,:,:currstep+1],axis=1),axis=1)                              
    
    #add T0 into T_domain
    T_domain = T_domain +T0   
    
    return T_domain

    