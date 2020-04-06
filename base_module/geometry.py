# -*- coding: utf-8 -*- {}
'''
Geometry of the analytical model
bhe order: upperleft to upper right, upper to down
'''
import numpy as np
from pandas import read_csv

#%% User setting
#input geometry of the model

bhe_num = 6 # virtual BHE number
bhe_num_real = 25 # total BHE number in model
BHE_wall_points_num = 4 #4 reference points on each BHE wall
br = 0.063 # borehole raduis in m

#End User setting
#%% real BHEs coordinates
#note: the sequence of the BHEs should be same as the BHEs sequence defined in the global_bhes_sequence_table.csv
real_df_nw = read_csv('./base_module/pre/global_bhes_sequence_table.csv',
                    delimiter=';')

real_bhe_pos_x = real_df_nw['coord_x'].values
real_bhe_pos_y = real_df_nw['coord_y'].values

#%% virtual BHEs coordinates
#note: the sequence of the BHEs should be same as the BHEs sequence defined in the tespy model
#Read the coordinates from the BHE network csv
v_df_nw = read_csv('./base_module/pre/bhe_network.csv',
                    delimiter=';',
                    index_col=[0],
                    dtype={'data_index': str})

v_bhe_pos_x = v_df_nw['coord_x'].values
v_bhe_pos_y = v_df_nw['coord_y'].values

#%% selected bohrehole wall location point 
# each selected bhe has 4 around points, the calculated soil tempererature on the borhehole
# wall is the average temperature of the temperature at the 4 points in each timestep

#create local variables of the components label and connections label in network
localVars = locals()

#create the 4 reference points on each BHE's wall 
for i in range(bhe_num):
    localVars['bhe_'+ str(i) + '_wall_pos_x' ]  =  np.array([-br + v_bhe_pos_x[i],
                                                         0 + v_bhe_pos_x[i],
                                                         br + v_bhe_pos_x[i],
                                                         0 + v_bhe_pos_x[i]],dtype=float)
    localVars['bhe_'+ str(i) + '_wall_pos_y' ]   = np.array([0 + v_bhe_pos_y[i],
                                                             br + v_bhe_pos_y[i],
                                                             0 + v_bhe_pos_y[i],
                                                            -br + v_bhe_pos_y[i]],dtype=float)

#combine all reference points into one reference array, same with BHE coordinates 
#array sequence, every 4 points refer to one BHE
bhe_wall_pos_x = np.array([])
bhe_wall_pos_y = np.array([])
for i in range(bhe_num):
    bhe_wall_pos_x = np.concatenate((bhe_wall_pos_x, localVars['bhe_'+ str(i) + '_wall_pos_x' ]))
    bhe_wall_pos_y = np.concatenate((bhe_wall_pos_y, localVars['bhe_'+ str(i) + '_wall_pos_y' ]))

#the BHEs weights array. (The weights means the related BHE represents how many symmetry BHEs in the entire BHEs array)
bhe_weights = df_nw['weights'].values