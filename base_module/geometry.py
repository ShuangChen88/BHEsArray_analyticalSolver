# -*- coding: utf-8 -*- {}
'''
Geometry of the analytical model
bhe order: upperleft to upper right, upper to down
'''
import numpy as np

#%% User setting
#input geometry of the model

bhe_num = 3 #BHE number
BHE_wall_points_num = 4 #4 reference points on each BHE wall
adj = 6 # borehole adjacent distance in m
coord_center_point = [0,100]#the absolute coordinate of the center point of the array


#%% BHEs coordinates

bhe_pos_x = np.array([-adj + coord_center_point[0],
                       coord_center_point[0], 
                      adj + coord_center_point[0]],dtype=float)

bhe_pos_y = np.array([coord_center_point[1],
                      coord_center_point[1],
                      coord_center_point[1]],dtype=float)

#%% selected bohrehole wall location point 
# each selected bhe has 4 around points, the calculated soil tempererature on the borhehole
# wall is the average temperature of the temperature at the 4 points in each timestep
br = 0.075 # borehole raduis in m

#create local variables of the components label and connections label in network
localVars = locals()

#create the 4 reference points on each BHE's wall 
for i in range(bhe_num):
    localVars['bhe_'+ str(i) + '_wall_pos_x' ]  =  np.array([-br + bhe_pos_x[i],
                                                         0 + bhe_pos_x[i],
                                                         br + bhe_pos_x[i],
                                                         0 + bhe_pos_x[i]],dtype=float)
    localVars['bhe_'+ str(i) + '_wall_pos_y' ]   = np.array([0 + bhe_pos_y[i],
                                                             br + bhe_pos_y[i],
                                                             0 + bhe_pos_y[i],
                                                            -br + bhe_pos_y[i]],dtype=float)

#combine all reference points into one reference array, same with BHE coordinates 
#array sequence, every 4 points refer to one BHE
bhe_wall_pos_x = np.array([])
bhe_wall_pos_y = np.array([])
for i in range(bhe_num):
    bhe_wall_pos_x = np.concatenate((bhe_wall_pos_x, localVars['bhe_'+ str(i) + '_wall_pos_x' ]))
    bhe_wall_pos_y = np.concatenate((bhe_wall_pos_y, localVars['bhe_'+ str(i) + '_wall_pos_y' ]))