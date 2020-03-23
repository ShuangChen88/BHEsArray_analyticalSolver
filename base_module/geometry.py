# -*- coding: utf-8 -*- {}
'''
Geometry of the analytical model
bhe order: upperleft to upper right, upper to down
'''
import numpy as np

#%% input geometry of the model
#source term
bhe_num = 9
adj = 6 # borehole adjacent distance in m

bhe_pos_x = np.array([-adj,0,adj,
                      -adj,0,adj,
                      -adj,0,adj],dtype=float)

bhe_pos_y = np.array([adj,adj,adj,
                      0,0,0,
                     -adj,-adj,-adj],dtype=float)

#%% selected bohrehole wall location point 
# each selected bhe has 4 around points, the calculated soil tempererature on the borhehole
# wall is the average temperature of the temperature at the 4 points in each timestep
br = 0.075 # borehole raduis in m

#create local variables of the components label and connections label in network
localVars = locals()

# BHE 1
localVars['bhe_'+ str(0) + '_wall_pos_x' ]  = np.array([-br + (-adj),
                                                         0 + (-adj),
                                                         br + (-adj),
                                                         0 + (-adj),],dtype=float)
localVars['bhe_'+ str(0) + '_wall_pos_y' ]   = np.array([0 + (adj),
                                                         br + (adj),
                                                         0 + (adj),
                                                        -br + (adj),],dtype=float)

# BHE 2
localVars['bhe_'+ str(1) + '_wall_pos_x' ]  = np.array([-br + (-0),
                                                         0 + (-0),
                                                         br + (-0),
                                                         0 + (-0),],dtype=float)
localVars['bhe_'+ str(1) + '_wall_pos_y' ]   = np.array([0 + (adj),
                                                         br + (adj),
                                                         0 + (adj),
                                                        -br + (adj),],dtype=float)

# BHE 3
localVars['bhe_'+ str(2) + '_wall_pos_x' ]  = np.array([-br + (adj),
                                                         0 + (adj),
                                                         br + (adj),
                                                         0 + (adj),],dtype=float)
localVars['bhe_'+ str(2) + '_wall_pos_y' ]   = np.array([0 + (adj),
                                                         br + (adj),
                                                         0 + (adj),
                                                        -br + (adj),],dtype=float)

# BHE 4
localVars['bhe_'+ str(3) + '_wall_pos_x' ]  = np.array([-br + (-adj),
                                                         0 + (-adj),
                                                         br + (-adj),
                                                         0 + (-adj),],dtype=float)
localVars['bhe_'+ str(3) + '_wall_pos_y' ]   = np.array([0 + (0),
                                                         br + (0),
                                                         0 + (0),
                                                        -br + (0),],dtype=float)

# BHE 5
localVars['bhe_'+ str(4) + '_wall_pos_x' ]  = np.array([-br + (-0),
                                                         0 + (-0),
                                                         br + (-0),
                                                         0 + (-0),],dtype=float)
localVars['bhe_'+ str(4) + '_wall_pos_y' ]   = np.array([0 + (0),
                                                         br + (0),
                                                         0 + (0),
                                                        -br + (0),],dtype=float)

# BHE 6
localVars['bhe_'+ str(5) + '_wall_pos_x' ]  = np.array([-br + (adj),
                                                         0 + (adj),
                                                         br + (adj),
                                                         0 + (adj),],dtype=float)
localVars['bhe_'+ str(5) + '_wall_pos_y' ]   = np.array([0 + (0),
                                                         br + (0),
                                                         0 + (0),
                                                        -br + (0),],dtype=float)

# BHE 7
localVars['bhe_'+ str(6) + '_wall_pos_x' ]  = np.array([-br + (-adj),
                                                         0 + (-adj),
                                                         br + (-adj),
                                                         0 + (-adj),],dtype=float)
localVars['bhe_'+ str(6) + '_wall_pos_y' ]   = np.array([0 + (-adj),
                                                         br + (-adj),
                                                         0 + (-adj),
                                                        -br + (-adj),],dtype=float)

# BHE 8
localVars['bhe_'+ str(7) + '_wall_pos_x' ]  = np.array([-br + (-0),
                                                         0 + (-0),
                                                         br + (-0),
                                                         0 + (-0),],dtype=float)
localVars['bhe_'+ str(7) + '_wall_pos_y' ]   = np.array([0 + (-adj),
                                                         br + (-adj),
                                                         0 + (-adj),
                                                        -br + (-adj),],dtype=float)

# BHE 9
localVars['bhe_'+ str(8) + '_wall_pos_x' ]  = np.array([-br + (adj),
                                                         0 + (adj),
                                                         br + (adj),
                                                         0 + (adj),],dtype=float)
localVars['bhe_'+ str(8) + '_wall_pos_y' ]   = np.array([0 + (-adj),
                                                         br + (-adj),
                                                         0 + (-adj),
                                                        -br + (-adj),],dtype=float)