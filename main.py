# -*- coding: utf-8 -*- {}
'''
2D BHEs array model for shifting behavior calculation procedure
Author: Shuang Chen
'''
import sys
print(sys.version)
import os
import numpy as np
import pandas as pd
import math
import time

from tespy.networks import load_network

#import base_module.ILS_analytical_solver
import base_module.Utpye_bhe_analytical_solver as mod_bhe
import base_module.bcs_tespy as mod_nw

#%%main parameters
time_tot = 12*30*24*60*60 #s
timestep = 86400 #s
step_tot = int(time_tot/timestep)

#soil
T0 = 15 #soil initial temperature
rho_s = 1120 #kg/m3
c_s = 1851.4 #J/kg*K
k_s = 1.8 #W/m*K
alpha = k_s/(rho_s*c_s) #m^2/s

#BHE
BHE_num = 3
BHE_f_r = 0.5 #kg/s

#%%initialize
###data containers initialise
row_name = []
for i in range(1, BHE_num+1):
    row_name.append('BHE' + str(i))
col_name = []
for i in range(0, step_tot+1):
    col_name.append('Time at ' + str(i) + 'step')
##soil
Result_df_soil = pd.DataFrame(np.zeros((BHE_num, step_tot+1)),index=row_name, columns=col_name)
##BHE
#BHE power is the driving force of the system, start from timestep = 1,
#all other results e.g. fluid and soil temperature are started from timestep = 1.
#power
Result_df_BHE_power = pd.DataFrame(np.zeros((BHE_num, step_tot+1)),index=row_name, columns=col_name)
#velocity
Result_df_BHE_f_r = pd.DataFrame(np.zeros((BHE_num, step_tot+1)),index=row_name, columns=col_name)
##fluid
#inflow T
Result_df_fluid_in = pd.DataFrame(index=row_name, columns=col_name)
#outflow T
Result_df_fluid_out = pd.DataFrame(index=row_name, columns=col_name)

#parameters initialise
#soil
Result_df_soil.iloc[:,0] = T0
#BHE
Result_df_BHE_power.iloc[:,1] = mod_nw.consumer_demand(0)/BHE_num
Result_df_BHE_f_r = Result_df_BHE_f_r + BHE_f_r

#'''    
#%%solve
#solver time counter starts
time_solver_start = time.clock()
print('Solve process')
max_iter = 100

for step in range(1, step_tot +1):
    t = step * timestep
    print('timestep = %d start' % step )
    #picard iteration until Tout achieves converge
    for i in range(max_iter):
        print('Iteration %d :' % i)
        if i == max_iter - 1:
            raise Exception('The iteration could not achieve converge within %d steps' % max_iter)
        else:
            #loop all BHEs
            for j in range(BHE_num):
                # equal power for all BHEs on the first thermal loading
                if t == 1:
                    Result_df_fluid_in.iloc[i,1] = mod_bhe.Type_1U_BHE_cal(
                            Result_df_BHE_power.iloc[i,1], T0)[0]
                    Result_df_fluid_out.iloc[i,1] = mod_bhe.Type_1U_BHE_cal(
                            Result_df_BHE_power.iloc[i,1], T0)[1]
                else:
                    break
                    #start with TESPy solver
                    
                    
                





#solver time counter end
time_solver_end = time.clock()
print('total solver execution took', (time_solver_end - time_solver_start)/60, 'min' )
print('The whole computation terminated on', time.ctime() )
#'''
