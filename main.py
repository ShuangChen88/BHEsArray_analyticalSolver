# -*- coding: utf-8 -*- {}
'''
2D BHEs array model for shifting behavior calculation procedure
Author: Shuang Chen
'''
import numpy as np
from scipy import special as sp
import os
import math
import pandas as pd
import time

#import base_module.ILS_analytical_solver
import base_module.bcs_tespy as mod_nw

#%%main parameters
time_tot = 12*30*24*60*60 #s
timestep = 3600 #s
step_tot = int(time_tot/timestep)

#soil
T0 = 15 #soil initial temperature
rho_s = 1120 #kg/m3
c_s = 1851.4 #J/kg*K
k_s = 1.8 #W/m*K
alpha = k_s/(rho_s*c_s) #m^2/s

#BHE
BHE_num = 3

#%%initialize
###data containers initialise
row_name = []
for i in range(1, BHE_num+1):
    row_name.append('BHE' + str(i))
col_name = []
for i in range(0, step_tot+1):
    col_name.append('Time at ' + str(i) + 'step')
##soil
Result_df_soil = pd.DataFrame(index=row_name, columns=col_name)
##BHE
Result_df_BHE_power = pd.DataFrame(index=row_name, columns=col_name)
##fluid
#inflow
Result_df_fluid_in = pd.DataFrame(index=row_name, columns=col_name)
#outflow
Result_df_fluid_out = pd.DataFrame(index=row_name, columns=col_name)

#parameter initialise
#soil
for i in range(BHE_num):
    Result_df_soil.iloc[i,0] = T0
#BHE
for i in range(BHE_num):
    Result_df_soil.iloc[i,0] = math.nan
    Result_df_soil.iloc[i,1] = mod_nw.consumer_demand(0)
    
#%%solve
#solver time counter start
time_solver_start = time.clock()
print('Solve process')
max_iter = 100

for t in range(step_tot):
    print('timestep = %d start' % t )
    for i in range(max_iter):
        print('Iteration %d :' % i)
        if i == max_iter - 1:
            raise Exception('The iteration could not achieve converge within %d steps' % max_iter)
    #base_module.bcs_tespy






#solver time counter end
time_solver_end = time.clock()
print('total solver execution took', (time_solver_end - time_solver_start)/60, 'min' )
print('The whole computation terminated on', time.ctime() )
    