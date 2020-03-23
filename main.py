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
###TESPy initialise
# initialize the tespy model of the bhe network
# load path of network model:
# loading the TESPy model
project_dir = os.getcwd()
print("Project dir is: ", project_dir)
nw = load_network('./base_module/pre/tespy_nw')
# set if print the network iteration info
nw.set_attr(iterinfo=False)

# create bhe dataframe of the network system from bhe_network.csv
df = mod_nw.create_dataframe()
n_BHE = np.size(df.iloc[:, 0])

# create local variables of the components label and connections label in
# network
localVars = locals()
data_index = df.index.tolist()
for i in range(n_BHE):
    for c in nw.conns.index:
        # bhe inlet and outlet conns
        if c.t.label == data_index[i]:  # inlet conns of bhe
            localVars['inlet_BHE' + str(i + 1)] = c
        if c.s.label == data_index[i]:  # outlet conns of bhe
            localVars['outlet_BHE' + str(i + 1)] = c

# time depended consumer thermal demand
if mod_nw.switch_dyn_demand == 'on':
    # import the name of bus from the network csv file
    bus_name = pd.read_csv('./base_module/pre/tespy_nw/comps/bus.csv',
                    delimiter=';',
                    index_col=[0]).index[0]

# time depended flowrate
if mod_nw.switch_dyn_frate == 'on':
    # import the name of inlet connection from the network csv file
    inlet_name = pd.read_csv('./base_module/pre/tespy_nw/conn.csv',
                    delimiter=';',
                    index_col=[0]).iloc[0,0]
    for c in nw.conns.index:
        # bhe inflow conns
        if c.s.label == inlet_name:  # inlet conns of bhe
            localVars['inlet_name'] = c

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
'''    
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
'''    