###
# Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
# Distributed under a Modified BSD License.
# See accompanying file LICENSE.txt or
# http://www.opengeosys.org/project/license
# TESPy network model module
###

import os 
import numpy as np
from pandas import read_csv
from tespy.networks import load_network

#%% User setting
# refrigerant parameters
#refrig_density = 992.92  # kg/m3
# switch for special boundary conditions
# 'on','off', switch of the function for dynamic thermal demand from consumer
switch_dyn_demand = 'on'
# 'on','off', switch of the function for dynamic flowrate in BHE
switch_dyn_frate = 'off'


# network status setting
def network_status(t):
    nw_status = 'on'
    # scenairo 1: month for closed network
    timerange_nw_off_month = [-9999]  # No month for closed network
    # t-1 to avoid the calculation problem at special time point,
    # e.g. t = 2592000.
    t_trans = int((t - 1) / 86400) + 1
    t_trans_month = t_trans
    if t_trans_month > 12:
        t_trans_month = t_trans - 12 * (int(t_trans / 12))
    if t_trans_month in timerange_nw_off_month:
        nw_status = 'off'

    # scenario 2: When dynamic system thermal load with 0 W
    if switch_dyn_demand == 'on':
        # consumer thermal load:
        cur_month_demand = consumer_demand(t)
        if cur_month_demand == 0:
            nw_status = 'off'
    return nw_status


# dynamic consumer thermal load
def consumer_demand(t):  # dynamic thermal demand from consumer
    # time conversion
    t_trans = int((t - 1) / 86400 /30 ) + 1
    if t_trans > 12:
        t_trans = t_trans - 12 * (int(t_trans / 12))
    # thermal demand in each month (assumed specific heat extraction rate*
    # length of BHE* number of BHE)
    month_demand = [
        -25 * 50 * 3, -0 * 50 * 3, -25 * 50 * 3, -25 * 50 * 3, -25 * 50 * 3,
        -25 * 50 * 3, -25 * 50 * 3, -25 * 50 * 3, -25 * 50 * 3, -25 * 50 * 3,
    ]
    return month_demand[t_trans - 1]


# dynamic hydraulic flow rate
def dyn_frate(t):  # dynamic flowrate in BHE
    # time conversion
    t_trans = int((t - 1) / 86400) + 1
    if t_trans > 12:
        t_trans = t_trans - 12 * (int(t_trans / 12))
    # flow rate in kg / s time curve in month
    month_frate = [-9999]
    return month_frate[t_trans - 1]


# End User setting
#%% create network dataframe
def create_dataframe():
    # return dataframe
    df_nw = read_csv('./base_module/pre/bhe_network.csv',
                        delimiter=';',
                        index_col=[0],
                        dtype={'data_index': str})
    return (df_nw)



# TESPy Thermal calculation process
def tespy_solver(t):
    # bhe network thermal re parametrization
    # if network exist dynamic flowrate
    if switch_dyn_frate == 'on':
        cur_frate = dyn_frate(t)
        localVars['inlet_name'].set_attr(m=cur_frate)
    # if network exist dynamic thermal load
    if switch_dyn_demand == 'on':
        # consumer thermal load:
        cur_month_demand = consumer_demand(t)
        # print('cur_month_demand', cur_month_demand)
        nw.busses[bus_name].set_attr(P=cur_month_demand)
    # T_out:
    for i in range(n_BHE):
        localVars['outlet_BHE' + str(i + 1)].set_attr(T=df_nw.loc[data_index[i],
                                                               'Tout_val'])
    # print('Tout=', df.loc[data_index[i], 'Tout_val'])
    # solving network
    nw.solve(mode='design')
    # get tespy solve result
    for i in range(n_BHE):
        # get flowrate # kg ^ 3 / s
        df_nw.loc[df_nw.index[i],
               'flowrate'] = localVars['inlet_BHE' +
                                      str(i + 1)].get_attr('m').val_SI
        # get Tin_val
        df_nw.loc[df_nw.index[i],
               'Tin_val'] = localVars['inlet_BHE' +
                                      str(i + 1)].get_attr('T').val_SI
    # print('Tin=', df_nw.loc[df_nw.index[i], 'Tin_val'])
    return (df_nw['flowrate'].tolist(), df_nw['Tin_val'].tolist())


#%% Tespy network model main solver function
def nw_solver(t, Tout_val):
    #system operation status control
    if_sys_off = False
    # network status:
    nw_status = network_status(t)
    # if network closed:
    if nw_status == 'off':
        if_sys_off = True
        #flowrate
        for i in range(n_BHE):
            df_nw.loc[df_nw.index[i], 'flowrate'] = 0
        cur_cal_f_r_val = df_nw['flowrate'].tolist()
        return (if_sys_off, cur_cal_f_r_val, Tout_val)
    else:
        # read Tout_val to dataframe
        for i in range(n_BHE):
            df_nw.loc[df_nw.index[i], 'Tout_val'] = Tout_val[i]
        # TESPy solver
        cur_cal_f_r_val, cur_cal_Tin_val = tespy_solver(t)
        # return to main
        return (if_sys_off, cur_cal_f_r_val, cur_cal_Tin_val)


#%%main
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
df_nw = create_dataframe()
n_BHE = np.size(df_nw.iloc[:, 0])
# create local variables of the components label and connections label in
# network
localVars = locals()
data_index = df_nw.index.tolist()
for i in range(n_BHE):
    for c in nw.conns.index:
        # bhe inlet and outlet conns
        if c.t.label == data_index[i]:  # inlet conns of bhe
            localVars['inlet_BHE' + str(i + 1)] = c
        if c.s.label == data_index[i]:  # outlet conns of bhe
            localVars['outlet_BHE' + str(i + 1)] = c

# time depended consumer thermal demand
if switch_dyn_demand == 'on':
    # import the name of bus from the network csv file
    bus_name = read_csv('./base_module/pre/tespy_nw/comps/bus.csv',
                    delimiter=';',
                    index_col=[0]).index[0]

# time depended flowrate
if switch_dyn_frate == 'on':
    # import the name of inlet connection from the network csv file
    inlet_name = read_csv('./base_module/pre/tespy_nw/conn.csv',
                    delimiter=';',
                    index_col=[0]).iloc[0,0]
    for c in nw.conns.index:
        # bhe inflow conns
        if c.s.label == inlet_name:  # inlet conns of bhe
            localVars['inlet_name'] = c
