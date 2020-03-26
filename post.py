# -*- coding: utf-8 -*- {}
'''
Post procedure to output the results csv files
'''
import pandas as pd

#%%output 
def output_csv(BHE_num, timestep_tot,
               df_Tsoil, df_power, df_flowrate, df_Tin, df_Tout):
    row_name = []
    for i in range(1, BHE_num + 1):
        row_name.append('BHE' + str(i))
        
    #Bohrhole wall soil temperature:
    col_name = []
    for i in range(0, timestep_tot + 1):
        col_name.append('Time at ' + str(i) + 'step [K]')
    rt = pd.DataFrame(index=row_name, columns=col_name)
    rt.iloc[:] = df_Tsoil
    rt.to_csv('Result_Soil_T.csv')
    #BHE power:
    col_name = []
    for i in range(0, timestep_tot + 1):
        col_name.append('Time at ' + str(i) + 'step [W]')
    rt = pd.DataFrame(index=row_name, columns=col_name)
    rt.iloc[:] = df_power
    rt.to_csv('Result_BHE_power.csv')
    #BHE' pipe flowrate:
    col_name = []
    for i in range(0, timestep_tot + 1):
        col_name.append('Time at ' + str(i) + 'step [kg/s]')
    rt = pd.DataFrame(index=row_name, columns=col_name)
    rt.iloc[:] = df_flowrate
    rt.to_csv('Result_pipe_flowrate.csv')
    #BHE' inflow temperature:
    col_name = []
    for i in range(0, timestep_tot + 1):
        col_name.append('Time at ' + str(i) + 'step [K]')
    rt = pd.DataFrame(index=row_name, columns=col_name)
    rt.iloc[:] = df_Tin
    rt.to_csv('Result_Tin.csv')
    #BHE' outflow temperature:
    col_name = []
    for i in range(0, timestep_tot + 1):
        col_name.append('Time at ' + str(i) + 'step [K]')
    rt = pd.DataFrame(index=row_name, columns=col_name)
    rt.iloc[:] = df_Tout
    rt.to_csv('Result_Tout.csv')
    
    return

