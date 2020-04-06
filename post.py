# -*- coding: utf-8 -*- {}
'''
Post procedure to output the results csv files
'''
import pandas as pd

#%%output
def output_csv(BHE_num, timestep_tot,
               df_Tsoil, df_power, df_flowrate, df_Tin, df_Tout):
    # get the feature BHE names from the BHE_network.csv
    df_nw = pd.read_csv('./base_module/pre/bhe_network.csv',
                     delimiter=';',
                     index_col=[0],
                     dtype={'data_index': str})
    row_name = df_nw.index.values.tolist()

    # Bohrhole wall soil temperature:
    col_name = []
    for i in range(0, timestep_tot + 1):
        col_name.append('Time at ' + str(i) + 'step [K]')
    rt = pd.DataFrame(index=row_name, columns=col_name)
    rt.iloc[:] = df_Tsoil
    rt.to_csv('Result_Soil_T.csv')
    # BHE power:
    col_name = []
    for i in range(0, timestep_tot + 1):
        col_name.append('Time at ' + str(i) + 'step [W]')
    rt = pd.DataFrame(index=row_name, columns=col_name)
    rt.iloc[:] = df_power
    rt.to_csv('Result_BHE_power.csv')
    # BHE' pipe flowrate:
    col_name = []
    for i in range(0, timestep_tot + 1):
        col_name.append('Time at ' + str(i) + 'step [kg/s]')
    rt = pd.DataFrame(index=row_name, columns=col_name)
    rt.iloc[:] = df_flowrate
    rt.to_csv('Result_pipe_flowrate.csv')
    # BHE' inflow temperature:
    col_name = []
    for i in range(0, timestep_tot + 1):
        col_name.append('Time at ' + str(i) + 'step [K]')
    rt = pd.DataFrame(index=row_name, columns=col_name)
    rt.iloc[:] = df_Tin
    rt.to_csv('Result_Tin.csv')
    # BHE' outflow temperature:
    col_name = []
    for i in range(0, timestep_tot + 1):
        col_name.append('Time at ' + str(i) + 'step [K]')
    rt = pd.DataFrame(index=row_name, columns=col_name)
    rt.iloc[:] = df_Tout
    rt.to_csv('Result_Tout.csv')

    return

