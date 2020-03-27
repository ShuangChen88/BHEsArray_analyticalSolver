# -*- coding: utf-8 -*- {}
'''
2D BHEs array model for shifting behavior calculation procedure
Author: Shuang Chen
'''
import sys
print(sys.version)
import numpy as np
import time

import base_module.ILS_analytical_solver as mod_ILS
import base_module.Utpye_bhe_analytical_solver as mod_bhe
import base_module.bcs_tespy as mod_nw
import post

#%% User setting
#main parameters
time_tot = 10*24*60*60 #s
delta_t = 86400 #s
timestep_tot = int(time_tot/delta_t)

#BHE
BHE_num = 3
BHE_length = 50 #m
#initial flowrate in each BHE, the global BHE flowrate curve is defined 
#in bcs_tespy.py or tespy model
BHE_f_r_ini = 0.2 #kg/s 

# End User setting
#%%initialize
###main data containers initialise,
#col: BHE id (:BHE_num)
#row: timestep (:timestep_tot+1), timestep = 0 means initial condition

#create the data containers
##soil: BHE' wall soil interface temperature
Result_df_soil = np.zeros((BHE_num, timestep_tot+1))
##BHE
#power
Result_df_BHE_power = np.zeros((BHE_num, timestep_tot+1))
#refrigerant flowrate
Result_df_BHE_f_r = np.zeros((BHE_num, timestep_tot+1))
#inflow T
Result_df_fluid_in = np.zeros((BHE_num, timestep_tot+1))
#outflow T
Result_df_fluid_out = np.zeros((BHE_num, timestep_tot+1))

#parameters initialise
#soil
T0 = mod_ILS.T0
Result_df_soil[:,0] = T0
##BHE
#refrigerant flowrate
Result_df_BHE_f_r += BHE_f_r_ini
#power
Result_df_BHE_power[:,0] = np.nan
##fluid
Result_df_fluid_in[:,0] = np.nan
#inflow T
Result_df_fluid_out[:,0] = np.nan

#%%solve
#solver time counter starts
time_solver_start = time.perf_counter()
print('Solve process')
max_iter = 100

for step in range(1, timestep_tot +1):
    t = step * delta_t
    time_step_start = time.perf_counter()
    print('timestep = %d start:' % step )
    # equal power for all BHEs on the first thermal loading
    if step == 1:
        #update global dataframe BHE power
        Result_df_BHE_power[:,step] = mod_nw.consumer_demand(0)/BHE_num
#        TODO: mod_bhe two calculation sort
        #update global dataframe BHE inflow and out flow
        first_step_Tin_and_Tout = mod_bhe.Type_1U_BHE_cal_singel(
                Result_df_BHE_power[0,step], T0, Result_df_BHE_f_r[0,step])
        Result_df_fluid_in[:,step] = first_step_Tin_and_Tout[0]
        Result_df_fluid_out[:,step] = first_step_Tin_and_Tout[1]
        #soil
        #global sourceterm dataframe in [W/m] in module ILS initialise ab timestep = 1
        mod_ILS.st_dataframe(1,0, Result_df_BHE_power[0,1]/BHE_length)
        #update global dataframe BHE' wall soil interface temperature
        Result_df_soil[:,step] = mod_ILS.ILS_solver(step)
        #sys time info output
        print('timestep %d took %.3f s' %(step, time.perf_counter() - time_step_start))
    else:
        #picard iteration until Tout achieves converge
        time_picard_iter_start = time.perf_counter()
        for i in range(max_iter):
            print('Picard iteration %d start:' % i)
            if i == max_iter - 1:
                raise Exception('The picard iteration could not achieve converge within %d steps' % max_iter)
                                
            #1st: TESPy solver
            time_mod_nw_start = time.perf_counter()
            if_converge, f_r, T_in = mod_nw.nw_solver(t,
                                      Result_df_fluid_in[:,step-1],
                                      Result_df_fluid_out[:,step-1])
            #sys time info output
            print('Solve tespy network took %.3f s' %(time.perf_counter() - time_mod_nw_start))
            #update flowrate and Tin
            Result_df_BHE_f_r[:,step] = f_r
            Result_df_fluid_in[:,step] = T_in
            
            #2nd: BHE solver
            time_mod_bhe_start = time.perf_counter()
            for j in range(BHE_num):#loop all BHEs
#                TODO: mod_bhe two calculation sort
                #get each BHE's Tout and power from the current timestep 
                #Tin, flowrate and Tsoil from last timestep
                cur_Tout_and_power = mod_bhe.Type_1U_BHE_cal(
                Result_df_fluid_in[j,step], 
                Result_df_soil[j,step-1],
                Result_df_BHE_f_r[j,step])
                #get each BHE's Tout
                Result_df_fluid_out[j,step] = cur_Tout_and_power[0]
                #get each BHE's power
                Result_df_BHE_power[j,step] = cur_Tout_and_power[1]
                #update the global sourceterm dataframe in module ILS
                mod_ILS.st_dataframe(step,j, Result_df_BHE_power[j,step]/BHE_length)
            #sys time info output
            print('Solve BHE analytical solution for all BHEs took %.3f s' 
                  %(time.perf_counter() - time_mod_bhe_start))
            #determin if the Tout is converged
            if (if_converge):
                #sys time info output
                print('Picard iteration achieves converge at %d steps, the total total iteration took %.3f s'
                      %(i, (time.perf_counter() - time_picard_iter_start)))
                break
            
        #3nd: ILS solver
        #update global dataframe soil temperature
        Result_df_soil[:,step] = mod_ILS.ILS_solver(step)
        
        #sys time info output
        print('Timestep %d took %.3f s' %(step, time.perf_counter() - time_step_start))
        
#solver time counter end
print('total solver execution took', (time.perf_counter() - time_solver_start)/60, 'min' )

#import post procedure
post.output_csv(BHE_num, timestep_tot,
                Result_df_soil, Result_df_BHE_power, Result_df_BHE_f_r,
                Result_df_fluid_in, Result_df_fluid_out)
#end main
print('The whole computation terminated on', time.ctime() )
