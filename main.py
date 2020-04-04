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
time_tot = 90*24*60*60 #s
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
##soil: BHE' wall soil interface temperature in [K]
Result_df_soil = np.zeros((BHE_num, timestep_tot+1))
##BHE
#power in [W]
Result_df_BHE_power = np.zeros((BHE_num, timestep_tot+1))
#refrigerant flowrate in [kg/s]
Result_df_BHE_f_r = np.zeros((BHE_num, timestep_tot+1))
#inflow T in [K]
Result_df_fluid_in = np.zeros((BHE_num, timestep_tot+1))
#outflow T in [K]
Result_df_fluid_out = np.zeros((BHE_num, timestep_tot+1))

#parameters initialise
#soil
T0 = mod_ILS.T0
Result_df_soil[:,0] = T0
##BHE
#refrigerant flowrate
# flow rate in 4 significant digits to avoid deviation in tespy solver
Result_df_BHE_f_r += np.around(BHE_f_r_ini, decimals = 4)
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
        #update global dataframe BHE inflow and out flow
        first_step_Tin_and_Tout = mod_bhe.Type_1U_BHE_cal_singel(
                Result_df_BHE_power[0,step], T0, Result_df_BHE_f_r[0,step])
        Result_df_fluid_in[:,step] = first_step_Tin_and_Tout[0]
        Result_df_fluid_out[:,step] = first_step_Tin_and_Tout[1]
        #soil
        #global sourceterm dataframe in [W/m] in module ILS initialise ab timestep = 1
        mod_ILS.st_dataframe(1, Result_df_BHE_power[:,1]/BHE_length)
        #update global dataframe BHE' wall soil interface temperature
        Result_df_soil[:,step] = mod_ILS.ILS_solver(step)
        #sys time info output
        print('timestep %d took %.3f s' %(step, time.perf_counter() - time_step_start))
    else:
        #sys time info output      
        time_picard_iter_start = time.perf_counter()
        #uses BHE Tout array from last step as the initial Tout in the current network calculation
        Result_df_fluid_out[:,step] = Result_df_fluid_out[:,step - 1]
        #picard iteration until Tout achieves converge
        for i in range(max_iter):
            print('Picard iteration %d start:' % i)
            #convergence criterion initialise
            if_converge = False
            if i == max_iter - 1:
                raise Exception('The picard iteration could not achieve converge within %d steps' % max_iter)
                                
            #1st: TESPy solver
            time_mod_nw_start = time.perf_counter()
            if_sys_off, f_r, T_in = mod_nw.nw_solver(t,
                                      Result_df_fluid_out[:,step])
            # if system is shut off, no need to calculation for the current timestep
            if (if_sys_off):
                #set Tin same as it in the last timestep, Tout has been already set
                Result_df_fluid_in[:, step] = Result_df_fluid_in[:, step - 1]
                # flow rate dataframe do not need to update for save the calculation cost
                # set BHEs power to 0
                Result_df_BHE_power[:, step] = 0
                # break the iteration
                print('The system is shut off during the current time step')
                break
            #sys time info output
            print('Solve tespy network took %.3f s' %(time.perf_counter() - time_mod_nw_start))
            #update Tin and flowrate
            Result_df_fluid_in[:,step] = T_in
            Result_df_BHE_f_r[:,step] = np.around(f_r, decimals = 4) # flow rate in 4 significant digits
                                                                   # to avoid deviation in tespy solver

            #2nd: BHE solver
            time_mod_bhe_start = time.perf_counter()
            #record the Tout array from last iteration for the next converge check
            pre_BHEs_Tout = Result_df_fluid_out[:,step].copy()
            for j in range(BHE_num):#loop all BHEs
                #get the jth BHE's Tout and power from the current timestep
                #Tin, flowrate and Tsoil from last timestep
                #transfer the last step's flow rate of the BHE to determine if
                #the hydraulic coefficients need to be updated in the selected BHE
                cur_Tout_and_power = mod_bhe.Type_1U_BHE_cal(j,
                Result_df_fluid_in[j,step], 
                Result_df_soil[j,step - 1],
                Result_df_BHE_f_r[j,step],
                Result_df_BHE_f_r[j,step - 1])
                #get each BHE's Tout
                Result_df_fluid_out[j,step] = cur_Tout_and_power[0]
                #get each BHE's power
                Result_df_BHE_power[j,step] = cur_Tout_and_power[1]
            #sys time info output
            print('Solve BHE analytical solution for all BHEs took %.3f s' 
                  %(time.perf_counter() - time_mod_bhe_start))
            #determin if the Tout is converged
            #check norm if Tout array achieves the converge
            norm_delta_x = np.linalg.norm(abs((Result_df_fluid_out[:,step]) - pre_BHEs_Tout))
            norm_x = np.linalg.norm(abs(Result_df_fluid_out[:,step]))
            #sys time info output
            print('Convergence criterion: |dx|=%.3e, |x|=%.3e, |dx|/|x|=%.3e'
                  %(norm_delta_x, norm_x, norm_delta_x/norm_x))
            if (norm_delta_x/norm_x) < 1e-6:
                if_converge = True
            if (if_converge):
                #sys time info output
                print('Picard iteration achieves converge at %d steps, the total total iteration took %.3f s'
                      %(i, (time.perf_counter() - time_picard_iter_start)))
                break
            
        #3nd: ILS solver
        #update the global sourceterm dataframe in module ILS
        mod_ILS.st_dataframe(step, Result_df_BHE_power[:,step]/BHE_length)
        #update global dataframe soil temperature
        Result_df_soil[:,step] = mod_ILS.ILS_solver(step)
        
        #sys time info output
        print('Timestep %d took %.3f s' %(step, time.perf_counter() - time_step_start))
        
#solver time counter end
print('total solver execution took %d s' %(time.perf_counter() - time_solver_start))

#import post procedure
post.output_csv(BHE_num, timestep_tot,
                Result_df_soil, Result_df_BHE_power, Result_df_BHE_f_r,
                Result_df_fluid_in, Result_df_fluid_out)
#end main
print('The whole computation terminated on', time.ctime() )
