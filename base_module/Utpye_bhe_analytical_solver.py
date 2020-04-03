# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:09:59 2020

@author: caiw
"""
import numpy as np
from scipy import special as sp
import math

# function used to calculate the analytical solution
# of temperature distribution in an around Borehole 
# Heat Exchangers (BHE) proposed by Beier2014 #
# Richard A. Beier (2014) Transient heat transfer in a 
# U-tube borehole heat exchanger, 
# Applied Thermal Engineering 62: 256-266. 
# http://dx.doi.org/10.1016/j.applthermaleng.2013.09.014 #
# Author: Haibing Shao
# Email:  haibing(dot)shao(at)gmail(dot)com

#main parameters
time_tot = 90*24*60*60 #s
delta_t = 86400 #s
timestep_tot = int(time_tot/delta_t)

BHE_num = 3 #BHEs number
# set length of the borehole
Len = 50  # unit meter
# set input heat rate
# positive Q means heating up the soil
# negative Q means cooling down.
# Q0 = -25*50 # unit power watt
# BHE inlet temperature
#T_in = 10
# initial soil temperature
#T_s = 12

# property of the ground
k_s = 2.4  # heat conductivity of soil. Unit W/m/K
rho_s = 1120 #kg/m3
c_s = 1780 * rho_s  # heat capacity of soil. Unit J/m3/K
alpha_s = k_s / c_s

# property of the BHE
# borehole radius
r_b = 0.063  # unit m ni
# set dimensions and properties of U-tube pipe
d_po = 0.0167  # outer pipe diameter, unit m
d_pi = 0.0137  # inner pipe diameter, unit m
k_p = 0.39  # heat conductivity of pipe, unit W/m/K
# property of grout
c_g = 2190*1735.160  # heat capacity of grout. unit J/K/m3
k_g = 0.806  # heat conductivity of grout. unit W/m/K

# set circulating water rate
#parameter of the refrigerant
rho_f = 996 # kg/m3
#w = 0.1 * 4.428 / 1052  # unit m3/sec rho
w = 0.2 / rho_f  # unit m3/sec
# heat capacity of circulating water
c_f = 4179 * rho_f  # unit J/K/m3
k_f = 0.62863  # unit W/m/K

# calculate the convective film resistance inside pipe
v_pi = w / (math.pi / 4 * d_pi * d_pi)
mu_f = 0.0067418  # unit kg/m/s

# Reynolds number
Re = rho_f * v_pi * d_pi / mu_f

# estimate film coefficients by the Gnielinski correlation
Pr = mu_f / rho_f / alpha_s
# Churchill correlation for friction factor Eq.(B.4)
f = 1 / (1 / ((8 / Re) ** 10 + (Re / 36500) ** 20) ** 0.5 + (2.21 * math.log(Re / 7)) ** 10) ** 0.2
# Gnielinski correlation is used to estimate convective film coefficient Eq.(B.1)
Nu_Gl = (f / 8 * (Re - 1000) * Pr) / (1 + 12.7 * ((f / 8) ** 0.5) * (Pr ** (2 / 3) - 1))
Nu_lam = 4.364

# evaluate Nusselt number based on value of Re choose appropriate correlation
if Re < 2000:
    Nu = Nu_lam
else:
    Nu = Nu_Gl

h_pi = k_f * Nu / d_pi  # unit: W/m2/K

# evaluate shunt resistance for the heat transfer between fluids in two pipes
shank = r_b / 2  # m
# grout resistance fro shunt path Eq.(35)
R_gshunt = np.arccosh(2 * shank * shank / (d_po / 2) ** 2 - 1) / (2 * math.pi * k_g)
# Assign convective film coefficient on 1/2 of inside pipe wall area as corresponding to shunt resistance See Eq.(36)
# Assign resistance of 1/2 of pipe wall area contributing to shunt resistance See Eq. (36)
R_w1shunt = 1 / (1 / 2 * math.pi * d_pi * h_pi) + 1 / (math.pi * k_p) * math.log(d_po / d_pi)
# evaluate the shunt resistance 
R_12 = 2 * R_w1shunt + R_gshunt
# calculate dimensionless shunt conductance Eq.(9)
N_12 = Len / (R_12 * w * c_f)
# calculate resistance at grout/pipe interface
R_w1 = 1 / (2 * math.pi * k_p) * math.log(d_po / d_pi) + 1 / (math.pi * d_pi * h_pi)

R_b = 0.13  # unit m K/W
# R_w1 = 0.0001 # case 1: neglect this Resistance term
R_bmod = R_b - R_w1 / 2

# calculate equivalent radium Eq(31)
# r_eq = r_b / exp(2*math.pi*k_g*R_b) # This is case 2 on page 261 of Beier2014
r_eq = r_b / np.exp(2 * math.pi * k_g * R_bmod)  # We use case 1 on page 261 of Beier2014

# calculate dimensionless flow areas
A_D = (d_pi / 2 / r_eq) ** 2
# calculate conductance at grout/pipe interface for outer pipe Eq.(10)
N_w1 = Len / (w * c_f * R_w1)
N_w2 = N_w1
# calculate the ratio of grout an ground thermal conductivities 
kappa = k_g / k_s
# calculate ratio of volumetric heat capacities
H_g = c_g / c_s
# calculate ratio of water and soil heat capacities
H_f = c_f / c_s
# calculate dimensionless soil conductance
N_s = 2 * math.pi * k_s * Len / (w * c_f)
# calculate dimensionless borehole radius
r_Db = r_b / r_eq

# adjust volumetric heat capacity of grout to take into account actual volume of grout compared to volume of grout in model
Vol_actual = math.pi * (r_b * r_b - 2 * (d_po / 2) * (d_po / 2))
Vol_model = math.pi * (r_b * r_b - r_eq * r_eq)
H_g = H_g * Vol_actual / Vol_model

# to evaluate the dimensionless fluid inlet and outlet temperatures,  set zD = 0 
z_D = 0.0
# values of r_D and z_Dsand correspond to the location chosen to evaluate the temperature in the grout. 
# values of r_D and z_Dsand correspond to the location chosen to evaluate the temperature in the ground (sand).
# set the location to evaluate temperature in soil
# z_Dsand = 0.5
# r_D = r_Db
# r_Dsand = (r_b + 0.024) / r_eq  # this is problematic

# functions

def result(N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, z_D, timeDstart):
    # ns = 12
    # ncycles = 5
    # ndiv = 10
    # # zero out the data space
    # rt = np.zeros((ncycles * ndiv, 8))
    # # start calculating result
    # for i in range(ncycles):
    #     for j in range(ndiv):
    #         tval = timeD[j] * 10 ** (i)
    #         k = i * ndiv + j
    #         rt[k, 0] = k
    #         rt[k, 1] = tval
    #         rt[k, 2] = Stehfest_inv_Lap(F_1, tval, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, r_D, z_D, ns)
    #         rt[k, 3] = Stehfest_inv_Lap(F_2, tval, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, r_D, z_D, ns)
    #         rt[k, 4] = Stehfest_inv_Lap(FS_1, tval, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, r_D, z_Dsand, ns)
    #         rt[k, 5] = Stehfest_inv_Lap(FG_1, tval, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, r_Db, z_Dsand,ns)
    #         rt[k, 6] = Stehfest_inv_Lap(FS_2, tval, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, r_D, z_Dsand, ns)
    #         rt[k, 7] = Stehfest_inv_Lap(FG_2, tval, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, r_Db, z_Dsand,ns)
    ns = 12
    tval = timeDstart
    rt = np.zeros((3))
    rt[0] = tval # Fourier time
    rt[1] = Stehfest_inv_Lap(F_1, tval, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, z_D, ns) # inlet temperature at given time
    rt[2] = Stehfest_inv_Lap(F_2, tval, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, z_D, ns) # outlet temperature at given time
    return rt


# This function is performing the numerical inversion of Laplace transformation using the Stehfest algorithm;
def Stehfest_inv_Lap(funcHandle, time, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, z_D, ns):
    # ns is the number of terms to be used. Number of terms must be even. 
    nh = ns / 2
    v = np.zeros((ns, 1))
    # loop over each of vi
    for s in range(ns):
        fro = math.floor((s + 2) / 2)
        unt = min(s + 1, nh)
        if fro == 1 and unt == 1:
            v[s, 0] = v[s, 0] + (math.factorial(2 * 1)* 1 ** (ns / 2)) / (
                        math.factorial(ns / 2 - 1) * math.factorial(1) * math.factorial(1 - 1) * math.factorial(
                    s + 1 - 1) * math.factorial(2 * 1 - s - 1))
            v[s, 0] = v[s, 0] * (-1) ** (ns / 2 + s + 1)
        elif fro == 6 and unt == 6:
            v[s, 0] = v[s, 0] + (math.factorial(2 * 6)* 6 ** (ns / 2)) / (
                    math.factorial(ns / 2 - 6) * math.factorial(6) * math.factorial(6 - 1) * math.factorial(
                s + 1 - 6) * math.factorial(2 * 6 - s - 1))
            v[s, 0] = v[s, 0] * (-1) ** (ns / 2 + s + 1)
        else:
            for k in range(int(fro), int(unt) + 1):
                v[s, 0] = v[s, 0] + (math.factorial(2 * k)*k ** (ns / 2)) / (
                        math.factorial(ns / 2 - k) * math.factorial(k) * math.factorial(k - 1) * math.factorial(
                    s + 1 - k) * math.factorial(2 * k - s - 1))
            v[s, 0] = v[s, 0] * (-1) ** (ns / 2 + s + 1)
    # calculating inverse of laplace
    A = math.log(2) / time
    TD = float(0)
    for i in range(1, ns + 1):
        TD = TD + float(v[i - 1, 0]) * funcHandle(float(i * A), N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, z_D)
    TD = float(A) * TD
    rt = TD
    return rt


# This function correspons to T_D1(z,s); Eq.(A.17)
def F_1(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, z_D):
    rt = C_1(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db) * math.exp(
        a1(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db) * z_D) + C_2(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa,
                                                                              A_D, r_Db) * math.exp(
        a2(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db) * z_D)
    return rt


# This function correspons to T_D2(zD,s); Eq.(A.18)
def F_2(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, z_D):
    rt = C_3(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db) * math.exp(
        a1(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db) * z_D) + C_4(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa,
                                                                              A_D, r_Db) * math.exp(
        a2(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db) * z_D)
    return rt


# This function correspons to T_Dg1(rD,zD,s); Eq.(32) with Eq.(A.31) for B11(zD) and Eq.(A.8) for B21/B11)
# def FG_1(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, r_D, z_D):
#     numerator = F_1(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, r_D, z_D) * 2 * N_w1 / (kappa * N_s) * (
#             sp.i0(r_D * math.sqrt(H_g * s / kappa)) + B2B1(s, H_g, kappa, r_Db) * sp.k0(r_D * math.sqrt(s)))
#     denominator = 2 * N_w1 / (kappa * N_s) * (
#             sp.i0(math.sqrt(H_g * s / kappa)) + B2B1(s, H_g, kappa, r_Db) * sp.k0(math.sqrt(H_g * s / kappa))) + (
#                       -1) * math.sqrt(H_g * s / kappa) * sp.i1(math.sqrt(H_g * s / kappa)) + B2B1(s, H_g, kappa,
#                                                                                                   r_Db) * math.sqrt(
#         H_g * s / kappa) * sp.k1(math.sqrt(H_g * s / kappa))
#     rt = numerator / denominator
#     return rt


# This function correspons to T_Dg2(rD,zD,s); Eq.(32) with Eq.(A.31) for B11(zD) and Eq.(A.8) for B22/B12)
# def FG_2(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, r_D, z_D):
#     numerator = F_2(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, r_D, z_D) * 2 * N_w2 / (kappa * N_s) * (
#             sp.i0(r_D * math.sqrt(H_g * s / kappa)) + B2B1(s, H_g, kappa, r_Db) * sp.k0(r_D * math.sqrt(s)))
#     denominator = 2 * N_w2 / (kappa * N_s) * (
#             sp.i0(math.sqrt(H_g * s / kappa)) + B2B1(s, H_g, kappa, r_Db) * sp.k0(math.sqrt(H_g * s / kappa))) + (
#                       -1) * math.sqrt(H_g * s / kappa) * sp.i1(math.sqrt(H_g * s / kappa)) + B2B1(s, H_g, kappa,
#                                                                                                   r_Db) * math.sqrt(
#         H_g * s / kappa) * sp.k1(math.sqrt(H_g * s / kappa))
#     rt = numerator / denominator
#     return rt


# This function correspons to T_Ds1(rD,zD,s); Eq. (33) combined with Eq.(A.31) for B11(zD);
# def FS_1(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, r_D, z_D):
#     T_bar_D1 = F_1(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, r_D, z_D)
#     B11 = 2 * N_w1 / kappa / N_s * T_bar_D1
#     B11 = B11 / (2 * N_w1 / (kappa * N_s) * (sp.i0(math.sqrt(H_g * s / kappa)) + B2B1(s, H_g, kappa, r_Db) * sp.k0(
#         math.sqrt(H_g * s / kappa))) - math.sqrt(H_g * s / kappa) * sp.i1(math.sqrt(H_g * s / kappa)) + B2B1(s, H_g,
#                                                                                                              kappa,
#                                                                                                              r_Db) * sp.k1(
#         math.sqrt(H_g * s / kappa)))
#     rt = B11 * B4B1(s, H_g, kappa, r_Db) * sp.k0(r_D * math.sqrt(s))
#     return rt


# This function correspons to T_Ds2(rD,zD,s); Eq. (33) combined with Eq.(A.31) for B12(zD);
# def FS_2(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, r_D, z_D):
#     T_bar_D2 = F_2(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, r_D, z_D)
#     B12 = 2 * N_w1 / kappa / N_s * T_bar_D2
#     B12 = B12 / (2 * N_w1 / (kappa * N_s) * (sp.i0(math.sqrt(H_g * s / kappa)) + B2B1(s, H_g, kappa, r_Db) * sp.k0(
#         math.sqrt(H_g * s / kappa))) - math.sqrt(H_g * s / kappa) * sp.i1(math.sqrt(H_g * s / kappa)) + B2B1(s, H_g,
#                                                                                                              kappa,
#                                                                                                              r_Db) * sp.k1(
#         math.sqrt(H_g * s / kappa)))
#     rt = B12 * B4B1(s, H_g, kappa, r_Db) * sp.k0(r_D * math.sqrt(s))
#     return rt


# This function is according to Eq. (A.12b) It gives the same expression used for i=1 an i=2.
def C0(s, N_s, N_w, H_g, kappa, r_Db):
    rt = 1.0 - 1.0 / (1.0 - kappa * N_s / 2 / N_w * math.sqrt(H_g * s / kappa) * (
            sp.i1(math.sqrt(H_g * s / kappa)) - B2B1(s, H_g, kappa, r_Db) * sp.k1(math.sqrt(H_g * s / kappa))) / (
                              sp.i0(math.sqrt(H_g * s / kappa)) + B2B1(s, H_g, kappa, r_Db) * sp.k0(
                          math.sqrt(H_g * s / kappa))))
    return rt


# This function is according to Eq. (A.28)
def C_1(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db):
    rt = N_s / s / ((1 - delta_1(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db)) * (
            1 - math.exp(a1(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db)) / math.exp(
        a2(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db))))
    return rt


# This function is according to Eq. (A.29)
def C_2(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db):
    rt = ((N_s / s) - (
        (1 - delta_1(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db)) * C_1(s, N_s, N_w1, N_w2, N_12, H_g,
                                                                                    H_f, kappa, A_D, r_Db))) / (
                 1 - delta_2(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db))
    return rt


# This function is according to Eq. (A.22)
def C_3(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db):
    rt = delta_1(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db) * C_1(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa,
                                                                             A_D, r_Db)
    return rt


# This function is according to Eq. (A.23)
def C_4(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db):
    rt = delta_2(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db) * C_2(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa,
                                                                             A_D, r_Db)
    return rt


# This function is according to Eq. (A.24)
def delta_1(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db):
    rt = (a1(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db) + N_s * H_f * A_D * s / 2 + N_12 + N_w1 * C0(s, N_s,
                                                                                                                N_w1,
                                                                                                                H_g,
                                                                                                                kappa,
                                                                                                                r_Db)) / N_12
    return rt


# This function is according to Eq. (A.25)
def delta_2(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db):
    rt = (a2(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db) + N_s * H_f * A_D * s / 2 + N_12 + N_w1 * C0(s, N_s,
                                                                                                                N_w1,
                                                                                                                H_g,
                                                                                                                kappa,
                                                                                                                r_Db)) / N_12
    return rt


# This function is according to Eq. (A.21)
def gamma_local(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db):
    rt = -1.0 * (N_s * H_f * A_D * s / 2 + N_12 + N_w2 * C0(s, N_s, N_w2, H_g, kappa, r_Db)) * (
            N_s * H_f * A_D * s / 2 + N_12 + N_w1 * C0(s, N_s, N_w1, H_g, kappa, r_Db)) + N_12 * N_12
    return rt


# This function is according to Eq. (A.19a)
def a1(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db):
    rt = -0.5 * beta(s, N_s, N_w1, N_w2, H_g, H_f, kappa, A_D, r_Db) + 0.5 * math.sqrt(
        (beta(s, N_s, N_w1, N_w2, H_g, H_f, kappa, A_D, r_Db)) ** 2.0 - 4 * gamma_local(s, N_s, N_w1, N_w2, N_12, H_g,
                                                                                        H_f, kappa, A_D, r_Db))
    return rt


# This function is according to Eq. (A.19b)
def a2(s, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db):
    rt = -0.5 * beta(s, N_s, N_w1, N_w2, H_g, H_f, kappa, A_D, r_Db) - 0.5 * math.sqrt(
        (beta(s, N_s, N_w1, N_w2, H_g, H_f, kappa, A_D, r_Db)) ** 2.0 - 4 * gamma_local(s, N_s, N_w1, N_w2, N_12, H_g,
                                                                                        H_f, kappa, A_D, r_Db))
    return rt


# This function is according to Eq. (A.08)
# It gives the same expression used for i=1 an i=2
def B2B1(s, H_g, kappa, r_Db):
    numerator = kappa * math.sqrt(H_g * s / kappa) * sp.i1(r_Db * math.sqrt(H_g * s / kappa)) * sp.k0(
        r_Db * math.sqrt(s)) + math.sqrt(s) * sp.i0(r_Db * math.sqrt(H_g * s / kappa)) * sp.k1(r_Db * math.sqrt(s))
    denominator = kappa * math.sqrt(H_g * s / kappa) * sp.k1(r_Db * math.sqrt(H_g * s / kappa)) * sp.k0(
        r_Db * math.sqrt(s)) - math.sqrt(s) * sp.k0(r_Db * math.sqrt(H_g * s / kappa)) * sp.k1(r_Db * math.sqrt(s))
    rt = numerator / denominator
    return rt


# This function is according to Eq. (A.09)
def B4B1(s, H_g, kappa, r_Db):
    numerator = kappa * math.sqrt(H_g * s / kappa) * sp.i1(r_Db * math.sqrt(H_g * s / kappa)) * sp.k0(
        r_Db * math.sqrt(H_g * s / kappa)) + kappa * math.sqrt(H_g * s / kappa) * sp.i0(
        r_Db * math.sqrt(H_g * s / kappa)) * sp.k1(r_Db * math.sqrt(H_g * s / kappa))
    denominator = kappa * math.sqrt(H_g * s / kappa) * sp.k1(r_Db * math.sqrt(H_g * s / kappa)) * sp.k0(
        r_Db * math.sqrt(s)) - math.sqrt(s) * sp.k0(r_Db * math.sqrt(H_g * s / kappa)) * sp.k1(r_Db * math.sqrt(s))
    rt = numerator / denominator
    return rt


# This function is according to Eq. (A.20)
def beta(s, N_s, N_w1, N_w2, H_g, H_f, kappa, A_D, r_Db):
    rt = -0.5 * (N_s * H_f * A_D) * s - N_w2 * C0(s, N_s, N_w2, H_g, kappa,
                                                  r_Db) + 0.5 * N_s * H_f * A_D * s + N_w1 * C0(s, N_s, N_w1, H_g,
                                                                                                kappa, r_Db)
    return rt

def zresult(nz, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, z_D, timeDstart):
    ns = 12
    rt = np.zeros((nz, 4))
    for i in range(nz):
        z_D = (i+1)/nz
        tval = timeDstart
        rt[i,0] = i
        rt[i,1] = z_D
        rt[i,2] = Stehfest_inv_Lap(F_1, tval, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, z_D, ns)  # inlet temperature at given time
        rt[i,3] = Stehfest_inv_Lap(F_2, tval, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, z_D, ns)  # outlet temperature at given time
    return rt

# calculate fluid temperature distribution at the axial direction of borehole
def Type_1U_BHE_tDist(Power, T_soil, nz, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, z_D, timeDstart):
    zRD = zresult(nz, N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, z_D, timeDstart)
    # TD1 and TD2 are the dimensionless temperature
    zD = zRD[:, 1]
    zTD1 = zRD[:, 2]
    zTD2 = zRD[:, 3]
    H = np.zeros((nz))
    T1z = np.zeros((nz))
    T2z = np.zeros((nz))
    # evaluate the fluid temperatures at the axial direction from correspnding dimensionless variables
    for iz in range(nz):
        H[iz] = Len * zD[iz]
        T1z[iz] = Power * zTD1[iz] / (2 * math.pi * k_s * Len) + T_soil
        T2z[iz] = Power * zTD2[iz] / (2 * math.pi * k_s * Len) + T_soil

    return (H, T1z, T2z)

# interface functions to main.py calculation procedure for the first timestep
def Type_1U_BHE_cal_singel(Power, Tsoil, f_r):
    Tin = Power * global_coeff_Tin / (2 * math.pi * k_s * Len) + Tsoil
    Tout = Power * global_coeff_Tout / (2 * math.pi * k_s * Len) + Tsoil
    return (Tin,Tout)

# interface functions to main.py calculation procedure after the first timestep
def Type_1U_BHE_cal(BHE_id, T_in, T_soil, f_r_cur, f_r_pre):
    #set flow rate, hydraulic coefficient as global variables
    global w, global_coeff_Tin, global_coeff_Tout
    #determine if the hydraulic status in the BHE is changed
    if f_r_cur == f_r_pre:# hydraulic status unchanged
        Power = ((2 * math.pi * k_s * Len)*(T_in - T_soil))/global_coeff_Tin[BHE_id]
        Tout = Power * global_coeff_Tout[BHE_id] / (2 * math.pi * k_s * Len) + T_soil
    else:
        #update flow rate
        w = f_r_cur / rho_f
        #update the current BHE's hydraulic coefficient in the global_coeff_Tin
        global_coeff_Tin[BHE_id] = result(N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, 
                        A_D, r_Db, z_D, timeDstart)[1]
        global_coeff_Tout[BHE_id] = result(N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, 
                        A_D, r_Db, z_D, timeDstart)[2]
        #get the BHE's power and Tout
        Power = ((2 * math.pi * k_s * Len)*(T_in - T_soil))/global_coeff_Tin[BHE_id]
        Tout = Power * global_coeff_Tout[BHE_id] / (2 * math.pi * k_s * Len) + T_soil

    return (Tout, Power)

#%% main
# calculate the dimensionless time tD1
# ndiv = 10
# ncycles = 5
# timeDstart = 0.02
# xincr = 1.0 / ndiv
# iid = np.transpose([np.arange(ndiv + 1)])
# expon_id = xincr * iid
# timeDfirst = np.zeros((ndiv + 1, 1))
# for i in range(ndiv + 1):
#     timeDfirst[i] = 10 ** expon_id[i] * timeDstart
# Fourier time at given time step
timeDstart = (k_s*delta_t)/(c_s*r_eq**2)
# calculate dimensionless temperature as a function of dimensionless time
# TD1 = dimensionless temperature
# DD1 = dimensionless temperature derivative
# tD1 = dimensionless time (Fourier number)
# RD = result(N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, r_Dsand, z_D, z_Dsand, timeDfirst)
RD = result(N_s, N_w1, N_w2, N_12, H_g, H_f, kappa, A_D, r_Db, z_D, timeDstart)
# print(RD)
# convert from dimensionless Temperature to the real temperature
# TD1 and TD2 are the dimensionless temperature
TD1 = RD[1]
TD2 = RD[2]
# tD1 is the dimensionless time
tD1 = RD[0]
# =============================================================================
# TDS1 = RD[:, 4]
# TGS1 = RD[:, 5]
# TDS2 = RD[:, 6]
# TGS2 = RD[:, 7]
# =============================================================================
#t = (tD1 * c_s * (r_eq) ** 2 / k_s) * 3600  # this is in seconds
#T1 = Q0 * TD1 / (2 * math.pi * k_s * Len) + T_s
#T2 = Q0 * TD2 / (2 * math.pi * k_s * Len) + T_s
# evaluate the inlet & outlet temperatures and power from correspnding dimensionless variables
#t = (tD1 * c_s * (r_eq) ** 2 / k_s) * 3600  # this is in seconds
#T1 = Q0 * TD1 / (2 * math.pi * k_s * Len) + T_s
#T2 = Q0 * TD2 / (2 * math.pi * k_s * Len) + T_s
#Q =  ((2 * math.pi * k_s * Len)*(T_in - T_s))/TD1[0]
#T_out = T_in - Q/(c_f*w)
#print('For the given initial power', Q0, 'W, the inlet  temperature is', T1[0],'and outlet temperature is', T2[0])
#print('For the given inlet temperature', T_in, ', the BHE power is', Q , 'W, the outlet temperature is ',T_out)
#%% initialise TD1 and TD2 coefficient for inflow and outflow termperature of each BHEs
global_coeff_Tin = np.zeros(BHE_num) + TD1
global_coeff_Tout= np.zeros(BHE_num) + TD2