import os, re
import numpy as np 
import matplotlib.pyplot as plt

debug = True

def compute_rcr_parameters(P_min, P_max_diastolic, P_mean, Q_mean, ratio_prox_to_distal_resistors, decay_time, C_prefactor=1.0):

    tol = 1e-14

    # total resistance is determined by mean pressure and mean flow 
    R_total = P_mean / Q_mean 

    # ratio of resistors is constant 
    # resistors sum to total resistance 
    R_distal = R_total / (1.0 + ratio_prox_to_distal_resistors)
    R_proximal = R_total - R_distal

    assert abs(R_distal + R_proximal - R_total) < tol 

    # timescale for pressure decrease during aortic valve closure 
    C = -C_prefactor * decay_time / (R_distal * np.log(P_min/P_max_diastolic))

    return R_proximal, C, R_distal, R_total


if __name__== "__main__":

    MMHG_TO_CGS = 1333.22368

    standard_case = True 
    if standard_case: 
        high_pressure = False 
        low_pressure = False
        assert not (high_pressure and low_pressure)

        R_proximal =  77.0
        C = 0.001154
        R_distal =  1185.0
        print ("lit values:", )
        print ("R_proximal = ", R_proximal)
        print ("C = ", C)
        print ("R_distal = ", R_distal)
        print ("R_total = ", R_proximal + R_distal, "\n\n")


        beat_time = 0.8
        diastolic_time = .6 * beat_time
        systolic_time = beat_time - diastolic_time

        diastolic_time_fraction = diastolic_time / beat_time
        systolic_time_fraction = systolic_time / beat_time

        print ("diastolic_time_fraction = ", diastolic_time_fraction)
        print ("systolic_time_fraction = ", systolic_time_fraction)

        P_systolic = 120.0
        P_max_diastolic = 100.0
        P_min = 80.0

        P_diastolic_mean = P_min # 0.5 * (P_systolic + P_min)
        print ("P_diastolic_mean = ", P_diastolic_mean)

        P_mean = systolic_time_fraction*P_systolic + diastolic_time_fraction*(P_diastolic_mean)

        print ("P_mean = ", P_mean)

        if high_pressure: 
            P_mean *= 2 
            P_min *= 2
            P_max_diastolic = P_min + 2*20 

        if low_pressure: 
            P_mean *= 0.5 
            P_min *= 0.5
            P_max_diastolic = P_min + 0.5*20 

        # P_mean = 1.0433437500000039e+02
        
        # 5.6 L/min
        Q_mean = 5600.0 / 60 # ml/s 

        ratio_prox_to_distal_resistors = 77.0 / 1185.0 

        print("ratio_prox_to_distal_resistors = ", ratio_prox_to_distal_resistors)

        decay_time = diastolic_time

        C_prefactor = 1.0

        R_p_mmHg, C_mmHg, R_d_mmHg, R_total_mmHg = compute_rcr_parameters(P_min, P_max_diastolic, P_mean, Q_mean, ratio_prox_to_distal_resistors, decay_time, C_prefactor)

        name = 'aorta'

        print ("Values mmHg")
        print (name, ",\t", R_p_mmHg, ",\t", C_mmHg, ",\t", R_d_mmHg, ",\t", R_total_mmHg)
        print ("\n\n\n")

        P_min *= MMHG_TO_CGS
        P_max_diastolic *= MMHG_TO_CGS
        P_mean *= MMHG_TO_CGS

        R_p, C, R_d, R_total = compute_rcr_parameters(P_min, P_max_diastolic, P_mean, Q_mean, ratio_prox_to_distal_resistors, decay_time, C_prefactor)

        print ("Values CGS")
        print (name, ",\t", R_p, ",\t", C, ",\t", R_d, ",\t", R_total)

        print ("R_proximal = ", R_p)
        print ("C = ", C)
        print ("R_distal = ", R_d)


    fish_case = False
    if fish_case: 

        t_initial = 0.13361
        t_final = 0.457261


        beat_time = t_final - t_initial
        diastolic_time = .6 * beat_time
        systolic_time = beat_time - diastolic_time

        # diastolic_time_fraction = diastolic_time / beat_time
        # systolic_time_fraction = systolic_time / beat_time

        # print ("diastolic_time_fraction = ", diastolic_time_fraction)
        # print ("systolic_time_fraction = ", systolic_time_fraction)

        P_max_diastolic = 2.081957761631667
        P_min = 0.834874912254971
        P_mean = 1.391159596104874

        print ("P_mean = ", P_mean)

        # Q_mean = 5600.0 / 60 # ml/s 
        Q_mean = 8.218729433865492e-04

        ratio_prox_to_distal_resistors = 77.0 / 1185.0 

        print("ratio_prox_to_distal_resistors = ", ratio_prox_to_distal_resistors)

        decay_time = diastolic_time

        C_prefactor = 1.0

        R_p_mmHg, C_mmHg, R_d_mmHg, R_total_mmHg = compute_rcr_parameters(P_min, P_max_diastolic, P_mean, Q_mean, ratio_prox_to_distal_resistors, decay_time, C_prefactor)

        name = 'aorta'

        print ("Values mmHg")
        print (name, ",\t", R_p_mmHg, ",\t", C_mmHg, ",\t", R_d_mmHg, ",\t", R_total_mmHg)
        print ("\n\n\n")

        P_min *= MMHG_TO_CGS
        P_max_diastolic *= MMHG_TO_CGS
        P_mean *= MMHG_TO_CGS

        R_p, C, R_d, R_total = compute_rcr_parameters(P_min, P_max_diastolic, P_mean, Q_mean, ratio_prox_to_distal_resistors, decay_time, C_prefactor)

        print ("Values CGS")
        print (name, ",\t", R_p, ",\t", C, ",\t", R_d, ",\t", R_total)

        print ("R_proximal = ", R_p)
        print ("C = ", C)
        print ("R_distal = ", R_d)

