import os, re
import numpy as np 
import matplotlib.pyplot as plt

debug = True

def compute_rcr_parameters(P_min, P_max, P_mean, Q_mean, R_total, ratio_prox_to_distal_resistors, decay_time, C_prefactor=1.0):

    tol = 1e-14 

    # ratio of resistors is constant 
    # resistors sum to total resistance 
    R_distal = R_total / (1.0 + ratio_prox_to_distal_resistors)
    R_proximal = R_total - R_distal

    assert abs(R_distal + R_proximal - R_total) < tol 

    # timescale for pressure decrease during aortic valve closure 
    C = -C_prefactor * decay_time / (R_distal * np.log(P_min/P_max))

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
        #print ("R_total = ", R_proximal + R_distal, "\n\n")

        HR = 146 # bpm
        MAP = 40 # mmHg
        CVP = 6 # mmHg

        P_systolic = 60 
        P_diastolic = 24

        LVEDV = 4.2 # ml
        LVESV = 1.8 # ml

        beat_time = 60/HR # s

        diastolic_time = (MAP - P_systolic)/(P_diastolic - P_systolic) * beat_time

        systolic_time = beat_time - diastolic_time

        diastolic_time_fraction = diastolic_time / beat_time
        systolic_time_fraction = systolic_time / beat_time

        print ("diastolic_time_fraction = ", diastolic_time_fraction)
        print ("systolic_time_fraction = ", systolic_time_fraction)

        P_diastolic_max = (P_systolic - P_diastolic)/2 + P_diastolic
        P_diastolic_min = P_diastolic

        P_mean = systolic_time_fraction*P_systolic + diastolic_time_fraction*(P_diastolic_min)

        print ("P_mean = ", P_mean)

        if high_pressure: 
            P_mean *= 2 
            P_min *= 2
            P_max = P_min + 2*20 

        if low_pressure: 
            P_mean *= 0.5 
            P_min *= 0.5
            P_max = P_min + 0.5*20 
        
        # CO
        SV = LVEDV - LVESV
        CO = 2 * (HR * SV)
        
        Q_inflow = CO / 60 # ml/s 

        QPQS = 4
 
        Q_aorta = 0.2 * Q_inflow
        Q_RPA = 2 * Q_aorta
        Q_LPA = 2 * Q_aorta

        # total resistance calculated from CO = (MAP-CVP)/R_total
        P_diff = MAP - CVP
        P_diff *= MMHG_TO_CGS
        R_total = P_diff/Q_inflow

        print("Q_inflow= ", Q_inflow, "Q_aorta= ", Q_aorta, "Q_RPA= ", Q_RPA, "Q_LPA= ", Q_LPA) 
        print("R_total= ", R_total)

        # resistances
        R_aorta = 5 * R_total
        R_LPA = 0.5 * R_aorta
        R_RPA = 0.5 * R_aorta

        ratio_prox_to_distal_resistors = 77.0 / 1185.0 

        print("ratio_prox_to_distal_resistors = ", ratio_prox_to_distal_resistors)

        decay_time = diastolic_time

        C_prefactor = 1.0

        R_p_aorta_mmHg, C_aorta_mmHg, R_d_aorta_mmHg, R_total_aorta_mmHg = compute_rcr_parameters(P_diastolic, P_systolic, P_mean, Q_aorta, R_aorta, ratio_prox_to_distal_resistors, decay_time, C_prefactor)
        R_p_LPA_mmHg, C_LPA_mmHg, R_d_LPA_mmHg, R_total_LPA_mmHg = compute_rcr_parameters(P_diastolic, P_systolic, P_mean, Q_LPA, R_LPA, ratio_prox_to_distal_resistors, decay_time, C_prefactor)
        R_p_RPA_mmHg, C_RPA_mmHg, R_d_RPA_mmHg, R_total_RPA_mmHg = compute_rcr_parameters(P_diastolic, P_systolic, P_mean, Q_RPA, R_RPA, ratio_prox_to_distal_resistors, decay_time, C_prefactor)

        name = 'aorta'

        print ("Values mmHg")
        print (name, ",\t", R_p_aorta_mmHg, ",\t", C_aorta_mmHg, ",\t", R_d_aorta_mmHg, ",\t", R_total_aorta_mmHg)
        print (name, ",\t", R_p_LPA_mmHg, ",\t", C_LPA_mmHg, ",\t", R_d_LPA_mmHg, ",\t", R_total_LPA_mmHg)
        print (name, ",\t", R_p_RPA_mmHg, ",\t", C_RPA_mmHg, ",\t", R_d_RPA_mmHg, ",\t", R_total_RPA_mmHg)
        print ("\n\n\n")

        P_diastolic *= MMHG_TO_CGS
        P_systolic *= MMHG_TO_CGS
        P_mean *= MMHG_TO_CGS

        R_p_aorta, C_aorta, R_d_aorta, R_total_aorta = compute_rcr_parameters(P_diastolic, P_systolic, P_mean, Q_aorta, R_aorta, ratio_prox_to_distal_resistors, decay_time, C_prefactor)
        R_p_LPA, C_LPA, R_d_LPA, R_total_LPA = compute_rcr_parameters(P_diastolic, P_systolic, P_mean, Q_LPA, R_LPA, ratio_prox_to_distal_resistors, decay_time, C_prefactor)
        R_p_RPA, C_RPA, R_d_RPA, R_total_RPA = compute_rcr_parameters(P_diastolic, P_systolic, P_mean, Q_RPA, R_RPA, ratio_prox_to_distal_resistors, decay_time, C_prefactor)

        print ("Values CGS")
        print (name, ",\t", R_p_aorta, ",\t", C_aorta, ",\t", R_d_aorta, ",\t", R_total_aorta)
        print (name, ",\t", R_p_LPA, ",\t", C_LPA, ",\t", R_d_LPA, ",\t", R_total_LPA)
        print (name, ",\t", R_p_RPA, ",\t", C_RPA, ",\t", R_d_RPA, ",\t", R_total_RPA)

        print ("R_proximal_aorta = ", R_p_aorta)
        print ("C_aorta = ", C_aorta)
        print ("R_distal_aorta = ", R_d_aorta)

        print ("R_proximal_LPA = ", R_p_LPA)
        print ("C_LPA = ", C_LPA)
        print ("R_distal_LPA = ", R_d_LPA)

        print ("R_proximal_RPA = ", R_p_RPA)
        print ("C_RPA = ", C_RPA)
        print ("R_distal_RPA = ", R_d_RPA)

