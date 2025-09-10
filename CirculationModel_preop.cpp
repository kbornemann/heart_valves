// Filename: CirculationModel_preop.cpp
// Created on 20 Aug 2007 by Boyce Griffith

// Modified 2019, Alexander D. Kaiser

#include "CirculationModel_preop.h"
#include "pnpoly.h"
/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <tbox/RestartManager.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>

#include <Eigen/Dense>
using namespace Eigen;

namespace
{
// Name of output file.
static const string DATA_FILE_NAME = "bc_data.m";

} 

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

CirculationModel_preop::CirculationModel_preop(Pointer<Database> input_db, 
                                               const fourier_series_data *fourier_lvot,
                                               const fourier_series_data *fourier_rvot, 
					       string lvot_vertices_file_name, 
                                               string rvot_vertices_file_name,
 					       string aorta_vertices_file_name,
                                               string rpa_vertices_file_name,
                                               string lpa_vertices_file_name,
                                               const double  cycle_duration,
                                               const double  t_offset_bcs_unscaled, 
                                               const double  initial_time,
                                               double P_initial_aorta,
                                               double P_initial_rpa,
					       double P_initial_lpa,
                                               bool rcr_bcs_on,
                                               bool lvot_0D_on,
                                               bool rvot_0D_on, 
                                               double rcr_on_time)
    : 
      d_object_name("circ_model_preop"),  // constant name here  
      d_registered_for_restart(true),      // always true
      d_fourier_lvot(fourier_lvot),
      d_fourier_rvot(fourier_rvot), 
      d_cycle_duration(cycle_duration),
      d_t_offset_bcs_unscaled(t_offset_bcs_unscaled),
      d_current_idx_series(0),
      d_Q_lvot(0.0),
      d_Q_rvot(0.0),
      d_Q_aorta(0.0), 
      d_Q_rpa(0.0),
      d_Q_lpa(0.0),
      d_time_lvot(initial_time),
      d_time_rvot(initial_time),
      d_aorta_P_Wk(P_initial_aorta),
      d_rpa_P_Wk(P_initial_rpa),
      d_lpa_P_Wk(P_initial_lpa),
      d_area_lvot(0.0), 
      d_area_rvot(0.0),
      d_area_aorta(0.0),
      d_area_rpa(0.0),
      d_area_lpa (0.0),
      d_area_initialized(false), 
      d_rcr_bcs_on(rcr_bcs_on), 
      d_lvot_0D_on(lvot_0D_on),
      d_rvot_0D_on(rvot_0D_on),
      d_rcr_on_time(rcr_on_time)
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from the input and restart databases.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    
    if (d_rcr_bcs_on){
        if (input_db){
            // left and right equal for now
            
            d_aorta_R_proximal = input_db->getDouble("aorta_R_proximal");
            d_aorta_R_distal   = input_db->getDouble("aorta_R_distal");
            d_aorta_C          = input_db->getDouble("aorta_C");

            d_rpa_R_proximal = input_db->getDouble("rpa_R_proximal");
            d_rpa_R_distal   = input_db->getDouble("rpa_R_distal");
            d_rpa_C          = input_db->getDouble("rpa_C");

            d_lpa_R_proximal  = input_db->getDouble("lpa_R_proximal");
            d_lpa_R_distal    = input_db->getDouble("lpa_R_distal");
            d_lpa_C           = input_db->getDouble("lpa_C");

            pout << "input db got values:\n";
            pout << "aorta: R_proximal = " << d_aorta_R_proximal << "\tR_distal = " << d_aorta_R_distal << "\tC = " << d_aorta_C << "\n";
            pout << "rpa: R_proximal = " << d_rpa_R_proximal << "\tR_distal = " << d_rpa_R_distal << "\tC = " << d_rpa_C << "\n";
            pout << "lpa : R_proximal = " << d_lpa_R_proximal << "\tR_distal = " << d_lpa_R_distal << "\tC = " << d_lpa_C << "\n";
        }
        else {
            TBOX_ERROR("Must provide valid input_db");
        }
    }

    if (d_lvot_0D_on){
        d_lvot_0D = new lvot_0D_model(input_db, d_cycle_duration, d_rcr_on_time);
    } 
    else{
        d_lvot_0D = NULL;
        if (d_fourier_lvot == NULL){
            TBOX_ERROR("Must provide valid series for left ventricle if d_lvot_0D_on is off\n");
        }
    }

    if (d_rvot_0D_on){
        d_rvot_0D = new rvot_0D_model(input_db, d_cycle_duration, d_rcr_on_time);
    } else{
        d_rvot_0D = NULL;
        if (d_fourier_rvot == NULL){
            TBOX_ERROR("Must provide valid series for right ventricle if d_rvot_0D_on is off\n");
        }
    }

    double coord_normal, coord_normal_prev;
    double tol = 1.0e-2;

    // Find vertices from lvot file
    ifstream lvot_file(lvot_vertices_file_name.c_str(), ios::in);

    if(!lvot_file){
        TBOX_ERROR("LVOT file not found\n");
    }

    lvot_file >> d_n_pts_lvot;


    double *lvot_points_idx0 = new double[d_n_pts_lvot];
    double *lvot_points_idx1 = new double[d_n_pts_lvot];
    double *lvot_points_idx2 = new double[d_n_pts_lvot];

    int coords_flat[3] = {1, 1, 1};
    
    double x_prev, y_prev, z_prev;

    for (int i=0; i<d_n_pts_lvot; i++){

        lvot_file >> lvot_points_idx0[i];
        lvot_file >> lvot_points_idx1[i];
        lvot_file >> lvot_points_idx2[i];

        if (i>0){
            if (fabs(x_prev - lvot_points_idx0[i]) > tol){
                coords_flat[0] = 0;
            }
            if (fabs(y_prev - lvot_points_idx1[i]) > tol){
                coords_flat[1] = 0;
            }
            if (fabs(z_prev - lvot_points_idx2[i]) > tol){
                coords_flat[2] = 0;
            }
        }

        x_prev = lvot_points_idx0[i];
        y_prev = lvot_points_idx1[i];
        z_prev = lvot_points_idx2[i];

    }
    if (coords_flat[0]){

        if (coords_flat[1] || coords_flat[2]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_lvot_axis = 0;

        double tol_domain_x = (input_db->getDouble("X_HIGH") - input_db->getDouble("X_LOW"))/4.0;

        if (fabs(x_prev - input_db->getDouble("X_LOW")) < tol_domain_x){
            d_lvot_side = 0;
        }
        else if (fabs(x_prev - input_db->getDouble("X_HIGH")) < tol_domain_x){
            d_lvot_side = 1;
        }
        else{
            TBOX_ERROR("Flat coordinate not near boundary");
        }
        
        delete[] lvot_points_idx0;
        d_lvot_points_idx1 = lvot_points_idx1;
        d_lvot_points_idx2 = lvot_points_idx2;

    }
    else if (coords_flat[1]){
        if (coords_flat[0] || coords_flat[2]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_lvot_axis = 1;

        double tol_domain_y = (input_db->getDouble("Y_HIGH") - input_db->getDouble("Y_LOW"))/4.0;

        if (fabs(y_prev - input_db->getDouble("Y_LOW")) < tol_domain_y){
            d_lvot_side = 0;
        }
        else if (fabs(y_prev - input_db->getDouble("Y_HIGH")) < tol_domain_y){
            d_lvot_side = 1;
        }
        else{
            TBOX_ERROR("Flat coordinate not near boundary");
        }
       
        d_lvot_points_idx1 = lvot_points_idx0;
        delete[] lvot_points_idx1;
        d_lvot_points_idx2 = lvot_points_idx2;

    }
    else if (coords_flat[2]){
        if (coords_flat[0] || coords_flat[1]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_lvot_axis = 2;

        double tol_domain_z = (input_db->getDouble("Z_HIGH") - input_db->getDouble("Z_LOW"))/4.0;

        if (fabs(z_prev - input_db->getDouble("Z_LOW")) < tol_domain_z){
            d_lvot_side = 0;
        }
        else if (fabs(z_prev - input_db->getDouble("Z_HIGH")) < tol_domain_z){
            d_lvot_side = 1;
        }
        else{
            TBOX_ERROR("Flat coordinate not near boundary");
        }
       
        d_lvot_points_idx1 = lvot_points_idx0;
        d_lvot_points_idx2 = lvot_points_idx1;
        delete[] lvot_points_idx2;

    }
    else{
        TBOX_ERROR("Flat coordinate not found in boundary condition");
    }

    pout << "to lvot file close\n";
    pout << "found d_lvot_side = " << d_lvot_side << ", d_lvot_axis = " << d_lvot_axis << "\n";
    lvot_file.close();


    // Find vertices from rvot file
    ifstream rvot_file(rvot_vertices_file_name.c_str(), ios::in);

    if(!rvot_file){
        TBOX_ERROR("RVOT file not found\n");
    }

    rvot_file >> d_n_pts_rvot;


    double *rvot_points_idx0 = new double[d_n_pts_rvot];
    double *rvot_points_idx1 = new double[d_n_pts_rvot];
    double *rvot_points_idx2 = new double[d_n_pts_rvot];

    for (int i=0; i<d_n_pts_rvot; i++){

        rvot_file >> rvot_points_idx0[i];
        rvot_file >> rvot_points_idx1[i];
        rvot_file >> rvot_points_idx2[i];

        if (i>0){
            if (fabs(x_prev - rvot_points_idx0[i]) > tol){
                coords_flat[0] = 0;
            }
            if (fabs(y_prev - rvot_points_idx1[i]) > tol){
                coords_flat[1] = 0;
            }
            if (fabs(z_prev - rvot_points_idx2[i]) > tol){
                coords_flat[2] = 0;
            }
        }

        x_prev = rvot_points_idx0[i];
        y_prev = rvot_points_idx1[i];
        z_prev = rvot_points_idx2[i];

    }
    if (coords_flat[0]){

        if (coords_flat[1] || coords_flat[2]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_rvot_axis = 0;

        double tol_domain_x = (input_db->getDouble("X_HIGH") - input_db->getDouble("X_LOW"))/4.0;

        if (fabs(x_prev - input_db->getDouble("X_LOW")) < tol_domain_x){
            d_rvot_side = 0;
        }
        else if (fabs(x_prev - input_db->getDouble("X_HIGH")) < tol_domain_x){
            d_rvot_side = 1;
        }
        else{
            TBOX_ERROR("Flat coordinate not near boundary");
        }
        
        delete[] rvot_points_idx0;
        d_rvot_points_idx1 = rvot_points_idx1;
        d_rvot_points_idx2 = rvot_points_idx2;

    }
    else if (coords_flat[1]){
        if (coords_flat[0] || coords_flat[2]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_rvot_axis = 1;

        double tol_domain_y = (input_db->getDouble("Y_HIGH") - input_db->getDouble("Y_LOW"))/4.0;

        if (fabs(y_prev - input_db->getDouble("Y_LOW")) < tol_domain_y){
            d_rvot_side = 0;
        }
        else if (fabs(y_prev - input_db->getDouble("Y_HIGH")) < tol_domain_y){
            d_rvot_side = 1;
        }
        else{
            TBOX_ERROR("Flat coordinate not near boundary");
        }
       
        d_rvot_points_idx1 = rvot_points_idx0;
        delete[] rvot_points_idx1;
        d_rvot_points_idx2 = rvot_points_idx2;

    }
    else if (coords_flat[2]){
        if (coords_flat[0] || coords_flat[1]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_rvot_axis = 2;

        double tol_domain_z = (input_db->getDouble("Z_HIGH") - input_db->getDouble("Z_LOW"))/4.0;

        if (fabs(z_prev - input_db->getDouble("Z_LOW")) < tol_domain_z){
            d_rvot_side = 0;
        }
        else if (fabs(z_prev - input_db->getDouble("Z_HIGH")) < tol_domain_z){
            d_rvot_side = 1;
        }
        else{
            TBOX_ERROR("Flat coordinate not near boundary");
        }
       
        d_rvot_points_idx1 = rvot_points_idx0;
        d_rvot_points_idx2 = rvot_points_idx1;
        delete[] rvot_points_idx2;

    }
    else{
        TBOX_ERROR("Flat coordinate not found in boundary condition");
    }

    pout << "to rvot file close\n";
    pout << "found d_rvot_side = " << d_rvot_side << ", d_rvot_axis = " << d_rvot_axis << "\n";
    rvot_file.close();


    // Find vertices from aorta file
    ifstream aorta_file(aorta_vertices_file_name.c_str(), ios::in);

    if(!aorta_file){
        TBOX_ERROR("aorta file not found\n");
    }

    aorta_file >> d_n_pts_aorta;

    double *aorta_points_idx0 = new double[d_n_pts_aorta];
    double *aorta_points_idx1 = new double[d_n_pts_aorta];
    double *aorta_points_idx2 = new double[d_n_pts_aorta];

    coords_flat[0] = 1;
    coords_flat[1] = 1;
    coords_flat[2] = 1;

    for (int i=0; i<d_n_pts_aorta; i++){

        aorta_file >> aorta_points_idx0[i];
        aorta_file >> aorta_points_idx1[i];
        aorta_file >> aorta_points_idx2[i];

        if (i>0){
            if (fabs(x_prev - aorta_points_idx0[i]) > tol){
                coords_flat[0] = 0;
            }
            if (fabs(y_prev - aorta_points_idx1[i]) > tol){
                coords_flat[1] = 0;
            }
            if (fabs(z_prev - aorta_points_idx2[i]) > tol){
                coords_flat[2] = 0;
            }
        }

        x_prev = aorta_points_idx0[i];
        y_prev = aorta_points_idx1[i];
        z_prev = aorta_points_idx2[i];

       }
    if (coords_flat[0]){

        if (coords_flat[1] || coords_flat[2]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_aorta_axis = 0;

        double tol_domain_x = (input_db->getDouble("X_HIGH") - input_db->getDouble("X_LOW"))/4.0;

        if (fabs(x_prev - input_db->getDouble("X_LOW")) < tol_domain_x){
            d_aorta_side = 0;
        }
        else if (fabs(x_prev - input_db->getDouble("X_HIGH")) < tol_domain_x){
            d_aorta_side = 1;
        }
        else{
            TBOX_ERROR("Flat coordinate not near boundary");
        }
        
        delete[] aorta_points_idx0;
        d_aorta_points_idx1 = aorta_points_idx1;
        d_aorta_points_idx2 = aorta_points_idx2;

    }
    else if (coords_flat[1]){
        if (coords_flat[0] || coords_flat[2]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_aorta_axis = 1;

        double tol_domain_y = (input_db->getDouble("Y_HIGH") - input_db->getDouble("Y_LOW"))/4.0;

        if (fabs(y_prev - input_db->getDouble("Y_LOW")) < tol_domain_y){
            d_aorta_side = 0;
        }
        else if (fabs(y_prev - input_db->getDouble("Y_HIGH")) < tol_domain_y){
            d_aorta_side = 1;
        }
        else{
            TBOX_ERROR("Flat coordinate not near boundary");
        }
       
        d_aorta_points_idx1 = aorta_points_idx0;
        delete[] aorta_points_idx1;
        d_aorta_points_idx2 = aorta_points_idx2;

    }
    else if (coords_flat[2]){
        if (coords_flat[0] || coords_flat[1]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_aorta_axis = 2;

        double tol_domain_z = (input_db->getDouble("Z_HIGH") - input_db->getDouble("Z_LOW"))/4.0;

        if (fabs(z_prev - input_db->getDouble("Z_LOW")) < tol_domain_z){
            d_aorta_side = 0;
        }
        else if (fabs(z_prev - input_db->getDouble("Z_HIGH")) < tol_domain_z){
            d_aorta_side = 1;
        }
        else{
            TBOX_ERROR("Flat coordinate not near boundary");
        }
       
        d_aorta_points_idx1 = aorta_points_idx0;
        d_aorta_points_idx2 = aorta_points_idx1;
        delete[] aorta_points_idx2;

    }
    else{
        TBOX_ERROR("Aorta: Flat coordinate not found in boundary condition");
    }

    pout << "to aorta file close\n";
    pout << "found d_aorta_side = " << d_aorta_side << ", d_aorta_axis = " << d_aorta_axis << "\n";
    aorta_file.close();

    // Find vertices from rpa file
    ifstream rpa_file(rpa_vertices_file_name.c_str(), ios::in);

    if(!rpa_file){
        TBOX_ERROR("rpa file not found\n");
    }

    rpa_file >> d_n_pts_rpa;


    double *rpa_points_idx0 = new double[d_n_pts_rpa];
    double *rpa_points_idx1 = new double[d_n_pts_rpa];
    double *rpa_points_idx2 = new double[d_n_pts_rpa];

    for (int i=0; i<d_n_pts_rpa; i++){

        rpa_file >> rpa_points_idx0[i];
        rpa_file >> rpa_points_idx1[i];
        rpa_file >> rpa_points_idx2[i];

        if (i>0){
            if (fabs(x_prev - rpa_points_idx0[i]) > tol){
                coords_flat[0] = 0;
            }
            if (fabs(y_prev - rpa_points_idx1[i]) > tol){
                coords_flat[1] = 0;
            }
            if (fabs(z_prev - rpa_points_idx2[i]) > tol){
                coords_flat[2] = 0;
            }
        }

        x_prev = rpa_points_idx0[i];
        y_prev = rpa_points_idx1[i];
        z_prev = rpa_points_idx2[i];

    }
    if (coords_flat[0]){

        if (coords_flat[1] || coords_flat[2]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_rpa_axis = 0;

        double tol_domain_x = (input_db->getDouble("X_HIGH") - input_db->getDouble("X_LOW"))/4.0;

        if (fabs(x_prev - input_db->getDouble("X_LOW")) < tol_domain_x){
            d_rpa_side = 0;
        }
        else if (fabs(x_prev - input_db->getDouble("X_HIGH")) < tol_domain_x){
            d_rpa_side = 1;
        }
        else{
            TBOX_ERROR("Flat coordinate not near boundary");
        }
        
        delete[] rpa_points_idx0;
        d_rpa_points_idx1 = rpa_points_idx1;
        d_rpa_points_idx2 = rpa_points_idx2;

    }
    else if (coords_flat[1]){
        if (coords_flat[0] || coords_flat[2]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_rpa_axis = 1;

        double tol_domain_y = (input_db->getDouble("Y_HIGH") - input_db->getDouble("Y_LOW"))/4.0;

        if (fabs(y_prev - input_db->getDouble("Y_LOW")) < tol_domain_y){
            d_rpa_side = 0;
        }
        else if (fabs(y_prev - input_db->getDouble("Y_HIGH")) < tol_domain_y){
            d_rpa_side = 1;
        }
        else{
            TBOX_ERROR("Flat coordinate not near boundary");
        }
       
        d_rpa_points_idx1 = rpa_points_idx0;
        delete[] rpa_points_idx1;
        d_rpa_points_idx2 = rpa_points_idx2;

    }
    else if (coords_flat[2]){
        if (coords_flat[0] || coords_flat[1]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_rpa_axis = 2;

        double tol_domain_z = (input_db->getDouble("Z_HIGH") - input_db->getDouble("Z_LOW"))/4.0;

        if (fabs(z_prev - input_db->getDouble("Z_LOW")) < tol_domain_z){
            d_rpa_side = 0;
        }
        else if (fabs(z_prev - input_db->getDouble("Z_HIGH")) < tol_domain_z){
            d_rpa_side = 1;
        }
        else{
            TBOX_ERROR("Flat coordinate not near boundary");
        }
       
        d_rpa_points_idx1 = rpa_points_idx0;
        d_rpa_points_idx2 = rpa_points_idx1;
        delete[] rpa_points_idx2;

    }
    else{
        TBOX_ERROR("Rpa: Flat coordinate not found in boundary condition");
    }

    pout << "to rpa file close\n";
    pout << "found d_rpa_side = " << d_rpa_side << ", d_rpa_axis = " << d_rpa_axis << "\n";
    rpa_file.close();

    // Find vertices from lpa file
    ifstream lpa_file(lpa_vertices_file_name.c_str(), ios::in);

    if(!lpa_file){
        TBOX_ERROR("lpa file not found\n");
    }

    lpa_file >> d_n_pts_lpa;


    double *lpa_points_idx0 = new double[d_n_pts_lpa];
    double *lpa_points_idx1 = new double[d_n_pts_lpa];
    double *lpa_points_idx2 = new double[d_n_pts_lpa];

    for (int i=0; i<d_n_pts_lpa; i++){

        lpa_file >> lpa_points_idx0[i];
        lpa_file >> lpa_points_idx1[i];
        lpa_file >> lpa_points_idx2[i];

        if (i>0){
            if (fabs(x_prev - lpa_points_idx0[i]) > tol){
                coords_flat[0] = 0;
            }
            if (fabs(y_prev - lpa_points_idx1[i]) > tol){
                coords_flat[1] = 0;
            }
            if (fabs(z_prev - lpa_points_idx2[i]) > tol){
                coords_flat[2] = 0;
            }
        }

        x_prev = lpa_points_idx0[i];
        y_prev = lpa_points_idx1[i];
        z_prev = lpa_points_idx2[i];

    }
    if (coords_flat[0]){

        if (coords_flat[1] || coords_flat[2]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_lpa_axis = 0;

        double tol_domain_x = (input_db->getDouble("X_HIGH") - input_db->getDouble("X_LOW"))/4.0;

        if (fabs(x_prev - input_db->getDouble("X_LOW")) < tol_domain_x){
            d_lpa_side = 0;
        }
        else if (fabs(x_prev - input_db->getDouble("X_HIGH")) < tol_domain_x){
            d_lpa_side = 1;
        }
        else{
            TBOX_ERROR("Flat coordinate not near boundary");
        }
        
        delete[] lpa_points_idx0;
        d_lpa_points_idx1 = lpa_points_idx1;
        d_lpa_points_idx2 = lpa_points_idx2;

    }
    else if (coords_flat[1]){
        if (coords_flat[0] || coords_flat[2]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_lpa_axis = 1;

        double tol_domain_y = (input_db->getDouble("Y_HIGH") - input_db->getDouble("Y_LOW"))/4.0;

        if (fabs(y_prev - input_db->getDouble("Y_LOW")) < tol_domain_y){
            d_lpa_side = 0;
        }
        else if (fabs(y_prev - input_db->getDouble("Y_HIGH")) < tol_domain_y){
            d_lpa_side = 1;
        }
        else{
            TBOX_ERROR("Flat coordinate not near boundary");
        }
       
        d_lpa_points_idx1 = lpa_points_idx0;
        delete[] lpa_points_idx1;
        d_lpa_points_idx2 = lpa_points_idx2;

    }
    else if (coords_flat[2]){
        if (coords_flat[0] || coords_flat[1]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_lpa_axis = 2;

        double tol_domain_z = (input_db->getDouble("Z_HIGH") - input_db->getDouble("Z_LOW"))/4.0;

        if (fabs(z_prev - input_db->getDouble("Z_LOW")) < tol_domain_z){
            d_lpa_side = 0;
        }
        else if (fabs(z_prev - input_db->getDouble("Z_HIGH")) < tol_domain_z){
            d_lpa_side = 1;
        }
        else{
            TBOX_ERROR("Flat coordinate not near boundary");
        }
       
        d_lpa_points_idx1 = lpa_points_idx0;
        d_lpa_points_idx2 = lpa_points_idx1;
        delete[] lpa_points_idx2;

    }
    else{
        TBOX_ERROR("Lpa: Flat coordinate not found in boundary condition");
    }

    pout << "to lpa file close\n";
    pout << "found d_lpa_side = " << d_lpa_side << ", d_lpa_axis = " << d_lpa_axis << "\n";
    lpa_file.close();


    if (!from_restart){
        if ((d_lvot_0D_on) && (d_rvot_0D_on)){
            d_lvot_P = d_lvot_0D->d_P_lvot;
            d_rvot_P = d_rvot_0D->d_P_rvot;
        }

        else if ((d_lvot_0D_on) && (!d_rvot_0D_on)){
            d_lvot_P = d_lvot_0D->d_P_lvot;
            d_rvot_P = d_lvot_0D->d_P_lvot;
        }

        else if ((d_rvot_0D_on) && (!d_lvot_0D_on)){
            d_rvot_P = d_rvot_0D->d_P_rvot;
            d_lvot_P = d_rvot_0D->d_P_rvot;
        }
        else {
            d_lvot_P = MMHG_TO_CGS * d_fourier_lvot->values[0];
            d_rvot_P = MMHG_TO_CGS * d_fourier_rvot->values[0];
        }
        d_aorta_P = P_initial_aorta;
        d_rpa_P = P_initial_rpa;
        d_lpa_P = P_initial_lpa;

        d_P_min_linear_interp = 0.0;
    }
    d_p_equal_fraction = 0.1;

    pout << "passed constructor\n";
    pout << "intial lv pressure = " << d_lvot_P << ", initial rv pressure = " << d_rvot_P << ", initial aorta pressure = " << d_aorta_P << ", initial rpa pressure = " << d_rpa_P << ", initial lpa pressure = " << d_lpa_P << "\n";

    return;
} // CirculationModel

CirculationModel_preop::~CirculationModel_preop()
{
    return;
} // ~CirculationModel_preop


void CirculationModel_preop::advanceTimeDependentData(const double dt,
                                                        const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                        const int U_idx,
                                                        const int /*P_idx*/,
                                                        const int /*wgt_cc_idx*/,
                                                        const int wgt_sc_idx)
{
    // Compute the mean flow rates in the vicinity of the inflow and outflow
    // boundaries.
 
    double Q_lvot_local = 0.0;   
    double Q_rvot_local = 0.0;
    double Q_aorta_local = 0.0; 
    double Q_rpa_local = 0.0; 
    double Q_lpa_local = 0.0; 

    double area_lvot_local = 0.0;
    double area_rvot_local = 0.0;
    double area_aorta_local = 0.0; 
    double area_rpa_local = 0.0; 
    double area_lpa_local = 0.0; 

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            if (pgeom->getTouchesRegularBoundary())
            {
                Pointer<SideData<NDIM, double> > U_data = patch->getPatchData(U_idx);
                Pointer<SideData<NDIM, double> > wgt_sc_data = patch->getPatchData(wgt_sc_idx);
                const Box<NDIM>& patch_box = patch->getBox();
                const double* const x_lower = pgeom->getXLower();
                const double* const dx = pgeom->getDx();
                double dV = 1.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    dV *= dx[d];
                }

                for(int axis=0; axis<3; axis++)
                {
                    for(int side=0; side<2; side++)
                    {
                        const bool is_lower = (side == 0);
                        if (pgeom->getTouchesRegularBoundary(axis, side))
                        {
                            
                            IBTK::Vector n;
                            for (int d = 0; d < NDIM; ++d)
                            {
                                n[d] = axis == d ? (is_lower ? -1.0 : +1.0) : 0.0;
                            }
                            Box<NDIM> side_box = patch_box;
                            if (is_lower)
                            {
                                side_box.lower(axis) = patch_box.lower(axis);
                                side_box.upper(axis) = patch_box.lower(axis);
                            }
                            else
                            {
                                side_box.lower(axis) = patch_box.upper(axis) + 1;
                                side_box.upper(axis) = patch_box.upper(axis) + 1;
                            }
                            for (Box<NDIM>::Iterator b(side_box); b; b++)
                            {
                                const SAMRAI::hier::Index<NDIM>& i = b();

                                double X[NDIM];
                                for (int d = 0; d < NDIM; ++d)
                                {
                                    X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_box.lower(d)) + (d == axis ? 0.0 : 0.5));
                                }

                                double X_in_plane_1 = 0.0; 
                                double X_in_plane_2 = 0.0; 
                                if (axis == 0)
                                {
                                    X_in_plane_1 = X[1]; 
                                    X_in_plane_2 = X[2]; 
                                }
                                else if (axis == 1)
                                {
                                    X_in_plane_1 = X[0]; 
                                    X_in_plane_2 = X[2]; 
                                }
                                else if (axis == 2)
                                {
                                    X_in_plane_1 = X[0]; 
                                    X_in_plane_2 = X[1]; 
                                }
                                else{
                                    TBOX_ERROR("Invalid value of axis\n"); 
                                }

                                const int in_lvot  = this->point_in_lvot(X_in_plane_1, X_in_plane_2, axis, side);
                                const int in_rvot  = this->point_in_rvot(X_in_plane_1, X_in_plane_2, axis, side);
                                const int in_aorta         = this->point_in_aorta       (X_in_plane_1, X_in_plane_2, axis, side);
                                const int in_rpa         = this->point_in_rpa       (X_in_plane_1, X_in_plane_2, axis, side);
                                const int in_lpa          = this->point_in_lpa        (X_in_plane_1, X_in_plane_2, axis, side);

                                if (in_lvot)
                                {
                                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                                    if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                                    {
                                        double dA = dV / dx[axis];
                                        Q_lvot_local += (*U_data)(i_s)* n[axis] * dA;

                                        if (!d_area_initialized){
                                            area_lvot_local += dA;
                                        }

                                    }
                                }

                                if (in_rvot)
                                {
                                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                                    if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                                    {
                                        double dA = dV / dx[axis];
                                        Q_rvot_local += (*U_data)(i_s)* n[axis] * dA;

                                        if (!d_area_initialized){
                                            area_rvot_local += dA;
                                        }

                                    }
                                }

                                if (in_aorta)
                                {
                                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                                    if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                                    {
                                        double dA = dV / dx[axis];
                                        Q_aorta_local += (*U_data)(i_s) * n[axis] * dA;

                                        if (!d_area_initialized){
                                            area_aorta_local += dA;
                                        }

                                    }
                                }

                                if (in_rpa)
                                {
                                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                                    if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                                    {
                                        double dA = dV / dx[axis];
                                        Q_rpa_local += (*U_data)(i_s) * n[axis] * dA;

                                        if (!d_area_initialized){
                                            area_rpa_local += dA;
                                        }

                                    }
                                }

                                if (in_lpa)
                                {
                                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                                    if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                                    {
                                        double dA = dV / dx[axis];
                                        Q_lpa_local += (*U_data)(i_s) * n[axis] * dA;

                                        if (!d_area_initialized){
                                            area_lpa_local += dA;
                                        }

                                    }
                                }

                            }
                        }
                    }
                }
            }
        }
    }

    d_Q_lvot = SAMRAI_MPI::sumReduction(Q_lvot_local);
    d_Q_rvot = SAMRAI_MPI::sumReduction(Q_rvot_local);
    d_Q_aorta        = SAMRAI_MPI::sumReduction(Q_aorta_local);
    d_Q_rpa        = SAMRAI_MPI::sumReduction(Q_rpa_local);
    d_Q_lpa         = SAMRAI_MPI::sumReduction(Q_lpa_local);

    if (!d_area_initialized){
        d_area_lvot = SAMRAI_MPI::sumReduction(area_lvot_local);
        d_area_rvot = SAMRAI_MPI::sumReduction(area_rvot_local);
        d_area_aorta        = SAMRAI_MPI::sumReduction(area_aorta_local);   
        d_area_rpa        = SAMRAI_MPI::sumReduction(area_rpa_local);  
        d_area_lpa         = SAMRAI_MPI::sumReduction(area_lpa_local);  
        d_area_initialized = true;       
    }

    if (d_rcr_bcs_on){ 
            d_aorta_P_Wk = ((d_aorta_C / dt) * d_aorta_P_Wk + d_Q_aorta) / (d_aorta_C / dt + 1.0 / d_aorta_R_distal);
            d_aorta_P = d_aorta_P_Wk + d_aorta_R_proximal * d_Q_aorta;

            d_rpa_P_Wk = ((d_rpa_C / dt) * d_rpa_P_Wk + d_Q_rpa) / (d_rpa_C / dt + 1.0 / d_rpa_R_distal);
            d_rpa_P = d_rpa_P_Wk + d_rpa_R_proximal * d_Q_rpa;

            d_lpa_P_Wk = ((d_lpa_C / dt) * d_lpa_P_Wk + d_Q_lpa) / (d_lpa_C / dt + 1.0 / d_lpa_R_distal);
            d_lpa_P = d_lpa_P_Wk + d_lpa_R_proximal * d_Q_lpa;
    }


    d_time += dt;

    if ((d_lvot_0D_on) && (d_rvot_0D_on)){
        d_lvot_0D->advanceTimeDependentData(dt, d_time, -d_Q_lvot);
        d_lvot_P = d_lvot_0D->d_P_lvot;

        d_rvot_0D->advanceTimeDependentData(dt, d_time, -d_Q_rvot);
        d_rvot_P = d_rvot_0D->d_P_rvot;
    }

    else if ((d_lvot_0D_on) && (!d_rvot_0D_on)){
        d_lvot_0D->advanceTimeDependentData(dt, d_time, -d_Q_lvot);
        d_lvot_P = d_lvot_0D->d_P_lvot;
        d_rvot_P = d_lvot_0D->d_P_lvot;
    }

    else if ((d_rvot_0D_on) && (!d_lvot_0D_on)){
        d_rvot_0D->advanceTimeDependentData(dt, d_time, -d_Q_rvot);
        d_rvot_P = d_rvot_0D->d_P_rvot;
        d_lvot_P = d_rvot_0D->d_P_rvot;
    }

    else{
        
        // compute which index in the Fourier series we need here 
        // always use a time in current cycle 
        double t_reduced = d_time - d_cycle_duration * floor(d_time/d_cycle_duration);

        // fourier series has its own period, scale to that 
        double t_scaled = t_reduced * (d_fourier_lvot->L  / d_cycle_duration); 

        // start offset some arbitrary time in the cardiac cycle, but this is relative to the series length 
        double t_scaled_offset = t_scaled + d_t_offset_bcs_unscaled;

        // Fourier data here
        // index without periodicity 
        unsigned int k = (unsigned int) floor(t_scaled_offset / (d_fourier_lvot->dt));
    
        // // take periodic reduction
        d_current_idx_series = k % (d_fourier_lvot->N_times);

        d_lvot_P = MMHG_TO_CGS * d_fourier_lvot->values[d_current_idx_series];
        d_rvot_P = MMHG_TO_CGS * d_fourier_rvot->values[d_current_idx_series];       

 
    }

    writeDataFile(); 

} // advanceTimeDependentData

void CirculationModel_preop::set_Q_valve(double Q_valve){
    d_Q_valve = Q_valve; 
}

void
CirculationModel_preop::putToDatabase(Pointer<Database> db)
{

    db->putInteger("d_current_idx_series", d_current_idx_series);
    db->putDouble("d_Q_lvot", d_Q_lvot);
    db->putDouble("d_Q_rvot", d_Q_rvot);
    db->putDouble("d_Q_aorta", d_Q_aorta); 
    db->putDouble("d_Q_rpa", d_Q_rpa);
    db->putDouble("d_Q_lpa", d_Q_lpa);
    db->putDouble("d_Q_valve", d_Q_valve);
    db->putDouble("d_lvot_P", d_lvot_P);
    db->putDouble("d_rvot_P", d_rvot_P);
    db->putDouble("d_aorta_P", d_aorta_P);
    db->putDouble("d_aorta_P_Wk", d_aorta_P_Wk);
    db->putDouble("d_rpa_P", d_rpa_P);
    db->putDouble("d_rpa_P_Wk", d_rpa_P_Wk);
    db->putDouble("d_lpa_P",d_lpa_P);
    db->putDouble("d_lpa_P_Wk",d_lpa_P_Wk);
    db->putDouble("d_time", d_time); 
    db->putBool("d_rcr_bcs_on", d_rcr_bcs_on);
    db->putBool("d_lvot_0D_on", d_lvot_0D_on);
    db->putBool("d_rvot_0D_on", d_rvot_0D_on);
    db->putDouble("d_P_min_linear_interp", d_P_min_linear_interp);
    db->putDouble("d_rcr_on_time", d_rcr_on_time); 
    return; 
} // putToDatabase

void CirculationModel_preop::print_summary(){

    double P_lvot = d_lvot_P / MMHG_TO_CGS;
    double P_rvot = d_rvot_P / MMHG_TO_CGS;
    double P_aorta = d_aorta_P / MMHG_TO_CGS; 
    double P_rpa = d_rpa_P / MMHG_TO_CGS; 
    double P_lpa = d_lpa_P / MMHG_TO_CGS; 

    pout << "rcr_bcs_on = " << d_rcr_bcs_on << "\n"; 
    pout << "time\t P_lvot (mmHg)\t P_rvot (mmHg)\t P_aorta (mmHg)\t P_rpa (mmHg)\t P_lpa (mmHg)\t Q_lvot (ml/s)\t Q_rvot (ml/s)\t Q_aorta (ml/s)\t Q_rpa (ml/s)\t Q_lpa (ml/s)\t Q_valve (ml/s)\t idx_lv\t idx_rv\n";  
    pout << "\n";
    pout << d_time << " " << P_lvot << " " << P_rvot <<  " " << P_aorta << " " << P_rpa << " " << P_lpa << " " << d_Q_lvot << " " << d_Q_rvot << " " << d_Q_aorta << " " << d_Q_rpa << " " << d_Q_lpa << " " << d_Q_valve << " " << d_current_idx_series << " ";  
    pout << "\n";

}

int CirculationModel_preop::point_in_lvot(double testx, double testy, int axis, int side){

    if ((axis != d_lvot_axis) || (side != d_lvot_side))
        return 0;

    return pnpoly(d_n_pts_lvot, d_lvot_points_idx1, d_lvot_points_idx2, testx, testy);
}

int CirculationModel_preop::point_in_rvot(double testx, double testy, int axis, int side){
    // checks whether given point is in right ventricle

    // quick exit for correct side and axis 
    if ((axis != d_rvot_axis) || (side != d_rvot_side))
        return 0; 

    return pnpoly(d_n_pts_rvot, d_rvot_points_idx1, d_rvot_points_idx2, testx, testy); 
}

int CirculationModel_preop::point_in_aorta(double testx, double testy, int axis, int side){

    if ((axis != d_aorta_axis) || (side != d_aorta_side))
        return 0;

    return pnpoly(d_n_pts_aorta, d_aorta_points_idx1, d_aorta_points_idx2, testx, testy);
}

int CirculationModel_preop::point_in_rpa(double testx, double testy, int axis, int side){

    if ((axis != d_rpa_axis) || (side != d_rpa_side))
        return 0;

    return pnpoly(d_n_pts_rpa, d_rpa_points_idx1, d_rpa_points_idx2, testx, testy);
}

int CirculationModel_preop::point_in_lpa(double testx, double testy, int axis, int side){

    if ((axis != d_lpa_axis) || (side != d_lpa_side))
        return 0;

    return pnpoly(d_n_pts_lpa, d_lpa_points_idx1, d_lpa_points_idx2, testx, testy);
}


void CirculationModel_preop::write_plot_code()
{
    static const int mpi_root = 0;
    if (SAMRAI_MPI::getRank() == mpi_root)
    {
        ofstream fout(DATA_FILE_NAME.c_str(), ios::app);
        fout.setf(ios_base::scientific);
        fout.setf(ios_base::showpos);
        fout.precision(10);

        fout << "];\n";
        fout << "MMHG_TO_CGS = 1333.22368;\n";
        fout << "fig = figure;\n";
        fout << "times   =  bc_vals(:,1);\n";
        fout << "p_lvot  =  bc_vals(:,2);\n";
        fout << "p_rvot  =  bc_vals(:,3);\n";
        fout << "p_aorta =  bc_vals(:,4);\n";
        fout << "p_rpa   =  bc_vals(:,5);\n";
        fout << "p_lpa   =  bc_vals(:,6);\n";
        fout << "q_lvot  =  -bc_vals(:,7);\n";
        fout << "q_rvot  =  -bc_vals(:,8); \n";
        fout << "q_aorta =  bc_vals(:,9);\n";
        fout << "q_rpa   =  bc_vals(:,10);\n";
        fout << "q_lpa   =  bc_vals(:,11);\n";
        fout << "q_valve =  bc_vals(:,12);\n";
        fout << "subplot(2,1,1)\n";
        fout << "plot(times, p_lvot, ':k')\n";
        fout << "plot(times, p_rvot, '--k')\n";
        fout << "hold on\n";
        fout << "plot(times, p_aorta, 'b')\n";
        fout << "plot(times, p_rpa, ':b')\n";
        fout << "plot(times, p_lpa, '--b')\n";
        fout << "legend('P_{lvot}', 'P_{rvot}', 'P_{aorta}', 'P_{rpa}', 'P_{lpa}', 'Location','NorthEastOutside');\n";
        fout << "xlabel('t (s)')\n";
        fout << "ylabel('P (mmHg)')\n";
        fout << "subplot(2,1,2)\n";
        fout << "plot(times, q_lvot, ':k')\n";
        fout << "plot(times, q_rvot, '--k')\n";
        fout << "hold on\n";
        fout << "plot(times, q_aorta, 'b')\n";
        fout << "plot(times, q_rpa, ':b')\n";
        fout << "plot(times, q_lpa, '--b')\n";
        fout << "legend('Q_{lvot}', 'Q_{rvot}', 'Q_{aorta}', 'Q_{rpa}', 'Q_{lpa}', 'Location', 'NorthEastOutside')\n";
        fout << "xlabel('t (s)')\n";
        fout << "ylabel('Flow (ml/s)')\n";
        fout << "set(fig, 'Position', [100, 100, 1000, 750])\n";
        fout << "set(fig,'PaperPositionMode','auto')\n";
        fout << "saveas(fig, 'preop_pressure_flow')\n";
        fout << "save bc_data.mat\n";
    }
    return;
}



/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
    CirculationModel_preop::writeDataFile() const
{
    static const int mpi_root = 0;
    if (SAMRAI_MPI::getRank() == mpi_root)
    {
        static bool file_initialized = false;
        const bool from_restart = RestartManager::getManager()->isFromRestart();
        if (!from_restart && !file_initialized)
        {
            ofstream fout(DATA_FILE_NAME.c_str(), ios::out);
            fout << "% time\t P_lvot (mmHg)\t P_rvot (mmHg)\t P_aorta (mmHg)\t P_rpa (mmHg)\t P_lpa (mmHg)\t Q_lvot (ml/s)\t Q_rvot (ml/s)\t Q_aorta (ml/s)\t Q_rpa (ml/s)\t Q_lpa (ml/s)\t Q_valve (ml/s)";  
            fout << "\n"
                 << "bc_vals = [";
            file_initialized = true;
        }

        ofstream fout(DATA_FILE_NAME.c_str(), ios::app);

        fout << d_time;
        fout.setf(ios_base::scientific);
        fout.setf(ios_base::showpos);
        fout.precision(10);

        double P_lvot = d_lvot_P / MMHG_TO_CGS;
        double P_rvot = d_rvot_P / MMHG_TO_CGS;

        double P_aorta = 0.0;
        double P_rpa = 0.0;
        double P_lpa = 0.0; 

        if (d_rcr_bcs_on){
            P_aorta        = d_aorta_P/MMHG_TO_CGS;
            P_rpa          = d_rpa_P/MMHG_TO_CGS;
            P_lpa          = d_lpa_P/MMHG_TO_CGS;
        }

        fout << " " << P_lvot << " " << P_rvot << " " << P_aorta << " " << P_rpa << " " << P_lpa;
        fout << " " << d_Q_lvot << " " << d_Q_rvot << " " << d_Q_aorta << " " << d_Q_rpa << " " << d_Q_lpa << " " << d_Q_valve;

        if (d_lvot_0D_on){
            fout << " " << d_lvot_0D->d_V_lvot;
            fout << " " << d_lvot_0D->d_V_rest_lvot;
            fout << " " << d_lvot_0D->d_Elas;
            fout << " " << d_lvot_0D->d_act_temp;
            fout << " " << d_lvot_0D->d_Q_in;
            fout << " " << d_lvot_0D->d_P_lvot_in/MMHG_TO_CGS;
            fout << " " << d_lvot_0D->d_P_lvot_upstream/MMHG_TO_CGS;
        }
        if (d_rvot_0D_on){
            fout << " " << d_rvot_0D->d_V_rvot;
            fout << " " << d_rvot_0D->d_V_rest_rvot;
            fout << " " << d_rvot_0D->d_Elas;
            fout << " " << d_rvot_0D->d_act_temp;
            fout << " " << d_rvot_0D->d_Q_in;
            fout << " " << d_rvot_0D->d_P_rvot_in/MMHG_TO_CGS;
            fout << " " << d_rvot_0D->d_P_rvot_upstream/MMHG_TO_CGS;
        }

        fout << "; \n";

    }

    return;
} // writeDataFile

void
CirculationModel_preop::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("Restart database corresponding to " << d_object_name << " not found in restart file.");
    }

    d_current_idx_series         = db->getInteger("d_current_idx_series");
    d_Q_lvot          = db->getDouble("d_Q_lvot");
    d_Q_rvot          = db->getDouble("d_Q_rvot");
    d_Q_aorta         = db->getDouble("d_Q_aorta"); 
    d_Q_rpa           = db->getDouble("d_Q_rpa");
    d_Q_lpa           = db->getDouble("d_Q_lpa");
    d_Q_valve                    = db->getDouble("d_Q_valve");
    d_lvot_P          = db->getDouble("d_lvot_P");
    d_rvot_P          = db->getDouble("d_rvot_P");
    d_aorta_P             = db->getDouble("d_aorta_P");
    d_aorta_P_Wk          = db->getDouble("d_aorta_P_Wk");
    d_rpa_P                 = db->getDouble("d_rpa_P");
    d_rpa_P_Wk              = db->getDouble("d_rpa_P_Wk");
    d_lpa_P                  = db->getDouble("d_lpa_P");
    d_lpa_P_Wk               = db->getDouble("d_lpa_P_Wk");
    d_time                       = db->getDouble("d_time");
    d_rcr_bcs_on                 = db->getBool("d_rcr_bcs_on");
    d_lvot_0D_on     = db->getBool("d_lvot_0D_on");
    d_rvot_0D_on          = db->getBool("d_rvot_0D_on");
    d_P_min_linear_interp = db->getDouble("d_P_min_linear_interp"); 
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
