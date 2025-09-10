// Filename: CirculationModel_preop.h
// Created on 04 May 2007 by Boyce Griffith

// Modified 2019, Alexander D. Kaiser


#ifndef included_CirculationModel_preop
#define included_CirculationModel_preop

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/ibtk_utilities.h>

// SAMRAI INCLUDES
#include <CellVariable.h>
#include <tbox/Database.h>
#include <tbox/Serializable.h>

// C++ STDLIB INCLUDES
#include <vector>

// NAMESPACE
#include <ibamr/app_namespaces.h>

#include <boundary_condition_util.h>

#include "pnpoly.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class CirculationModel_preop
 */
class CirculationModel_preop : public Serializable
{
public:
    /*!
     * \brief The object name.
     */
    string d_object_name;

    /*!
     * \brief Whether the object is registered with the restart manager.
     */
    bool d_registered_for_restart;

    /*!
     * \brief model data.
     */
    const fourier_series_data *d_fourier_lvot;
    int     d_n_pts_lvot;
    double* d_lvot_points_idx1;
    double* d_lvot_points_idx2;
    int     d_lvot_axis; 
    int     d_lvot_side; 
    double  d_lvot_P;

    const fourier_series_data *d_fourier_rvot;
    int     d_n_pts_rvot;
    double* d_rvot_points_idx1;
    double* d_rvot_points_idx2;
    int     d_rvot_axis;     
    int     d_rvot_side;     
    double  d_rvot_P;

    lvot_0D_model *d_lvot_0D;
    rvot_0D_model *d_rvot_0D;     

    bool d_rcr_bcs_on; 
    bool d_lvot_0D_on; 
    bool d_rvot_0D_on;

    int     d_n_pts_aorta;
    double* d_aorta_points_idx1;
    double* d_aorta_points_idx2;
    int     d_aorta_axis; 
    int     d_aorta_side;
    double  d_aorta_P;
    double  d_aorta_P_Wk;
    double  d_aorta_R_proximal; 
    double  d_aorta_R_distal; 
    double  d_aorta_C;

    int     d_n_pts_rpa;
    double* d_rpa_points_idx1;
    double* d_rpa_points_idx2;
    int     d_rpa_axis;
    int     d_rpa_side;
    double  d_rpa_P;
    double  d_rpa_P_Wk;
    double  d_rpa_R_proximal;
    double  d_rpa_R_distal;
    double  d_rpa_C;

    int     d_n_pts_lpa;
    double* d_lpa_points_idx1;
    double* d_lpa_points_idx2;
    int     d_lpa_axis;
    int     d_lpa_side;
    double  d_lpa_P;
    double  d_lpa_P_Wk;
    double  d_lpa_R_proximal;
    double  d_lpa_R_distal;
    double  d_lpa_C;

    double  d_cycle_duration;
    double  d_t_offset_bcs_unscaled;
    double  d_t_offset_bcs_unscaled_lvot;
    double  d_t_offset_bcs_unscaled_rvot;
    unsigned int  d_current_idx_series; 
    double        d_Q_lvot;
    double        d_Q_rvot; 
    double        d_Q_aorta;
    double        d_Q_rpa;
    double        d_Q_lpa;
    double        d_Q_valve;
    double        d_time;
    double        d_time_lvot;
    double        d_time_rvot; 
    double        d_area_lvot;
    double        d_area_rvot; 
    double        d_area_aorta;
    double        d_area_rpa; 
    double        d_area_lpa; 
    bool          d_area_initialized;

    double d_p_extender_mean; 
    double d_p_extender_point;  
    double d_p_equal_fraction; 
    double d_P_min_linear_interp; 
    double d_rcr_on_time; 

    /*!
     * \brief The level of the patch hierarchy on which the Lagrangian
     * structures that interface the boundary are located.
     */
    int d_bdry_interface_level_number;

    /*!
     * \brief Constructor
     */
    CirculationModel_preop(Pointer<Database> input_db, 
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
                                               double rcr_on_time);
    /*!
     * \brief Destructor.
     */
    virtual ~CirculationModel_preop();

    /*!
     * \brief Advance time-dependent data.
     */
    void advanceTimeDependentData(const double dt,
                                  const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                  const int U_idx,
                                  const int P_idx,
                                  const int wgt_cc_idx,
                                  const int wgt_sc_idx);

    void set_Q_valve(double Q_valve); 

    void set_extender_pressures(double p_extender_mean, double p_extender_point); 

    /*!
     * \name Implementation of Serializable interface.
     */

    /*!
     * Write out object state to the given database.
     *
     * When assertion checking is active, database point must be non-null.
     */
    void putToDatabase(Pointer<Database> db);

    // basic data summary to stdout 
    void print_summary(); 
    void print_bc_debug();

    int point_in_lvot(double testx, double testy, int axis, int side);
    int point_in_rvot(double testx, double testy, int axis, int side);
    int point_in_aorta(double testx, double testy, int axis, int side);
    int point_in_rpa(double testx, double testy, int axis, int side);
    int point_in_lpa(double testx, double testy, int axis, int side);

    void write_plot_code(); 

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CirculationModel_preop(const CirculationModel_preop& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CirculationModel_preop& operator=(const CirculationModel_preop& that);

    /*!
     * Write out source/sink state data to disk.
     */
    void writeDataFile() const;

    /*!
     * Read object state from the restart file and initialize class data
     * members.  The database from which the restart data is read is determined
     * by the object_name specified in the constructor.
     *
     * Unrecoverable Errors:
     *
     *    -   The database corresponding to object_name is not found in the
     *        restart file.
     *
     *    -   The class version number and restart version number do not match.
     *
     */
    void getFromRestart();
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <CirculationModel.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CirculationModel_preop
