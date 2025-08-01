// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// Modified 2019, Alexander D. Kaiser

// Config files
// #include <IBAMR_config.h>
// #include <IBTK_config.h>
// #include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBSpringForceSpec.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/StaggeredStokesOpenBoundaryStabilizer.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>



// added this header
#include <ibamr/IBTargetPointForceSpec.h>
#include <ibamr/IBInstrumentPanel.h>
#include <ibamr/IBStandardSourceGen.h>
#include <vector>
#include <queue>
#include <cmath>
#include <timing.h>
#include <boundary_condition_util.h>
#include <CirculationModel.h>
#include <CirculationModel_with_lv.h>
#include <CirculationModel_aorta.h>
#include <FeedbackForcer.h>
#include <FourierBodyForce.h>


// #define IMPLICIT_SOLVER
#ifdef IMPLICIT_SOLVER
     #include <ibamr/IBImplicitStaggeredHierarchyIntegrator.h>
#endif


// #if defined(IBAMR_HAVE_SILO)
// #include <silo.h>
// #endif


inline double spring_function_collagen(double R, const double* params, int lag_mastr_idx, int lag_slave_idx);
inline double deriv_spring_collagen(double R, const double* params, int lag_mastr_idx, int lag_slave_idx);

inline double spring_function_aortic_circ(double R, const double* params, int lag_mastr_idx, int lag_slave_idx);
inline double deriv_spring_aortic_circ(double R, const double* params, int lag_mastr_idx, int lag_slave_idx);
inline double spring_function_aortic_rad(double R, const double* params, int lag_mastr_idx, int lag_slave_idx);
inline double deriv_spring_aortic_rad(double R, const double* params, int lag_mastr_idx, int lag_slave_idx);

inline double spring_function_compressive_only_linear_spring(double R, const double* params, int lag_mastr_idx, int lag_slave_idx);
inline double deriv_spring_compressive_only_linear_spring(double R, const double* params, int lag_mastr_idx, int lag_slave_idx);


#define DEBUG_OUTPUT 0 
#define ENABLE_INSTRUMENTS
// #define EXTENDER_INST_ON

#define USE_CIRC_MODEL_AORTA

#define MMHG_TO_CGS 1333.22368
#define CGS_TO_MMHG 0.000750061683
#define MPa_TO_CGS 1.0e7

namespace{
    inline double smooth_kernel(const double r){
        return std::abs(r) < 1.0 ? 0.5 * (cos(M_PI * r) + 1.0) : 0.0;
    } // smooth_kernel
}



/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int main(int argc, char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, NULL, NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    #ifdef USE_CIRC_MODEL_AORTA
        pout << "USE_CIRC_MODEL_AORTA is defined\n"; 
    #endif

    { // cleanup dynamically allocated objects prior to shutdown
    
        // time step used in various places throughout 
        double dt; 
        
        // make some timers
        timestamp_type time1_total, time2_total;    // For total time
        timestamp_type time1, time2;                // For step time
        get_timestamp(&time1_total);                // Get initialization time too

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
                
        // check if have a restarted run 
        string restart_read_dirname;
        // int restore_num = 0;
        bool from_restart = false;
        if (argc >= 4)
        {
            // Check whether this appears to be a restarted run.
            FILE* fstream = (SAMRAI_MPI::getRank() == 0 ? fopen(argv[2], "r") : NULL);
            if (SAMRAI_MPI::bcast(fstream != NULL ? 1 : 0, 0) == 1)
            {
                restart_read_dirname = argv[2];
                // restore_num = atoi(argv[3]);
                from_restart = true;
            }
            if (fstream != NULL)
            {
                fclose(fstream);
            }
        }

    
        #ifdef IMPLICIT_SOLVER
            // Read default Petsc options
            if (input_db->keyExists("petsc_options_file"))
            {
                std::string PetscOptionsFile = input_db->getString("petsc_options_file");
                #if (!PETSC_VERSION_RELEASE)
                    PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, PetscOptionsFile.c_str(), PETSC_TRUE);
                #else
                    PetscOptionsInsertFile(PETSC_COMM_WORLD, PetscOptionsFile.c_str(), PETSC_TRUE);
                #endif
            }
        #endif

        int n_restarts_written = 0;
        int max_restart_to_write = 20;
        if (input_db->keyExists("MAX_RESTART_TO_WRITE")){
            max_restart_to_write = input_db->getInteger("MAX_RESTART_TO_WRITE");
        }

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type =
            app_initializer->getComponentDatabase("Main")->getStringWithDefault("solver_type", "STAGGERED");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                                                   << "Valid options are: COLLOCATED, STAGGERED");
        }
        Pointer<IBMethod> ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
        
        
        #ifdef IMPLICIT_SOLVER
            Pointer<IBHierarchyIntegrator> time_integrator =
                new IBImplicitStaggeredHierarchyIntegrator("IBHierarchyIntegrator",
                                                            app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                                            ib_method_ops,
                                                            navier_stokes_integrator);
        #else
            Pointer<IBHierarchyIntegrator> time_integrator =
                new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                                    app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                                    ib_method_ops,
                                                    navier_stokes_integrator);
        #endif
        
        
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Configure the IB solver.
        Pointer<IBStandardInitializer> ib_initializer = new IBStandardInitializer(
            "IBStandardInitializer", app_initializer->getComponentDatabase("IBStandardInitializer"));
        ib_method_ops->registerLInitStrategy(ib_initializer);
        Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
        
        // adding custom function for collagen springs
        ib_force_fcn->registerSpringForceFunction(1, &spring_function_collagen,    &deriv_spring_collagen);
        ib_force_fcn->registerSpringForceFunction(2, &spring_function_aortic_circ, &deriv_spring_aortic_circ);
        ib_force_fcn->registerSpringForceFunction(3, &spring_function_aortic_rad,  &deriv_spring_aortic_rad);
        ib_force_fcn->registerSpringForceFunction(4, &spring_function_compressive_only_linear_spring,  &deriv_spring_compressive_only_linear_spring);

        ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);
        
        
        #ifdef ENABLE_INSTRUMENTS
         
            std::vector<double> flux_valve_ring;
            std::vector<double> p_extender_mean;
            std::vector<double> p_extender_point;

            // set below in loop
            int u_data_idx = 0; 
            int p_data_idx = 0; 
            
            Pointer<IBInstrumentPanel> instruments          = new IBInstrumentPanel("meter_0", input_db);
            #ifdef EXTENDER_INST_ON
                Pointer<IBInstrumentPanel> instruments_extender = new IBInstrumentPanel("aorta_0", input_db);
            #endif 
        #endif
        

        LDataManager* l_data_manager = ib_method_ops->getLDataManager();
        pout << "passed LDataManager creation\n" ; 

        pout << "to boundary and initial condition creation\n" ; 

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }
                
        // This is needed to pull variables later
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

        #ifndef USE_CIRC_MODEL_AORTA
            // Create Eulerian boundary condition specification objects (when necessary).
            vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, static_cast<RobinBcCoefStrategy<NDIM>*>(NULL));        
            const bool periodic_domain = grid_geometry->getPeriodicShift().min() > 0;
            if (!periodic_domain)
            {   
                pout << "Using b.c. from file, no series for boundary conditions (body force may still have series).\n";
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    ostringstream bc_coefs_name_stream;
                    bc_coefs_name_stream << "u_bc_coefs_" << d;
                    const string bc_coefs_name = bc_coefs_name_stream.str();
                    ostringstream bc_coefs_db_name_stream;
                    bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                    const string bc_coefs_db_name = bc_coefs_db_name_stream.str();
                    u_bc_coefs[d] = new muParserRobinBcCoefs(
                        bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
                }
                navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
             
                if (solver_type == "STAGGERED" && input_db->keyExists("BoundaryStabilization"))
                {
                    time_integrator->registerBodyForceFunction(new StaggeredStokesOpenBoundaryStabilizer(
                        "BoundaryStabilization",
                        app_initializer->getComponentDatabase("BoundaryStabilization"),
                        navier_stokes_integrator,
                        grid_geometry));
                }
            }
        #endif

        // generic body force
        // Create Eulerian body force function specification objects.
        if (input_db->keyExists("ForcingFunction"))
        {
            if (input_db->keyExists("BoundaryStabilization"))
            {
                TBOX_ERROR("Cannot currently use boundary stabilization with additional body forcing");
            }
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        #ifdef FOURIER_SERIES_BODY_FORCE
            #ifdef USE_CIRC_MODEL_AORTA
                TBOX_ERROR("Cannot have both FOURIER_SERIES_BODY_FORCE and USE_CIRC_MODEL_AORTA set")
            #endif  
        #endif 


        #ifdef USE_CIRC_MODEL_AORTA

            bool damping_outside = true; 
             
            dt = input_db->getDouble("DT");

            // End systolic / beginning diastolic PA pressure
            double P_aorta_0;
            if (input_db->keyExists("P_aorta_0_MMHG")){
                P_aorta_0 = input_db->getDouble("P_aorta_0_MMHG") * MMHG_TO_CGS;
            }
            else{
                P_aorta_0 = 100.0 * MMHG_TO_CGS;
            }

            std::string fourier_coeffs_name_rv; 
            if (input_db->keyExists("FOURIER_COEFFS_FILENAME_VENTRICLE")){
                fourier_coeffs_name_rv = input_db->getString("FOURIER_COEFFS_FILENAME_VENTRICLE");
            }
            else {
                fourier_coeffs_name_rv = "fourier_coeffs_ventricle.txt"; 
            }
            pout << "attempting to build series with from file: " << fourier_coeffs_name_rv << "\n"; 
            fourier_series_data *fourier_series_ventricle = new fourier_series_data(fourier_coeffs_name_rv.c_str(), dt);

            // boundary vertex files 
            std::string ventricle_vertices_file_name; 
            if (input_db->keyExists("BOUNDARY_FILENAME_VENTRICLE")){
                ventricle_vertices_file_name = input_db->getString("BOUNDARY_FILENAME_VENTRICLE");
            }
            else {
                ventricle_vertices_file_name = "ventricle_bdry.vertex"; 
            }

            std::string aorta_vertices_file_name; 
            if (input_db->keyExists("BOUNDARY_FILENAME_AORTA")){
                aorta_vertices_file_name = input_db->getString("BOUNDARY_FILENAME_AORTA");
            }
            else {
                aorta_vertices_file_name = "aorta_bdry.vertex"; 
            }

            std::string vessel_file_name; 
            if (input_db->keyExists("VESSEL_FILENAME")){
                vessel_file_name = input_db->getString("VESSEL_FILENAME");
            }
            else {
                pout << "VESSEL_FILENAME not found, damping_outside vessel is off\n"; 
                damping_outside = false; 
                vessel_file_name = ""; 
            }

            // scaled cycle length for this patient 
            double t_cycle_length = input_db->getDouble("CYCLE_DURATION");

            // start this far into the Fourier series 
            double t_offset_start_bcs_unscaled; 
            if (input_db->keyExists("T_OFFSET_START_BCS_UNSCALED"))
            {
                t_offset_start_bcs_unscaled = input_db->getDouble("T_OFFSET_START_BCS_UNSCALED"); // starts this 
            }
            else {
                t_offset_start_bcs_unscaled = 0.0; 
            }

            double rcr_on_time; 
            if (input_db->keyExists("RCR_ON_TIME")){
                rcr_on_time = input_db->getDouble("RCR_ON_TIME"); // starts this 
            }
            else {
                rcr_on_time = 0.2; 
            }

            // start in physical time with relation to Fourier series 
            double t_offeset_start = t_offset_start_bcs_unscaled * (t_cycle_length / fourier_series_ventricle->L);

            bool rcr_bcs_on = true;

            bool ventricle_0D_on = input_db->getBoolWithDefault("VENTRICLE_0D_ON", false);

            // start with a linear ramp up in pressure 
            bool P_initial_aorta_equal_to_ventricle = true; 
            // double rcr_on_time = 0.2; 

            CirculationModel_aorta *circ_model_aorta = new CirculationModel_aorta(input_db,
                                                                             fourier_series_ventricle, 
                                                                             ventricle_vertices_file_name,
                                                                             aorta_vertices_file_name,
                                                                             t_cycle_length,
                                                                             t_offset_start_bcs_unscaled, 
                                                                             time_integrator->getIntegratorTime(), 
                                                                             P_aorta_0,
                                                                             rcr_bcs_on, 
                                                                             ventricle_0D_on,
                                                                             P_initial_aorta_equal_to_ventricle, 
                                                                             rcr_on_time); 

            // Create Eulerian boundary condition specification objects.
            vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
            for (int d = 0; d < NDIM; ++d){
                u_bc_coefs[d] = new VelocityBcCoefs_aorta(d, circ_model_aorta);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);

            // flow straightener at boundary 
            Pointer<FeedbackForcer> feedback_forcer = new FeedbackForcer(navier_stokes_integrator, patch_hierarchy, NULL, NULL, circ_model_aorta, damping_outside, vessel_file_name, aorta_vertices_file_name);
            time_integrator->registerBodyForceFunction(feedback_forcer);


        #endif // #ifdef USE_CIRC_MODEL_AORTA


        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();
        if (uses_visit)
        {
            ib_initializer->registerLSiloDataWriter(silo_data_writer);
            ib_method_ops->registerLSiloDataWriter(silo_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        pout << "before initializePatchHierarchy\n" ; 
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        pout << "passed initializePatchHierarchy\n" ; 

        // Finest level does not change throughout 
        const int finest_hier_level = patch_hierarchy->getFinestLevelNumber();
        

        #ifdef ENABLE_INSTRUMENTS
            // do this after initialize patch hierarchy
            instruments->initializeHierarchyIndependentData(patch_hierarchy, l_data_manager);

            #ifdef EXTENDER_INST_ON
                instruments_extender->initializeHierarchyIndependentData(patch_hierarchy, l_data_manager);
            #endif
        #endif

        // Deallocate initialization objects.
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();
        app_initializer.setNull();

        // Setup Silo writer.
        if (silo_data_writer)
        {
            //const int finest_hier_level = patch_hierarchy->getFinestLevelNumber();
            Pointer<LData> F_data = l_data_manager->getLData("F", finest_hier_level);
            silo_data_writer->registerVariableData("F", F_data, finest_hier_level);
        }

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        
        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            silo_data_writer->writePlotData(iteration_num, loop_time);
            
            #ifdef ENABLE_INSTRUMENTS
                
                if (instruments->isInstrumented())
                    pout << "is instrumented passed\n" ;
                else
                    pout << "is instrumented returned FALSE\n" ;  
                    
                instruments->initializeHierarchyDependentData(patch_hierarchy, l_data_manager, iteration_num, loop_time); 
                instruments->readInstrumentData(u_data_idx, p_data_idx, patch_hierarchy, l_data_manager, iteration_num, loop_time); 
                
                #ifdef EXTENDER_INST_ON
                    instruments_extender->initializeHierarchyDependentData(patch_hierarchy, l_data_manager, iteration_num, loop_time); 
                    instruments_extender->readInstrumentData(u_data_idx, p_data_idx, patch_hierarchy, l_data_manager, iteration_num, loop_time); 
                #endif
            #endif
            
        }
 
        
        double dt_original = time_integrator->getMaximumTimeStepSize();


        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        dt = 0.0;
        double dt_prev = 0.0;
        bool prev_step_initialized = false;
        
        get_timestamp(&time2_total);
        double total_init_time = timestamp_diff_in_seconds(time1_total, time2_total);
        pout << "\n\nTotal initialization time = " << total_init_time << "\n\n"; 
        
        // add some timers         
        get_timestamp(&time1_total); 
        double step_time; 

        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            
            get_timestamp(&time1);          // start step clock 
            
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();


            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            
            // save the last time step 
            dt_prev = dt; 
            
            // get current step 
            dt = time_integrator->getMaximumTimeStepSize();
        
            if((!prev_step_initialized) && from_restart){
                    dt_prev = dt;
                    prev_step_initialized = true;
            }
                                    
            // step the whole thing
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            // stop step clock 
            get_timestamp(&time2); 
            step_time = timestamp_diff_in_seconds(time1, time2); 
            
            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "Wallclock time elapsed = " << step_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                silo_data_writer->writePlotData(iteration_num, loop_time);
                
            }
            
            // Write this every step
            #ifdef ENABLE_INSTRUMENTS
                if (instruments->isInstrumented()){
                    Pointer<hier::Variable<NDIM> > U_var = navier_stokes_integrator->getVelocityVariable();
                    Pointer<hier::Variable<NDIM> > P_var = navier_stokes_integrator->getPressureVariable();
                    Pointer<VariableContext> current_ctx = navier_stokes_integrator->getCurrentContext();
                    const int U_current_idx = var_db->mapVariableAndContextToIndex(U_var, current_ctx);
                    const int P_current_idx = var_db->mapVariableAndContextToIndex(P_var, current_ctx);
                    
                    instruments->initializeHierarchyDependentData(patch_hierarchy, l_data_manager, iteration_num, loop_time); 
                    instruments->readInstrumentData(U_current_idx, P_current_idx, patch_hierarchy, l_data_manager, iteration_num, loop_time); 
                    flux_valve_ring = instruments->getFlowValues(); 

                    #ifdef EXTENDER_INST_ON
                        instruments_extender->initializeHierarchyDependentData(patch_hierarchy, l_data_manager, iteration_num, loop_time); 
                        instruments_extender->readInstrumentData(U_current_idx, P_current_idx, patch_hierarchy, l_data_manager, iteration_num, loop_time); 
                        p_extender_mean = instruments_extender->getMeanPressureValues(); 
                        p_extender_point = instruments_extender->getPointwisePressureValues(); 
                    #endif

                    // pout << "flux at t = " << loop_time << ", Q = " << flux_valve_ring[0] << "\n";
                                    
                    #ifdef USE_CIRC_MODEL_AORTA
                        circ_model_aorta->set_Q_valve(flux_valve_ring[0]); 
                        #ifdef EXTENDER_INST_ON
                            circ_model_aorta->set_extender_pressures(p_extender_mean[0], p_extender_point[0]); 
                        #endif
                    #endif 
                }
            #endif
            
            // Update the circulation model if used 
            #ifdef USE_CIRC_MODEL_AORTA
                {
                    Pointer<hier::Variable<NDIM> > U_var = navier_stokes_integrator->getVelocityVariable();
                    Pointer<hier::Variable<NDIM> > P_var = navier_stokes_integrator->getPressureVariable();
                    Pointer<VariableContext> current_ctx = navier_stokes_integrator->getCurrentContext();
                    const int U_current_idx = var_db->mapVariableAndContextToIndex(U_var, current_ctx);
                    const int P_current_idx = var_db->mapVariableAndContextToIndex(P_var, current_ctx);
                    Pointer<HierarchyMathOps> hier_math_ops = navier_stokes_integrator->getHierarchyMathOps();
                    const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();
                    const int wgt_sc_idx = hier_math_ops->getSideWeightPatchDescriptorIndex();
                    circ_model_aorta->advanceTimeDependentData(dt, patch_hierarchy, U_current_idx, P_current_idx, wgt_cc_idx, wgt_sc_idx);
                }
            #endif
            
            
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
                
                n_restarts_written++;
                
                if (n_restarts_written > max_restart_to_write){
                    if (SAMRAI_MPI::getRank() == 0){
                        std::ofstream controlled_stop_stream;
                        controlled_stop_stream.open("controlled_stop.txt", ios_base::out | ios_base::trunc);
                        controlled_stop_stream << "stopped\n" ;
                        controlled_stop_stream.close();
                    }
                    SAMRAI_MPI::barrier();
                    SAMRAI_MPI::abort();
                }
            }


            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            
        }

        get_timestamp(&time2_total); 
        double total_time = timestamp_diff_in_seconds(time1_total, time2_total); 
        double average_time = total_time / ((double) iteration_num); 
        
        pout << "total run time = " << total_time << " s. \n" ; 
        pout << "average run time = " << average_time << " s. \n" ;
                
        #ifdef USE_CIRC_MODEL_AORTA
            circ_model_aorta->write_plot_code(); 
        #endif 
        
        if (SAMRAI_MPI::getRank() == 0){
            std::ofstream done_stream;
            done_stream.open("done.txt", ios_base::out | ios_base::trunc);
            done_stream << "done\n" ;
            done_stream.close(); 
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
} // main


inline double spring_function_collagen(double R, const double* params, int lag_mastr_idx, int lag_slave_idx){
    /* 
    Compute force for collagen springs
    
    Zero force under compression
    
    Exponential growth until full recruitment,
    experimentally determined strain at which microscopic fibers align
     
    Linear after full recruitment
    
    All collagen assumed to have same modulus
    General parameters params[0] are normalized,
    but if a stiffer spring is requested then it is used here.
    
    params[0] must include the width element ds for converting forces to areas 
    */
    
    static const double a                    = 4643.4;       // Coeff of exponential term
    static const double b                    = 49.9643;      // Exponential rate
    static const double full_recruitment     = 0.145;             // Linear at strains larger than this
    static const double eta_collagen         = 32.5 * MPa_TO_CGS; // Linear slope, in barye = dynes/cm^2 = g cm/(s cm^2)
    static const double collagen_x_intercept = 0.125;             // Linear collagen part intercepts x axis at this strain
    static const double collagen_y_intercept = -collagen_x_intercept * eta_collagen; // Linear collagen part intercepts y axis at this stress
    //static const double thickness            = 0.1; // cm, Thickness of leaflet tissue, for converting strains to forces
    
    // static const double ds, width element of the spring, also for converting strains to forces
    // included in params[0]
    
    const double kappa    = params[0];
    const double rest_len = params[1];
    
    // Strain, dimension 1
    const double E = R/rest_len - 1.0;
    
    // Compute the force
    if (E > full_recruitment){
        /*if ((lag_mastr_idx % 2500) == 0){
            std::cout << "Affine. (idx,nbr) = (" << lag_mastr_idx << ", " <<  lag_slave_idx
                      << "\tE = " << E
                      << "\tF = " << kappa * (eta_collagen*E + collagen_y_intercept)
                      << "\tEffective slope = " << kappa * eta_collagen
                      << "\tRest len = " << rest_len
                      << "\n";
        }*/ 
        return kappa * (eta_collagen*E + collagen_y_intercept);
    }
    else if (E > 0.0){
        /*if ((lag_mastr_idx % 2500) == 0){
            std::cout << "Exp.   (idx,nbr) = (" << lag_mastr_idx << ", " <<  lag_slave_idx << ")"
                      << "\tE = " << E
                      << "\tF = " << kappa * a * (exp(b*E) - 1)
                      << "\tEffective slope = " << kappa * a * b // taylor series coefficient on first term
                      << "\tRest len = " << rest_len
                      << "\n";
        }*/ 
        return kappa * a * (exp(b*E) - 1);
    }
    else{
        return 0.0;
    }
    return 0.0;
} // spring_function_collagen


inline double deriv_spring_collagen(double R, const double* params, int lag_mastr_idx, int lag_slave_idx){
    // not implemented

    SAMRAI_MPI::abort();
    return 0.0;
} // deriv_spring_collagen



inline double spring_function_aortic_circ(double R, const double* params, int lag_mastr_idx, int lag_slave_idx){
    // function idx 2
    static const double b = 57.456509400487398; // Exponential rate

    const double kappa    = params[0];
    const double rest_len = params[1];
    
    // Strain, dimension 1
    const double E = R/rest_len - 1.0;
    
    // Compute the force
    if (E > 0.0){
        return kappa * (exp(b*E) - 1);
    }
    else{
        // linear for compressive strains with continuous slope at origin 
        return 0.0; //kappa * b *  E;
    }
    return 0.0;
} // spring_function_aortic_circ


inline double deriv_spring_aortic_circ(double R, const double* params, int lag_mastr_idx, int lag_slave_idx){
    // not implemented

    SAMRAI_MPI::abort();
    return 0.0;
} // deriv_spring_aortic_circ



inline double spring_function_aortic_rad(double R, const double* params, int lag_mastr_idx, int lag_slave_idx){
    // function idx 3
    static const double b = 22.397200094241359; // Exponential rate

    const double kappa    = params[0];
    const double rest_len = params[1];
    
    // Strain, dimension 1
    const double E = R/rest_len - 1.0;
    
    // Compute the force
    if (E > 0.0){
        return kappa * (exp(b*E) - 1);
    }
    else{
        // linear for compressive strains with continuous slope at origin 
        return 0.0; // kappa * b *  E;
    }
    return 0.0;
} // spring_function_aortic_rad


inline double deriv_spring_aortic_rad(double R, const double* params, int lag_mastr_idx, int lag_slave_idx){
    // not implemented

    SAMRAI_MPI::abort();
    return 0.0;
} // deriv_spring_aortic_rad


inline double spring_function_compressive_only_linear_spring(double R, const double* params, int lag_mastr_idx, int lag_slave_idx){
    // function idx 4
    const double kappa    = params[0];
    const double rest_len = params[1];
    
    // Strain, dimension 1
    const double E = R/rest_len - 1.0;
    
    // Compute the force
    if (E > 0.0){
        // zero under extension 
        return 0.0; 
    }
    else{
        // linear for compressive strains
        return kappa * E; 
    }
    return 0.0;
} // spring_function_compressive_only_linear_spring


inline double deriv_spring_compressive_only_linear_spring(double R, const double* params, int lag_mastr_idx, int lag_slave_idx){
    // not implemented

    SAMRAI_MPI::abort();
    return 0.0;
} // deriv_spring_compressive_only_linear_spring
