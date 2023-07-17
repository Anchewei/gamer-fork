#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double ParM;
static double BackDen;
static double RanVelL;
static double RanVelH;

static RandomNumber_t *RNG = NULL;

// =======================================================================================


//-------------------------------------------------------------------------------------------------------
// Function    :  RanVec_FixRadius
// Description :  Compute a random 3D vector with a fixed radius
//
// Note        :  Uniformly random sample in theta and phi does NOT give a uniformly random sample in 3D space
//                --> Uniformly random sample in a 3D sphere and then normalize all vectors to the given radius
//
// Parameter   :  r      : Input radius
//                RanVec : Array to store the random 3D vector
//
// Return      :  RanVec
//-------------------------------------------------------------------------------------------------------
void RanVec_FixRadius( const double r, double RanVec[] )
{

   double Norm, RanR2;

   do
   {
      RanR2 = 0.0;

      for (int d=0; d<3; d++)
      {
         RanVec[d]  = RNG->GetValue( 0, -1.0, +1.0 );
         RanR2     += SQR( RanVec[d] );
      }
   } while ( RanR2 > 1.0 );

   Norm = r / sqrt( RanR2 );

   for (int d=0; d<3; d++)    RanVec[d] *= Norm;

} // FUNCTION : RanVec_FixRadius



//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );

#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

// examples
/*
// errors
#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

#  ifdef GRAVITY
   if ( OPT__BC_FLU[0] == BC_FLU_PERIODIC  ||  OPT__BC_POT == BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "do not use periodic BC for this test !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
#     ifndef DUAL_ENERGY
         Aux_Message( stderr, "WARNING : it's recommended to enable DUAL_ENERGY for this test !!\n" );
#     endif

      if ( FLAG_BUFFER_SIZE < 5 )
         Aux_Message( stderr, "WARNING : it's recommended to set FLAG_BUFFER_SIZE >= 5 for this test !!\n" );
   } // if ( MPI_Rank == 0 )
*/


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//                3. Must call EoS_Init() before calling any other EoS routine
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "ParM",           &ParM,               1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "BackDen",        &BackDen,          1e-27,           Eps_double,       NoMax_double      );
   ReadPara->Add( "RanVelL",        &RanVelL,            0.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "RanVelH",        &RanVelH,            1.0,           Eps_double,       NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   ParM *= Const_Msun/UNIT_M;
   RanVelL *= Const_km/UNIT_V;
   RanVelH *= Const_km/UNIT_V;

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = __FLT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID           = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  ParM                      = %13.7e\n",        ParM );
      Aux_Message( stdout, "  BackDen                   = %13.7e\n",     BackDen );
      Aux_Message( stdout, "  RanVelL                   = %13.7e\n",     RanVelL );
      Aux_Message( stdout, "  RanVelH                   = %13.7e\n",     RanVelH );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//                4. For MHD, do NOT add magnetic energy (i.e., 0.5*B^2) to fluid[ENGY] here
//                   --> It will be added automatically later
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

// HYDRO example
   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;

   Dens = BackDen;
   MomX = 0.0;
   MomY = 0.0;
   MomZ = 0.0;
   Pres = 2.0;
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );   // assuming EoS requires no passive scalars
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here

// set the output array
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : SetGridIC



#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  SetBFieldIC
// Description :  Set the problem-specific initial condition of magnetic field
//
// Note        :  1. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  magnetic : Array to store the output magnetic field
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  magnetic
//-------------------------------------------------------------------------------------------------------
void SetBFieldIC( real magnetic[], const double x, const double y, const double z, const double Time,
                  const int lv, double AuxArray[] )
{

   /*
// example
   magnetic[MAGX] = 1.0;
   magnetic[MAGY] = 2.0;
   magnetic[MAGZ] = 3.0;
   */

} // FUNCTION : SetBFieldIC
#endif // #ifdef MHD
#endif // #if ( MODEL == HYDRO )

#  ifdef PARTICLE
//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction
// Description :  User-specified function to initialize particle attributes
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Par_Init_ByFunction_Ptr"
//                   --> This function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Periodicity should be taken care of in this function
//                   --> No particles should lie outside the simulation box when the periodic BC is adopted
//                   --> However, if the non-periodic BC is adopted, particles are allowed to lie outside the box
//                       (more specifically, outside the "active" region defined by amr->Par->RemoveCell)
//                       in this function. They will later be removed automatically when calling Par_Aux_InitCheck()
//                       in Init_GAMER().
//                3. Particles set by this function are only temporarily stored in this MPI rank
//                   --> They will later be redistributed when calling Par_LB_Init_RedistributeByRectangular()
//                       and LB_Init_LoadBalance()
//                   --> Therefore, there is no constraint on which particles should be set by this function
//
// Parameter   :  NPar_ThisRank : Number of particles to be set by this MPI rank
//                NPar_AllRank  : Total Number of particles in all MPI ranks
//                ParMass       : Particle mass     array with the size of NPar_ThisRank
//                ParPosX/Y/Z   : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z   : Particle velocity array with the size of NPar_ThisRank
//                ParTime       : Particle time     array with the size of NPar_ThisRank
//                ParType       : Particle type     array with the size of NPar_ThisRank
//                AllAttribute  : Pointer array for all particle attributes
//                                --> Dimension = [PAR_NATT_TOTAL][NPar_ThisRank]
//                                --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                    to access the data

//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, ParType, AllAttribute
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_ParTest( const long NPar_ThisRank, const long NPar_AllRank,
                          real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                          real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                          real *ParType, real *AllAttribute[PAR_NATT_TOTAL] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

   real *ParData_AllRank[PAR_NATT_TOTAL];
   for (int v=0; v<PAR_NATT_TOTAL; v++)   ParData_AllRank[v] = NULL;

   const double BoxSize[3]   = { amr->BoxSize[0], amr->BoxSize[1], amr->BoxSize[2] };
   const real   ParDist      = BoxSize[0]/(2*9);                          // distance betwwen adjacent particles
   const double Radius       = BoxSize[0]/(2*10);

   double RanPos, RanV, RanVec[3];

// only the master rank will construct the initial condition
   if ( MPI_Rank == 0 )
   {
      ParData_AllRank[PAR_MASS] = new real [NPar_AllRank];
      ParData_AllRank[PAR_POSX] = new real [NPar_AllRank];
      ParData_AllRank[PAR_POSY] = new real [NPar_AllRank];
      ParData_AllRank[PAR_POSZ] = new real [NPar_AllRank];
      ParData_AllRank[PAR_VELX] = new real [NPar_AllRank];
      ParData_AllRank[PAR_VELY] = new real [NPar_AllRank];
      ParData_AllRank[PAR_VELZ] = new real [NPar_AllRank];

//    initialize the random number generator
      RNG = new RandomNumber_t( 1 );
      RNG->SetSeed( 0, 123 );

//    set particle attributes
      ParData_AllRank[PAR_POSX][0] = BoxSize[0]/4; // position of the first particle
      ParData_AllRank[PAR_POSY][0] = BoxSize[1]/4;
      ParData_AllRank[PAR_POSZ][0] = BoxSize[2]/4;
      for (long p=0; p<NPar_AllRank; p++)
      {
         ParData_AllRank[PAR_MASS][p] = ParM;

         RanPos    = RNG->GetValue( 0, 0.0, Radius );
         RanVec_FixRadius( RanPos, RanVec );

         for (int i=0; i<3; i++) ParData_AllRank[PAR_POSX+i][p] = RanVec[i] + BoxSize[i]/2;

         RanV    = RNG->GetValue( 0, RanVelL, RanVelH );
         RanVec_FixRadius( RanV, RanVec );
         for (int i=0; i<3; i++) ParData_AllRank[PAR_VELX+i][p] = RanVec[i];
      } // for (long p=0; p<NPar_AllRank; p++)
   } // if ( MPI_Rank == 0 )

// send particle attributes from the master rank to all ranks
   Par_ScatterParticleData( NPar_ThisRank, NPar_AllRank, _PAR_MASS|_PAR_POS|_PAR_VEL, ParData_AllRank, AllAttribute );

// synchronize all particles to the physical time on the base level,
// set generic particle type
   for (long p=0; p<NPar_ThisRank; p++) {
      ParTime[p] = Time[0];
      ParType[p] = PTYPE_GENERIC_MASSIVE;
   }

// free resource
   if ( MPI_Rank == 0 )
   {
      delete RNG;
      for (int v=0; v<PAR_NATT_TOTAL; v++)   delete [] ParData_AllRank[v];
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction
#  endif


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_ParTest
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_ParTest()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// procedure to enable a problem-specific function:
// 1. define a user-specified function (example functions are given below)
// 2. declare its function prototype on the top of this file
// 3. set the corresponding function pointer below to the new problem-specific function
// 4. enable the corresponding runtime option in "Input__Parameter"
//    --> for instance, enable OPT__OUTPUT_USER for Output_User_Ptr
   Init_Function_User_Ptr            = SetGridIC;
#  ifdef MHD
   Init_Function_BField_User_Ptr     = SetBFieldIC;
   Init_BField_ByVecPot_User_Ptr     = NULL; // option: OPT__INIT_BFIELD_BYVECPOT=2;  example: Model_Hydro/MHD_Init_BField_ByVecPot_Function.cpp
#  endif
// comment out Init_ByFile_User_Ptr to use the default
// Init_ByFile_User_Ptr              = NULL; // option: OPT__INIT=3;                  example: Init/Init_ByFile.cpp -> Init_ByFile_Default()
   Init_Field_User_Ptr               = NULL; // set NCOMP_PASSIVE_USER;               example: TestProblem/Hydro/Plummer/Init_TestProb_Hydro_Plummer.cpp --> AddNewField()
   Flag_Region_Ptr                   = NULL; // option: OPT__FLAG_REGION;             example: Refing/Flag_Region.cpp
   Flag_User_Ptr                     = NULL; // option: OPT__FLAG_USER;               example: Refine/Flag_User.cpp
   Mis_GetTimeStep_User_Ptr          = NULL; // option: OPT__DT_USER;                 example: Miscellaneous/Mis_GetTimeStep_User.cpp
   Mis_UserWorkBeforeNextLevel_Ptr   = NULL; //                                       example: Miscellaneous/Mis_UserWorkBeforeNextLevel.cpp
   Mis_UserWorkBeforeNextSubstep_Ptr = NULL; //                                       example: Miscellaneous/Mis_UserWorkBeforeNextSubstep.cpp
   BC_User_Ptr                       = NULL; // option: OPT__BC_FLU_*=4;              example: TestProblem/ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> BC()
#  ifdef MHD
   BC_BField_User_Ptr                = NULL; // option: OPT__BC_FLU_*=4;
#  endif
   Flu_ResetByUser_Func_Ptr          = NULL; // option: OPT__RESET_FLUID;             example: Fluid/Flu_ResetByUser.cpp
   Init_DerivedField_User_Ptr        = NULL; // option: OPT__OUTPUT_USER_FIELD;       example: Fluid/Flu_DerivedField_User.cpp
   Output_User_Ptr                   = NULL; // option: OPT__OUTPUT_USER;             example: TestProblem/Hydro/AcousticWave/Init_TestProb_Hydro_AcousticWave.cpp --> OutputError()
   Output_UserWorkBeforeOutput_Ptr   = NULL; // option: none;                         example: Output/Output_UserWorkBeforeOutput.cpp
   Aux_Record_User_Ptr               = NULL; // option: OPT__RECORD_USER;             example: Auxiliary/Aux_Record_User.cpp
   Init_User_Ptr                     = NULL; // option: none;                         example: none
   End_User_Ptr                      = NULL; // option: none;                         example: TestProblem/Hydro/ClusterMerger_vs_Flash/Init_TestProb_ClusterMerger_vs_Flash.cpp --> End_ClusterMerger()
#  ifdef GRAVITY
   Init_ExtAcc_Ptr                   = NULL; // option: OPT__EXT_ACC;                 example: SelfGravity/CPU_Gravity/CPU_ExtAcc_PointMass.cpp
   End_ExtAcc_Ptr                    = NULL;
   Init_ExtPot_Ptr                   = NULL; // option: OPT__EXT_POT;                 example: SelfGravity/CPU_Poisson/CPU_ExtPot_PointMass.cpp
   End_ExtPot_Ptr                    = NULL;
   Poi_AddExtraMassForGravity_Ptr    = NULL; // option: OPT__GRAVITY_EXTRA_MASS;      example: none
   Poi_UserWorkBeforePoisson_Ptr     = NULL; // option: none;                         example: SelfGravity/Poi_UserWorkBeforePoisson.cpp
#  endif
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr           = Par_Init_ByFunction_ParTest; // option: PAR_INIT=1;                   example: Particle/Par_Init_ByFunction.cpp
   Par_Init_Attribute_User_Ptr       = NULL; // set PAR_NATT_USER;                    example: TestProblem/Hydro/AGORA_IsolatedGalaxy/Init_TestProb_Hydro_AGORA_IsolatedGalaxy.cpp --> AddNewParticleAttribute()
#  endif
#  if ( EOS == EOS_USER )
   EoS_Init_Ptr                      = NULL; // option: EOS in the Makefile;          example: EoS/User_Template/CPU_EoS_User_Template.cpp
   EoS_End_Ptr                       = NULL;
#  endif
#  endif // #if ( MODEL == HYDRO )
   Src_Init_User_Ptr                 = NULL; // option: SRC_USER;                     example: SourceTerms/User_Template/CPU_Src_User_Template.cpp
#  ifdef FEEDBACK
   FB_Init_User_Ptr                  = NULL; // option: FB_USER;                      example: TestProblem/Hydro/Plummer/FB_Plummer.cpp
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_ParTest
