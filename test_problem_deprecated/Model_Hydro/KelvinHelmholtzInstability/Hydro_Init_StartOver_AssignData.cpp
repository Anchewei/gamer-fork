#include "GAMER.h"

#if ( MODEL == HYDRO )

static void Init_Function_User( real fluid[], const double x, const double y, const double z, const double Time );
void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time ) = Init_Function_User;




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Function_User
// Description :  Function to initialize the fluid field
//
// Note        :  1. Invoked by "Hydro_Init_StartOver_AssignData"
//                2. This function will be invoked by multiple OpenMP threads
//                   --> Must ensure everything here is thread-safe
//
// Parameter   :  fluid : Fluid field to be initialized
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void Init_Function_User( real fluid[], const double x, const double y, const double z, const double Time )
{

   const double Gamma2  = 1.0/GAMMA/(GAMMA-1.0);
   const double C1[3] = { 0.5*amr->BoxSize[0]+100.0,
                          0.5*amr->BoxSize[1]+200.0,
                          0.5*amr->BoxSize[2]+300.0 };
   const double C2[3] = { 20.0, 40.0, 10.0 };

   const double Cs      =   1.0;
   const double Height1 = 100.0;
   const double Height2 = 400.0;
   const double Width1  = 640.0;
   const double Width2  = 512.0;

// set active variables
   fluid[DENS] = 1.0 + Height1*exp(  -( SQR(x-C1[0])+ SQR(y-C1[1]) + SQR(z-C1[2]) ) / SQR(Width1)  );
   fluid[DENS] +=      Height2*exp(  -( SQR(x-C2[0])+ SQR(y-C2[1]) + SQR(z-C2[2]) ) / SQR(Width2)  );
   fluid[MOMX] = 1.0;
   fluid[MOMY] = 2.0;
   fluid[MOMZ] = 3.0;
   fluid[ENGY] = Cs*Cs*fluid[DENS]*Gamma2 + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];

// set passive scalars

} // FUNCTION : Init_Function_User



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Init_StartOver_AssignData
// Description :  Construct the initial condition in HYDRO
//
// Note        :  1. Work for the option "OPT__INIT == INIT_STARTOVER"
//                2. The initialization function should be specified in "Init_Function_Ptr", which is a function
//                   pointer pointing to either "Init_Function_User" or the test problem specified function
//                   (e.g., Hydro_TestProbSol_Riemann) set in "Init_TestProb"
//
// Parameter   :  lv : Targeted refinement level
//-------------------------------------------------------------------------------------------------------
void Hydro_Init_StartOver_AssignData( const int lv )
{

   const int    NSub     = INIT_SUBSAMPLING_NCELL;
   const double dh       = amr->dh[lv];
   const double dh_sub   = dh / NSub;
   const double _NSub3   = 1.0/(NSub*NSub*NSub);
   const real   Gamma_m1 = GAMMA - (real)1.0;
   const real  _Gamma_m1 = (real)1.0 / Gamma_m1;

   real   fluid[NCOMP_TOTAL], fluid_sub[NCOMP_TOTAL];
   double x, y, z, x0, y0, z0;


   if ( NSub > 1 )   // with sub-sampling
   {
//#     pragma omp parallel for private( fluid, fluid_sub, x, y, z, x0, y0, z0 ) schedule( runtime )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      for (int k=0; k<PS1; k++)  {  z0 = amr->patch[0][lv][PID]->EdgeL[2] + k*dh + 0.5*dh_sub;
      for (int j=0; j<PS1; j++)  {  y0 = amr->patch[0][lv][PID]->EdgeL[1] + j*dh + 0.5*dh_sub;
      for (int i=0; i<PS1; i++)  {  x0 = amr->patch[0][lv][PID]->EdgeL[0] + i*dh + 0.5*dh_sub;

         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] = 0.0;

         for (int kk=0; kk<NSub; kk++)    {  z = z0 + kk*dh_sub;
         for (int jj=0; jj<NSub; jj++)    {  y = y0 + jj*dh_sub;
         for (int ii=0; ii<NSub; ii++)    {  x = x0 + ii*dh_sub;

            Init_Function_Ptr( fluid_sub, x, y, z, Time[lv] );

            for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] += fluid_sub[v];

         }}}

         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] *= _NSub3;

//       check minimum density and pressure
         fluid[DENS] = FMAX( fluid[DENS], (real)MIN_DENS );
         fluid[ENGY] = CPU_CheckMinPresInEngy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                               Gamma_m1, _Gamma_m1, MIN_PRES );

//       calculate the dual-energy variable (entropy or internal energy)
#        if   ( DUAL_ENERGY == DE_ENPY )
         fluid[ENPY] = CPU_Fluid2Entropy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY], Gamma_m1 );
#        elif ( DUAL_ENERGY == DE_EINT )
#        error : DE_EINT is NOT supported yet !!
#        endif

//       floor and normalize passive scalars
#        if ( NCOMP_PASSIVE > 0 )
         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  fluid[v] = FMAX( fluid[v], TINY_NUMBER );

         if ( OPT__NORMALIZE_PASSIVE )
            CPU_NormalizePassive( fluid[DENS], fluid+NCOMP_FLUID, PassiveNorm_NVar, PassiveNorm_VarIdx );
#        endif

         for (int v=0; v<NCOMP_TOTAL; v++)   amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i] = fluid[v];

      }}}
   } // if ( NSub > 1 )

   else // without sub-sampling
   {
//#     pragma omp parallel for private( fluid, x, y, z ) schedule( runtime )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      for (int k=0; k<PS1; k++)  {  z = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh;
      for (int j=0; j<PS1; j++)  {  y = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh;
      for (int i=0; i<PS1; i++)  {  x = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh;

         Init_Function_Ptr( fluid, x, y, z, Time[lv] );

//       check minimum density and pressure
         fluid[DENS] = FMAX( fluid[DENS], (real)MIN_DENS );
         fluid[ENGY] = CPU_CheckMinPresInEngy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                               Gamma_m1, _Gamma_m1, MIN_PRES );

//       calculate the dual-energy variable (entropy or internal energy)
#        if   ( DUAL_ENERGY == DE_ENPY )
         fluid[ENPY] = CPU_Fluid2Entropy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY], Gamma_m1 );
#        elif ( DUAL_ENERGY == DE_EINT )
#        error : DE_EINT is NOT supported yet !!
#        endif

//       floor and normalize passive scalars
#        if ( NCOMP_PASSIVE > 0 )
         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  fluid[v] = FMAX( fluid[v], TINY_NUMBER );

         if ( OPT__NORMALIZE_PASSIVE )
            CPU_NormalizePassive( fluid[DENS], fluid+NCOMP_FLUID, PassiveNorm_NVar, PassiveNorm_VarIdx );
#        endif

         for (int v=0; v<NCOMP_TOTAL; v++)   amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i] = fluid[v];

      }}}
   } // if ( NSub > 1 ) ... else ...

} // FUNCTION : Hydro_Init_StartOver_AssignData



#endif // #if ( MODEL == HYDRO )