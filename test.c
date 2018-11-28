
 #include "udf.h"
 #include "common.h"

 #define C2 100.0

 
/* Initialize dem information
Executes only once at the begining of program*/
DEFINE_INIT(my_init, d)
{ 
  double xmax = 0.0, ymax=0.0, zmax=0.0;
  double xmin = 0.0, ymin=0.0, zmin=0.0;
 
  cell_t c;
  real x[ND_ND];
  Thread *t = Lookup_Thread(d,1);
  begin_c_loop(c,t)
    C_CENTROID(x,c,t);
    if(x[0] > xmax) xmax = x[0];
    if(x[1] > ymax) ymax = x[1];
    if(x[2] > zmax) zmax = x[2];
    if(x[0] < xmin) xmin = x[0];
    if(x[1] < ymin) ymin = x[1];
    if(x[2] < zmin) zmin = x[2];
    //printf("XMAX%lf\n",x[1]);
  end_c_loop(c,t)
  printf("MAX%lf,%lf,%lf\n",xmax,ymax,zmax);
  printf("MIN%lf,%lf,%lf\n",xmin,ymin,zmin);
  // thread_loop_c(t,d)
  // {

  // }

  demInit();
  buildDEMCFDCellMap();
  fflush(stdout);
}

DEFINE_SOURCE(xmom_source,c,t,dS,eqn)
{
  real x[ND_ND];
  real con, source;
  C_CENTROID(x,c,t);
  con = C2*0.5*C_R(c,t)*x[1];
  source = -con*fabs(C_U(c, t))*C_U(c,t);
  dS[eqn] = -2.*con*fabs(C_U(c,t));
  //printf("SOURCE\n");
  return source;
} 

DEFINE_SOURCE(ymom_source,c,t,dS,eqn)
{
  real x[ND_ND];
  real con, source;
  C_CENTROID(x,c,t);
  con = C2*0.5*C_R(c,t)*x[1];
  source = -con*fabs(C_U(c, t))*C_U(c,t);
  dS[eqn] = -2.*con*fabs(C_U(c,t));
  //printf("SOURCE\n");
  return source;
} 

DEFINE_SOURCE(zmom_source,c,t,dS,eqn)
{
  real x[ND_ND];
  real con, source;
  C_CENTROID(x,c,t);
  con = C2*0.5*C_R(c,t)*x[1];
  source = -con*fabs(C_U(c, t))*C_U(c,t);
  dS[eqn] = -2.*con*fabs(C_U(c,t));
  //printf("SOURCE\n");
  return source;
} 

/*Execute at the begining of each iteration before conservation equations are solved
  Porosity and drag force information from DEM is transferred to CFD*/
DEFINE_ADJUST(my_adjust, d)
{
  //Update porosity and momentum source term in DEM cells
  copyDEMInfo();

  Thread *t;
  cell_t c;

  //Loop over FLuent cells and update porosity and momentum term
  thread_loop_c(t,d)
  {
    begin_c_loop(c,t)
      //printf("CELL\n");
    end_c_loop(c,t)
  }

  
  fflush(stdout);
}

/* Execute at end of each iteration transfer CFD information to DEM
    NOTE: Each iteration may have number of timesteps
*/
DEFINE_EXECUTE_AT_END(execute_at_end)
{
  Domain *d1, *d2;
  Thread *t1, *t2;
  cell_t c;

  d1 = Get_Domain(1);
  d2 = Get_Domain(2);
  /*Array for storing cell centroid
    real xc[ND_ND];

  
    face_t f;
    Node *node;
    cell_t c;
    int n;
    real xc[ND_ND];

    d = Get_Domain(1);*/  /* mixture domain if multiphase */
  /*  t = Lookup_Thread(d,6);
    printf("DEFINE_EXECUTE_AT_END\n");
 
    int i=0;
*/
  t1 = Lookup_Thread(d1,1);
  t2 = Lookup_Thread(d2,1);
  thread_loop_c(t2,d2)
  {
    begin_c_loop(c,t2)
      //printf("Cell Vel %lf,%lf,%lf\n", C_U(c, t), C_V(c, t),C_W(c, t));
      //printf("Cell Press %lf\n", C_P(c, t));
      C_VOF(c,t2) = 0.4;
      //printf("Cell C_VOF %lf\n", C_VOF(c, t));
      
    end_c_loop(c,t2)
  }

  printf("TIME STEP%lf\n",CURRENT_TIMESTEP);
  thread_loop_c(t1,d1)
  {
    begin_c_loop(c,t1)
      //printf("D1 Cell Vel %lf,%lf,%lf\n", C_U(c, t1), C_V(c, t1),C_W(c, t1));
      //printf("D1 Cell Press %lf\n", C_P(c, t1));
    end_c_loop(c,t1)
  }
  /*
    begin_c_loop(c,t)
    {
      // c_node_loop(c,t,n)
      // {

      // }

      C_CENTROID(xc,c,t);

      printf("Cell loop %lf, %lf, %lf\n",xc[0], xc[1], xc[2]);
      i++;
    }
    end_c_loop(c,t)
  */
    /* DEM iteration*/
    //demLoop();
    printf("No of par %d\n", np);


    //run();
    fflush(stdout);
}

DEFINE_EXECUTE_AT_EXIT(execute_at_exit)
{
  // Delete dynamic memeory
  free(parPosX);
  free(parPosY);
  free(parPosZ);
  free(sortedList);
  free(sortedParIndex);
  free(cellSE);
    //free(parIndex);
  free(parDia);
  free(parNb);
  free(parNoOfNb);
  free(parCIndex);
  free(parInert);
  free(parMass);
  free(uVec);
  free(ipRVec);
  free(jpRVec);
  printf("EXECUTE AT EXIT\n");
  fflush(stdout);
}

/* Rotate the capsule for a given speed*/
// DEFINE_CG_MOTION(piston,dt,vel,omega,time,dtime)
//  {
//     omega[2]=50;
//  } 
