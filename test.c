
 
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
  end_c_loop(c,t)

  xmin = xmin*multif1_2;
  xmax = xmax*multif1_2;
  ymin = ymin*multif1_2;
  ymax = ymax*multif1_2;
  zmin = zmin*multif1_2;
  zmax = zmax*multif1_2;
  int xD = floor((xmax-xmin)/(largestParDia*conversion*multif3));
  int yD = floor((ymax-ymin)/(largestParDia*conversion*multif3));
  int zD = floor((zmax-zmin)/(largestParDia*conversion*multif3));

  double dx = (xmax-xmin)/xD;
  double dy = (ymax-ymin)/yD;
  double dz = (zmax-zmin)/zD;

  demInit(xD, yD, zD, dx, dy, dz);
  fflush(stdout);
  //exit(0);
}

DEFINE_SOURCE(xmom_source,c,t,dS,eqn)
{
  real x[ND_ND];
  real con, source;
  C_CENTROID(x,c,t);
  con = C2*0.5*C_R(c,t)*x[1];
  source = -con*fabs(C_U(c, t))*C_U(c,t);
  dS[eqn] = -2.*con*fabs(C_U(c,t));
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
  return source;
} 

/*Execute at the begining of each iteration before conservation equations are solved
  Porosity and drag force information from DEM is transferred to CFD*/
DEFINE_ADJUST(my_adjust, d)
{
  //reset no of fluid cells in bounding box cells

  //copyDEMInfo();
  //printf("DEF ADJUST ");
  fflush(stdout);
}

/* Execute at end of each CFD timestep and update background cells which is required for DEM
    NOTE: Each timestep may have number of iterations
*/
DEFINE_EXECUTE_AT_END(execute_at_end)
{
  Domain *d;
  Thread *t;
  cell_t c;
  real x[ND_ND]; //Fluent real array for storing centroid coordinates

  d = Get_Domain(1);
  t = Lookup_Thread(d,1);

  //Update boundary cells with fluid velocity
  thread_loop_c(t,d)
  {
     begin_c_loop(c,t)
      C_CENTROID(x,c,t);
      
      int iIndex = ceil(x[0]/domainDx);
      int jIndex = ceil(x[1]/domainDy);
      int kIndex = ceil(x[2]/domainDz);
      int cellIndex = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;

      if(cellIndex > xDiv*yDiv*zDiv){
        printf("ERRO CELL I %d\n",cellIndex);
        exit(0);
      }
      bdBox[cellIndex].noOfFluidCells+=1;
      int newSize = bdBox[cellIndex].noOfFluidCells;
      double preVX = bdBox[cellIndex].fluidVelX;
      double preVY = bdBox[cellIndex].fluidVelY;
      double preVZ = bdBox[cellIndex].fluidVelZ;
      
      bdBox[cellIndex].fluidVelX = (newSize-1)*preVX/newSize + C_U(c,t)*velocityFactor/newSize;
      bdBox[cellIndex].fluidVelY = (newSize-1)*preVY/newSize + C_V(c,t)*velocityFactor/newSize;
      bdBox[cellIndex].fluidVelZ = (newSize-1)*preVZ/newSize + C_W(c,t)*velocityFactor/newSize;

      //bdBox[cellIndex].fluidVolF = C_VOF(c,t);
      bdBox[cellIndex].fluidVolF = FVOLF; //For testing FVOLF=0.5
     end_c_loop(c,t)
  }

  //start DEM simulation
  demLoop();
  //save DEM output data 
  demSave();

  //reset bounding box size to zero
  for(int i=0; i<xDiv*yDiv*zDiv; i++){
    bdBox[i].noOfFluidCells = 0;
    bdBox[i].fluidVelX = 0.0;
    bdBox[i].fluidVelY = 0.0;
    bdBox[i].fluidVelZ = 0.0;
  }
  fflush(stdout);
}

DEFINE_EXECUTE_AT_EXIT(execute_at_exit)
{
  // Delete dynamic memeory
  
  free(sortedList);
  free(sortedParIndex);
  free(cellSE);
    //free(parIndex);
  
  free(parNb);
  free(parNoOfNb);
  free(parCIndex);
  
  
  free(uVec);
  free(ipRVec);
  free(jpRVec);
  free(bdBox);
  free(particle);
  printf("EXECUTE AT EXIT\n");
  fflush(stdout);
}

/* Rotate the capsule for a given speed*/
// DEFINE_CG_MOTION(piston,dt,vel,omega,time,dtime)
//  {
//     omega[2]=50;
//  } 
