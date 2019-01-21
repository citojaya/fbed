
#include "common.h"
#include "surf.h"
#include "mem.h"

#define C2 100.0

/* Initialize dem information
Executes only once at the begining of FLUEMT*/
DEFINE_INIT(my_init, d)
{ 
  xmax = 0.0, ymax=0.0, zmax=0.0;
  xmin = 0.0, ymin=0.0, zmin=0.0;
 
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

  xmin = xmin - offset*conversion;
  xmax = xmax + offset*conversion;
  ymin = ymin - offset*conversion;
  ymax = ymax + offset*conversion;
  zmin = zmin - offset*conversion;
  zmax = zmax + offset*conversion;
  xDiv = floor((xmax-xmin)/(largestParDia*conversion*multif3));
  yDiv = floor((ymax-ymin)/(largestParDia*conversion*multif3));
  zDiv = floor((zmax-zmin)/(largestParDia*conversion*multif3));

  domainDx = (xmax-xmin)/xDiv;
  domainDy = (ymax-ymin)/yDiv;
  domainDz = (zmax-zmin)/zDiv;

  demInit();
  fflush(stdout);
}

/*
Initialize DPM particle information
Executes only once on injected particles
*/
DEFINE_DPM_INJECTION_INIT(solid_paritcles, I)
{
  printf("TRACKED PARTICLE\n");
  Injection *I2;
  Injection *Ilist = Get_dpm_injections();
  //Get_dmp_injections();
  np = 0;
  loop(I2,Ilist)
  {
    Particle *p;
    loop(p,I2->p_init)
    {
      //Set particle mass according to DEM density input
      //demPart[p->part_id].cfdp = p;
      P_MASS(p) = (4.0/3.0)*PI*pow((0.5*P_DIAM(p)),3.0)*dens;
      printf("Particle Mass %lf\n",P_MASS(p));
      printf("INITIAL PAR POS %lf,%lf,%lf\n",P_POS(p)[0],P_POS(p)[1],P_POS(p)[2]);
      np++;
    }
  }

  //Setup DEM scaling 
  setReduceUnits();

  printf("XMIN, YMIN, ZMIN %lf,%lf,%lf\n", xmin/lengthFactor,ymin/lengthFactor,zmin/lengthFactor);
  printf("xDiv,yDiv,zDiv %d,%d,%d\n", xDiv,yDiv,zDiv);
  printf("DOMAIN DX, DY, DY %lf,%lf,%lf\n", domainDx/lengthFactor,domainDy/lengthFactor,domainDz/lengthFactor);

  //Assign to bDBox cells
  //addToBdBox();

  //addFaceToBdBox();
  

  //Set DEM gravitational force
  //assignGravity();
  
  //Update neighbour list
  //updateNeighbourList();
  fflush(stdout);
}

/*
Executes at every DPM time step
This macro is used to call DEM functions which calculate particle-particle contact forces 
*/
DEFINE_DPM_SOURCE(particle_source, c, t, S, strength, p)
{ 
  long int ip = p->part_id;
  //Reset force vector
  demPart[p->part_id].force[0] = 0.0;
  demPart[p->part_id].force[1] = -P_MASS(p)*massFactor; //gravitational force
  char str[10];
  
  //writeLog("logfile2.log","DEFINE_DPM_SOURCE Particle Time \n",5.0);
  demPart[p->part_id].force[2] = 0.0;

  int i;
  int n;
  
  Node *node;
  Thread *tf;
  face_t f;
  cell_t cc = P_CELL(p);
  //Thread *ft = 
 
  //start DEM calculations
  demLoop(p);

  //Find contact with surface
  c_face_loop(c,t,i)//Loop over faces
  {
    int t_id = THREAD_ID(C_FACE_THREAD(c,t,i));
    for(int n=0; n<noOfWalls; n++)
    {
      if(t_id == walls[n])
      {
        f = C_FACE(c,t,i);
        tf = C_FACE_THREAD(c,t,i);

        //printf("FACES -------------%d\n",i);
        //int face_t_id = THREAD_ID(C_FACE_THREAD(cc,t,i));
        int j=0;
        f_node_loop(f,tf,n)//Loop over face nodes
        {
          node = F_NODE(f,tf,n);
          if(j==0)
          {
            demPart[p->part_id].faceNode1[0] = NODE_X(node)*lengthFactor;
            demPart[p->part_id].faceNode1[1] = NODE_Y(node)*lengthFactor;
            demPart[p->part_id].faceNode1[2] = NODE_Z(node)*lengthFactor;
          }
          if(j==1)
          {
            demPart[p->part_id].faceNode2[0] = NODE_X(node)*lengthFactor;
            demPart[p->part_id].faceNode2[1] = NODE_Y(node)*lengthFactor;
            demPart[p->part_id].faceNode2[2] = NODE_Z(node)*lengthFactor;
          }
          if(j==2)
          {
            demPart[p->part_id].faceNode3[0] = NODE_X(node)*lengthFactor; 
            demPart[p->part_id].faceNode3[1] = NODE_Y(node)*lengthFactor;
            demPart[p->part_id].faceNode3[2] = NODE_Z(node)*lengthFactor;
          }
          j++;
       }

        real *uVec = allocateDoubleArray(DIM);
        getUnitVector(demPart[ip].faceNode1,demPart[ip].faceNode2,demPart[ip].faceNode3,uVec);
        real *parPos = allocateDoubleArray(DIM);
        parPos[0] = P_POS(p)[0]*lengthFactor;
        parPos[1] = P_POS(p)[1]*lengthFactor;
        parPos[2] = P_POS(p)[2]*lengthFactor;
        real gap = getOverlap(parPos,P_DIAM(p)*lengthFactor, demPart[ip].faceNode1,  
                              demPart[ip].faceNode2, demPart[ip].faceNode3,  uVec);
        
        if(gap < 0) //If contact exists calculate contact force
        {
          surfaceContactForce(p, -gap, uVec);
        }
        free(uVec);
        free(parPos);

      }
    }
  }
  //Calculate interparticle contact force
  //partContactForce(p, nrmDisp)
 
  //Update particle position
  //move(p);

  
  fflush(stdout);
}

DEFINE_DPM_SCALAR_UPDATE(scalar_update, c, t, initialize, p)
{
  //writeLog("logfile2.log","DEFINE_DPM_SCALAR_UPDATE\n",1.0);
  int iIndex = ceil((P_POS(p)[0]*lengthFactor-xmin)/domainDx);
  int jIndex = ceil((P_POS(p)[1]*lengthFactor-ymin)/domainDy);
  int kIndex = ceil((P_POS(p)[2]*lengthFactor-zmin)/domainDz);
  int cellIndex = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;
  deleteParticle(p, cellIndex);//delete exist particle 
  insertToBdBox(p, cellIndex);//add new particle
}

/*
When particles are in contact with surface this macro is executed by FLUENT
Loop through all particles and update contact surface normal vector and contact surface node positions
Only applicable when contact surface is triangular
*/
DEFINE_DPM_BC(bc_reflect, p, t, f, f_normal, dim)
{
  long int ip = p->part_id;
  //printf("DEFINE_DPM_BC **********************************************\n");
  //printf("PARTICLE ID %ld\n",ip);
  //printf("PAR DIA %llf \n",P_DIAM(p));
  //printf("POS %llf, %llf, %llf\n",P_POS(p)[0],P_POS(p)[1], P_POS(p)[2]);

  //test(p,f,t, f_normal);
     
  fflush(stdout);
  //return PATH_ABORT;
  return PATH_ACTIVE;
}

DEFINE_EXECUTE_AT_END(execute_at_end)
{
  //writeLog("logfile2.log","EXECUTE AT END\n",1.0);
  demSave();

}

DEFINE_ADJUST(define_adjust, d)
{
  //writeLog("logfile2.log","DEFINE_ADJUST\n",2.0);
  //demLoop();
}

DEFINE_EXECUTE_AT_EXIT(execute_at_exit)
{
  // Delete dynamic memeory
  fclose(LogFile);
  free(walls);
  free(uVec);
  free(ipRVec);
  free(jpRVec);

  free(bdBox);
    
  free(ijVec);
    //free(tempVec);
  free(rotVel);
  free(ipCntPntVel);
  free(jpCntPntVel);
  free(cntPntVel);

  for(int i=0; i<np; i++){
    free(demPart[i].pos);
    free(demPart[i].angVel);
    free(demPart[i].vel);
    free(demPart[i].hisDisp);
    free(demPart[i].force);
    free(demPart[i].momentum);
    free(demPart[i].faceNode1);
    free(demPart[i].faceNode2);
    free(demPart[i].faceNode3);
    free(demPart[i].surfNorm);
  }
  free(demPart);
  printf("EXECUTE AT EXIT\n");
  fflush(stdout);
}

/* Rotate the capsule for a given speed*/
// DEFINE_CG_MOTION(piston,dt,vel,omega,time,dtime)
//  {
//     omega[2]=50;
//  } 
