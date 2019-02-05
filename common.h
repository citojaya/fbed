
#ifndef COMMON_H
#define COMMON_H

#include "udf.h"
#include "dpm.h"
//#include "dem.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*--------------Common definitions---------------*/
#define NUM_MAT 2    // two types material
#define PI 3.1415926f
#define gravity 9.8f //ms^-2
#define CUTGAP 0.5f //mm
#define conversion 1.0e-3f //convert length values to meters
#define overlapLimit 0.005f //ovelap cuttoff value

/*------------ Particle information -------------*/
#define largestParDia 2.5f //mm
#define largestParDensity 2500.0f //kgm^-3
#define multif3 2 //multification factor used in bounding box divisions
#define multif1_2 1.2f //mulficiation factor used to adjust boundary size
#define offset 10.0f //boundary is adjusted by offset (mm) in order to include additional cells
//#define arrSize 5
#define DIM 3 // 3D problem
#define NO_OF_PARTICLES 100 //Particle array size
#define NBSIZE 50 //size of neighbourlist
#define NO_OF_FLUID_CELLS 20 //number of fluid cells in a bounding box
#define NO_OF_FACES 3 //number of faces contacting with particles
#define NO_OF_PARTICLES_IN_BDCELL 50
//Testing parameters
#define FVOLF 0.5f //fluid volume fraction
#define BETA 0.4f //interphase momentum exchange coefficient
#define yMinBottom  0.0f//Y coordinate of bottom surface

/*----------- Fluid properties -------------------*/
#define VISCOSITY 1.7893e-5f //kg/m-s
#define POROSITY 0.5f
//typedef unsigned int uint;

// err detection utility
#define FPRINTF(a) fprintf a


//Domain *fd; //fluid domain
//Thread *t; //fluid cell thread
/*-----------------------------------------------*/
/*--------------Global declaration---------------*/ 
/*-----------------------------------------------*/
//extern char* genfile;
//static int np;
//real *sortedList; //sorted array for recording particle start and end positions
//real *nrmDisp;
FILE *LogFile;
int time_count;
int particle_counter;
int noOfWalls;
int *walls;
int np, parArraySize; //number of particles, particle array size 
unsigned int initialized;
//int *sortedParIndex; //sorted array for recording particle index 
//int *cellSE; //keeps a record of start and end positions for all particle (start=1, end=2)
//int *parNb; //partilce neighbourlist
//int *parCIndex; //particle start and end positons of sortedList
//int *parNoOfNb; //no of neighbours in a particle
//real *parPosX, *parPosY, *parPosZ, *parDia, *parInert, *parMass, cutGap; //particle parameters
real cutGap; //particle parameters
real cellRadius; //Radius defined by bounding box cell given by 0.5*sqrt(dx^2+dy*2)
real refLength,refDensity,lengthFactor,volumeFactor,massFactor,timeFactor,
	densityFactor, forceFactor, pressureFactor, StressFactor, energyFactor, momentFactor,
	powerFactor, velocityFactor, accFactor, angVelFactor, angAccFactor, freqFactor, inertiaFactor;
int xDiv, yDiv, zDiv; //Number of divisions in orthogonal domain cells
real domainDx, domainDy, domainDz; //Domain cell size
real xmax, ymax, zmax; //Max values of domain boundary 
real xmin, ymin, zmin; //Min values of domain boundary 

/*
rIn = largestParDia
rOut = 1.55*rIn
allowedDisp = (rOut-rIn)/2
*/
real rIn, rOut, allowedDisp; //if particle displacement > allowedDisp -> update neighbourlist

real *uVec, *ipRVec, *jpRVec, *ijVec, *rotVel;
real *ipCntPntVel, *jpCntPntVel, *cntPntVel;
real dens, ymod, pois, sfc, rf, rec, dmpn, elasticMod; //particle material property
real dsmaxCff, dti, dd, dsmax, fdt; //used in contact force calculation
real cyldia; //cylinder diameter
real timeStep, demTime, maxTime;
struct BdBox *bdBox;
struct demParticle *demPart;
//Tracked_Particle *dpmList[2];
int updateDPM; //flag for updating DPM particles when DEFINE_ADJUST is called

//*--- Boundary Condition ---*// 
struct CylinderBC { 
	real cir;     // the position of axial and radius (X, Y)
	real R;        // radius
	real Tw;       // top position
	real Bw;       // bottom position
	real topv, btmv;            
};

// material properties
struct MatType {
	real density;
	real emod, ymod, pois;
	real dmpn;
	real sfrc, rfrc;
	real yldp;        
};

//Bounding box which holds CFD cells and DEM particles
struct BdBox{
	real fluidVelX, fluidVelY, fluidVelZ;
	real pGradX, pGradY, pGradZ;
	real fluidVolF;
	real dragFX, dragFY, dragFZ;

	//int noOfFluidCells; //number of fluid cells
	//int noOfFaces; //number of contact faces
	int noOfParticles; //number of DEM particles
	int parts[NO_OF_PARTICLES_IN_BDCELL];
	face_t surfaces[NO_OF_FACES];
	Thread *surfaceThread;
	//Tracked_Particle *parts[20];
};


//Particle
struct demParticle{
	real dt; //particle time step
	real currentTime; //current time
	real displacement; //if displacement > rMax update neighbourlist
	real dia, inert, mass, nrmDisp;
	real *pos, *angVel, *vel, *hisDisp, *force, *momentum;
	real *faceNode1, *faceNode2, *faceNode3, *surfNorm;
	int neigh[NBSIZE];
	int noOfNeigh;
	Particle *cfdp;
};

//parIndex		- Array index of particles
//parIndexNL 	- Array for storing particle indices for corresponding neighbour list
//parCord 		- cordinates of particles
//parDia		- particle diameter

//static int num_of_mats = 5;
//static int np; //Total number of particles in the system
real solidFraction(int ip);
void writeLogNum(char *infile, char *line, real num);
void test(Tracked_Particle *p, Thread *t);
void allocateMat(struct MatType *mt);
int *allocateIntArray(int size);
real *allocateDoubleArray(int size);

struct BdBox *allocateBdBoxArray(int size);
struct demParticle *allocatePar(int np);
//void insertCellToBBox(cell_t c, real x[]);
real partVol(int p);
void insertToBdBox(int p, int cI);
void addToBdBox();
void addFaceToBdBox();
void readInput(char *infile, int *np, real *dens, real *ymod, 
			real *pois, real *sfc, real *rec, real *dmpn, real *rf, real *cyldia, real *dt, int *nW, int *updateDPM);
void findRec(FILE *inFile, char* strDest);
void diaInput(char *diaFile, struct demParticle *par, int *np);
void readWalls(char *infile, int *walls);

void writeTec();
//void demInit(int *xD, int *yD, int *zD, real xmin,real xmax,
//					real ymin,real ymax,real zmin,real zmax);
void demInit();
void buildDEMCFDCellMap();
void copyDEMInfo();
//void demLoop(Tracked_Particle *p);
void demSave();
void allocate();

void updateForce(Tracked_Particle *p);
void partContactForce(int ip, int jp, real nrmDsp);
void boundaryContactForce(int pI, real *n1, real *n2, real *n3, real *uVec);
//void dragForce(int pI);
void assignGravity();
real getOverlap(real *parPos, real dia, real *n1, real *n2, real *n3, real *uVec);

void neighbourContactForce(int pI);
void surfaceContactForce(int p, real nrmDisp, real *uVec);
void vecAdd(real *v1, real *v2, real *vec);
void crossProd(real *v1, real *v2, real *vec);
void vecSub(real *v1, real *v2, real *vec);
void unitVec(real *v, real *vec);
void sclMult(real scl, real *vec);
void sclVecMult(real scl, real *inVec, real *outVec);
real vecMag(real *vec);
void projVec(real *v1, real *n, real *vec, int type);
real dotProduct(real *v1, real *v2);
void getUnitVector(real *v1, real *v2, real *v3, real *uVec);

void initialize(real *sortedList, int *sortedParIndex, int *cellSE, int np,
    real *pos, real *parDia);
//void insertionSort(real *nbArray, int size, int *parIArray, int *cellSEArray, int firstTime);
//void assignNeighbours(real *sortedList, int *sortedParIndex, int *cellSE, int size);

void addNeighbour(int  ip, int jp);
void updateNeighbourList(int p);
void deleteNeighbour(int ip, int jp);
void update(real *pX, int np);

real getCenterDist(int ip, int jp);
void setReduceUnits();
void updateBBFluidVel();
void insertToBdBox(int p, int cI);
void deleteParticle(int p, int cI);
void forceCalculation(Particle *p);
void updatePosition(Particle *p);
void run();

#endif 
