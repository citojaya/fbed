
#ifndef COMMON_H
#define COMMON_H

#include "udf.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*--------------Common definitions---------------*/
#define NUM_MAT 2    // two types material
#define PI  3.1415926f
#define gravity 9.8f //ms^-2
#define CUTGAP 0.1f //mm
#define conversion 1.0e-3f //convert length values to meters


/*------------ Particle information -------------*/
#define largestParDia 1.5f //mm
#define largestParDensity 2800.0f //kgm^-3
#define multif3 3 //multification factor used in bounding box divisions
#define multif1_2 1.2f //mulficiation factor used to adjust boundary size
#define arrSize 5
#define dim 3 // 3D problem
#define nbSize 5 //size of neighbourlist
#define NO_OF_FLUID_CELLS 20 //number of fluid cells in a bounding box

//Testing parameters
#define FVOLF 0.5f //fluid volume fraction
#define BETA 0.4f //interphase momentum exchange coefficient
#define yMinBottom  0.0f//Y coordinate of bottom surface

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
double *sortedList; //sorted array for recording particle start and end positions
int np; //number of particles 
unsigned int initialized;
int *sortedParIndex; //sorted array for recording particle index 
int *cellSE; //keeps a record of start and end positions for all particle (start=1, end=2)
int *parNb; //partilce neighbourlist
int *parCIndex; //particle start and end positons of sortedList
int *parNoOfNb; //no of neighbours in a particle
//double *parPosX, *parPosY, *parPosZ, *parDia, *parInert, *parMass, cutGap; //particle parameters
double cutGap; //particle parameters
double refLength,refDensity,lengthFactor,volumeFactor,massFactor,timeFactor,
	densityFactor, forceFactor, pressureFactor, StressFactor, energyFactor, momentFactor,
	powerFactor, velocityFactor, accFactor, angVelFactor, angAccFactor, freqFactor, inertiaFactor;
int xDiv, yDiv, zDiv; //Number of divisions in orthogonal domain cells
double domainDx, domainDy, domainDz; //Domain cell size
//double xmax, ymax, zmax; //Max values of domain boundary 
//double xmin, ymin, zmin; //Min values of domain boundary 

double *uVec, *ipRVec, *jpRVec;
double dens, ymod, pois, sfc, rec, dmpn; //particle material property
double dsmaxCff; //used in contact force calculation
double cyldia; //cylinder diameter
double timeStep, demTime, maxTime;
struct BdBox *bdBox;
struct Particle *particle;

//*--- Boundary Condition ---*// 
struct CylinderBC { 
	double cir;     // the position of axial and radius (X, Y)
	double R;        // radius
	double Tw;       // top position
	double Bw;       // bottom position
	double topv, btmv;            
};

// material properties
struct MatType {
	double density;
	double emod, ymod, pois;
	double dmpn;
	double sfrc, rfrc;
	double yldp;        
};

//Bounding box which holds CFD cells and DEM particles
struct BdBox{
	double fluidVelX, fluidVelY, fluidVelZ;
	double pGradX, pGradY, pGradZ;
	double fluidVolF;
	double dragFX, dragFY, dragFZ;

	int noOfFluidCells; //number of fluid cells
	int noOfParticles; //number of DEM particles
};

//Particle
struct Particle{
	double posX, posY, posZ, dia, inert, mass, vX, vY, vZ, forceX, forceY, forceZ;
};

//parIndex		- Array index of particles
//parIndexNL 	- Array for storing particle indices for corresponding neighbour list
//parCord 		- cordinates of particles
//parDia		- particle diameter

//static int num_of_mats = 5;
//static int np; //Total number of particles in the system

void allocateMat(struct MatType *mt);
int *allocateIntArray(int size);
double *allocateDoubleArray(int size);
struct BdBox *allocateBdBoxArray(int size);
struct Particle *allocatePar(int np);
//void insertCellToBBox(cell_t c, real x[]);
void insertParToBBox(int pNo, double x[]);
void readData(char *infile, int *np, double *dens, double *ymod, 
			double *pois, double *sfc, double *rec, double *dmpn, double *cyldia, double *dt);
void findRec(FILE *inFile, char* strDest);
void diaInput(char *diaFile, struct Particle *par, int *np);

void writeTec();
//void demInit(int *xD, int *yD, int *zD, double xmin,double xmax,
//					double ymin,double ymax,double zmin,double zmax);
void demInit(int xD, int yD, int zD, double dx, double dy, double dz);
void buildDEMCFDCellMap();
void copyDEMInfo();
void demLoop();
void demSave();
void allocate();
void move();
void updateForce();
void neighbourContactForce(int pI);
void boundaryContactForce(int pI);
void dragForce(int pI);
void vecAdd(double *v1, double *v2, double *vec);
void crossProd(double *v1, double *v2, double *vec);
void vecSub(double *v1, double *v2, double *vec);
void unitVec(double *v, double *vec);
void sclVecMult(double scl, double *vec);
double vecMag(double *vec);
void projVec(double *v1, double *n, double *vec, int type);

void initialize(double *sortedList, int *sortedParIndex, int *cellSE, int np,
    double *pos, double *parDia);
void insertionSort(double *nbArray, int size, int *parIArray, int *cellSEArray, int firstTime);
void assignNeighbours(double *sortedList, int *sortedParIndex, int *cellSE, int size);
void addNeighbour(int ip, int jp);
void deleteNeighbour(int ip, int jp);
void update(double *pX, int np);

double getCenterDist(int ip, int jp);
void setScaleFactors();
//void boundingBox(Domain *d);
void updateBBFluidVel();

void test();
void run();

#endif 
