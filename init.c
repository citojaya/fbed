#include "common.h"

//allocated = 0;
/*Read particle information and domain information
*/
void demInit(int xD, int yD, int zD, double dx, double dy, double dz){
    xDiv = xD;
    yDiv = yD;
    zDiv = zD;

    domainDx = dx;
    domainDy = dy;
    domainDz = dz;
    demTime = 0.0;

    //struct Particle *par;
   
    //Read material data
    readData("infile", &np, &dens, &ymod, &pois, &sfc, &rec, &dmpn, &cyldia, &timeStep);
    allocate();
 
    //printf("Initialized %d\n",initialized);
    printf("%d\n",np);
 
    // Read particle information
    diaInput("pardia", particle, &np);   
    
    //Setup scaling
    setScaleFactors();

    //generate bounding box for CFD and DEM sharing
    //boundingBox(d);

    for (int i=0; i<np; i++){
        printf("PAR %lf,%lf,%lf\n",particle[i].posX/(lengthFactor),particle[i].posY/(lengthFactor),particle[i].posZ/(lengthFactor));
    }
}

/* Allocate arrays */
void allocate(){
    sortedList = allocateDoubleArray(np*2); //sorted neighbourlist
    sortedParIndex = allocateIntArray(np*2); //sorted particle indices 
    // parPosX = allocateDoubleArray(np);
    // parPosY = allocateDoubleArray(np);
    // parPosZ = allocateDoubleArray(np);
    // parMass = allocateDoubleArray(np);
    // parInert = allocateDoubleArray(np);

    uVec = allocateDoubleArray(dim);
    ipRVec = allocateDoubleArray(dim);
    jpRVec = allocateDoubleArray(dim);

    cellSE = allocateIntArray(np*2);
//    parDia = allocateDoubleArray(np);
    cellSE = allocateIntArray(np*3); //to identify cell start and end positions
    parNb = allocateIntArray(np*nbSize); //neighbourlist of each particle
    parNoOfNb = allocateIntArray(np); //no of neighbours in each particle
    parCIndex = allocateIntArray(np*2);// start and end of particle position in sortedList array
    bdBox = allocateBdBoxArray(xDiv*yDiv*zDiv); //bounding box array
    
    for(int i=0; i<xDiv*yDiv*zDiv; i++){
        bdBox[i].noOfFluidCells = 0;
        bdBox[i].fluidVelX = 0.0;
        bdBox[i].fluidVelY = 0.0;
        bdBox[i].fluidVelZ = 0.0;
        bdBox[i].fluidVolF = 0.0;
        bdBox[i].dragFX = 0.0;
        bdBox[i].dragFY = 0.0;
        bdBox[i].dragFZ = 0.0;
    }
    //par = allocatePar(np);
    particle = allocatePar(np);
    // for(int i=0; i<np; i++){
    //     particle[i].posX = 0.0;
    //     particle[i].posY = 0.0;
    //     particle[i].posZ = 0.0;

    // }

}

/* Allocate array for domain bounding box and allocate CFD cells to bounding box cells
   param: *d Fluent Domain pointer
*/ 
// void boundingBox(Domain *d){
//     cell_t c; //Fluent cell
//     real x[ND_ND]; //Fluent real array for storing centroid coordinates
//     Thread *t = Lookup_Thread(d,1); //Fluent cell thread
//     begin_c_loop(c,t)
//       C_CENTROID(x,c,t);

//       insertCellToBBox(c, x);
//     end_c_loop(c,t)
// }

/*
Allocate CFD cells to bounding box
param: c Fluent fluid cell, x[] fluid cell coordinates
*/
// void insertCellToBBox(cell_t c, real x[]){
//     int iIndex = ceil(x[0]/domainDx);
//     int jIndex = ceil(x[1]/domainDy);
//     int kIndex = ceil(x[2]/domainDz);
//     int cellIndex = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;
//     if(cellIndex > xDiv*yDiv*zDiv){
//         printf("ERROR: cellIndex>xDiv*yDiv*zDiv %d\n",cellIndex);
//         exit(0);
//     }
//     else{
//         bdBox[cellIndex].fluid_cell[bdBox[cellIndex].noOfFluidCells] = c;
//         bdBox[cellIndex].noOfFluidCells+=1;
//     }
//     if(bdBox[cellIndex].noOfFluidCells>NO_OF_FLUID_CELLS){
//         printf("Cell %d\n",bdBox[cellIndex].noOfFluidCells);
//         exit(0);
//     }
// }

/*
Update average fluid velocity of bounding box cells
*/
// void updateBBFluidVel(){
//     for (int i=0; i<xDiv*yDiv*zDiv){

//     }
// }



/* Generate cells for the problem domain which is used by particle and fluent cells. 
These cells are used to reduce the number of searches when particle and fluid cells exchange 
information 
*/
void buildDEMCFDCellMap(){
    printf("buildDEMCFDCellMap\n");
}

/*
Initially neighbourlist array is empty. Fill neighbourlist. 
*/
void initialize(double *nbList, int *parIndex, int *cellSE, int np,
    double *pos, double *parDia)
{
    int j=0;
    for (int i=0; i<np; i++){
        sortedList[j] = pos[i] - 0.5*parDia[i];
        sortedList[j+1] = pos[i] + 0.5*parDia[i];
        parIndex[j] = i; //particle index
        parIndex[j+1] = i; //particle index
        cellSE[j] = 1; //cell start
        cellSE[j+1] = 2; //cell end
        parCIndex[j] = j; //cell start index for each particle
        parCIndex[j+1] = j+1; //cell end index for each particle
        j += 2;
    }
}

/*
Assign scale factors
*/
void setScaleFactors()
{
    //Scale factors for reduced units
	refLength = largestParDia*conversion;
	refDensity = largestParDensity;
	lengthFactor = 1.0/refLength;
	volumeFactor = pow(lengthFactor,3);
	massFactor = 6.0/(PI*pow(refLength,3)*refDensity);
	timeFactor = sqrt(gravity/refLength);
	densityFactor = 6.0/(PI*refDensity);
	forceFactor = 6.0/(gravity*PI*pow(refLength,3)*refDensity);
	pressureFactor = 6.0/(gravity*PI*refLength*refDensity);
	StressFactor = pressureFactor;
	energyFactor = 6.0/(gravity*PI*pow(refLength,4)*refDensity);
	momentFactor = energyFactor;
	powerFactor = 6.0/(pow(gravity,1.5)*PI*pow((double)refLength,3.5)*refDensity);
	velocityFactor = 1.0/sqrt(refLength*gravity);
	accFactor = 1.0/gravity;
	angVelFactor = sqrt(refLength/gravity);
	angAccFactor = refLength/gravity;
	freqFactor = sqrt(refLength/gravity);
	inertiaFactor = 6.0/(PI*pow(refLength,5)*refDensity);

    cutGap = CUTGAP*conversion*lengthFactor;

    dsmaxCff = sfc*(2.0-pois)/(2.0*(1.0-pois));

    // scale particle properties
    dens = dens*densityFactor;
    ymod = ymod*pressureFactor;
    cyldia = cyldia*lengthFactor;
	timeStep = timeStep*timeFactor;
	maxTime	= maxTime*timeFactor;

    for (int i=0; i<np; i++){
        particle[i].posX = particle[i].posX*conversion*lengthFactor;
        particle[i].posY = particle[i].posY*conversion*lengthFactor;
        particle[i].posZ = particle[i].posZ*conversion*lengthFactor;
        particle[i].dia = particle[i].dia*conversion*lengthFactor;
        particle[i].mass = (4.0/3.0)*PI*pow((0.5*particle[i].dia),3.0)*dens;
        particle[i].inert = 2.0*particle[i].mass*pow(0.5*particle[i].dia,2)/5.0;        
    }
    //printf("refLength  %lf\n ",refLength);
}

/*
void allocateMat(struct MatType *m){
    printf("allocation.c\n");
    //struct MatType *mt = (struct MatType *)malloc(num_of_mats*sizeof(struct MatType));
    //m = mt;    
}*/

/*
Allocate an integer type array
return: int* 
*/
int *allocateIntArray(int size)
{
    int *val = (int*)malloc(size*sizeof(int));
    memset(val,0,size*sizeof(int));
    return val;
}

/*
Allocate a double* type array
return: double* 
*/
double *allocateDoubleArray(int size)
{
    double *val = (double*)malloc(size*sizeof(double));
    memset(val,0.0,size*sizeof(double));
    return val;
}

/*
Allocate a char* type array
return: char* 
*/
char *allocateCharArray(int size)
{
    char *val = (char*)malloc(size*sizeof(char));
    return val;
}

/*
Allocate bounding box type array
return: BdBox*
*/
struct BdBox *allocateBdBoxArray(int size)
{
    //printf("BD BOX SIZE %d\n",size);
    struct BdBox *bdB = (struct BdBox*)malloc(size*sizeof(struct BdBox));
    return  bdB;
}

/*
Allocate particle array
*/
struct Particle *allocatePar(int np)
{
    struct Particle *par = (struct Particle*)malloc(np*sizeof(struct Particle));
    return par;
}

