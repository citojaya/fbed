#include "common.h"

//allocated = 0;
/*Read particle information and domain information
*/
void demInit(){
    updateDPM = 0;
    time_count = 0;// time counter used for saving output file
    particle_counter = 0; //counter for keeping the track of number of particles
    demTime = 0.0;

    if(LogFile){
        fclose(LogFile);
    }
    LogFile = fopen("logfile.log", "a");
    if(sizeof(walls) != 0){
        //printf("SIZE OF WALLS %d\n",sizeof(walls));
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
    } 


    //Read material data
    readInput("infile", &parArraySize, &dens, &ymod, &pois, &sfc, &rec, &dmpn, &rf, &cyldia, &timeStep, 
            &noOfWalls, &updateDPM);
    allocate();
    //Read particle-wall contact surfaces
    readWalls("infile", walls);
    //printf("TIME STEP %lf\n",timeStep/timeFactor);

    writeLogNum("logfile2.log","Update DPM ",updateDPM);
    writeLogNum("logfile2.log","Particle Array Size ",parArraySize);
    writeLogNum("logfile2.log","density ",dens);
    writeLogNum("logfile2.log","Youngs Modulus ",ymod);
    writeLogNum("logfile2.log","Timestep ",timeStep);
    // Read particle information
    //diaInput("pardia", demPart, &np);  
}

/* Allocate arrays */
void allocate(){
    uVec = allocateDoubleArray(DIM);
    ipRVec = allocateDoubleArray(DIM);
    jpRVec = allocateDoubleArray(DIM);
    ijVec = allocateDoubleArray(DIM);
    //tempVec = allocateDoubleArray(dim);
    rotVel = allocateDoubleArray(DIM);
    ipCntPntVel = allocateDoubleArray(DIM);
    jpCntPntVel = allocateDoubleArray(DIM);
    cntPntVel = allocateDoubleArray(DIM);

    //parNb = allocateIntArray(np*nbSize); //neighbourlist of each particle
    //parNoOfNb = allocateIntArray(np); //no of neighbours in each particle
    bdBox = allocateBdBoxArray(xDiv*yDiv*zDiv); //bounding box array
    
    for(int i=0; i<xDiv*yDiv*zDiv; i++){
        //bdBox[i].noOfFaces = 0;
        bdBox[i].noOfParticles = 0;
        //bdBox[i].noOfFluidCells = 0;
        bdBox[i].fluidVelX = 0.0;
        bdBox[i].fluidVelY = 0.0;
        bdBox[i].fluidVelZ = 0.0;
        bdBox[i].fluidVolF = 0.0;
        bdBox[i].dragFX = 0.0;
        bdBox[i].dragFY = 0.0;
        bdBox[i].dragFZ = 0.0;
    }
    
    demPart = allocatePar(parArraySize);
    for(int i=0; i<parArraySize; i++){
        demPart[i].dt = 0.0;
        demPart[i].currentTime = 0.0;
        demPart[i].nrmDisp = 0.0;
        demPart[i].noOfNeigh = 0;
        demPart[i].pos = allocateDoubleArray(DIM);
        demPart[i].angVel = allocateDoubleArray(DIM);
        demPart[i].vel = allocateDoubleArray(DIM);
        demPart[i].hisDisp = allocateDoubleArray(DIM);

        demPart[i].hisDisp[0]  = 0.0;
        demPart[i].hisDisp[1]  = 0.0;
        demPart[i].hisDisp[2]  = 0.0;   

        demPart[i].angVel[0] = 0.0;
        demPart[i].angVel[1] = 0.0;
        demPart[i].angVel[2] = 0.0;

        demPart[i].force = allocateDoubleArray(DIM);
        demPart[i].momentum = allocateDoubleArray(DIM);
        demPart[i].momentum[0] = 0.0;
        demPart[i].momentum[1] = 0.0;
        demPart[i].momentum[2] = 0.0;
        demPart[i].faceNode1 = allocateDoubleArray(DIM);
        demPart[i].faceNode2 = allocateDoubleArray(DIM);
        demPart[i].faceNode3 = allocateDoubleArray(DIM);
        demPart[i].surfNorm = allocateDoubleArray(DIM);
        demPart[i].displacement = 0.0;

    }


    walls = allocateIntArray(noOfWalls);
    //dpmList = (Tracked_Particle)malloc(2*sizeof(Tracked_Particle));

}

/* Generate cells for the problem domain which is used by particle and fluent cells. 
These cells are used to reduce the number of searches when particle and fluid cells exchange 
information 
*/
void buildDEMCFDCellMap(){
    printf("buildDEMCFDCellMap\n");
}


/*
Add contact faces to bounding box
*/
void addFaceToBdBox(){
    // for(int i=0; i<xDiv*yDiv*zDiv; i++){
    //     bdBox[i].noOfFaces = 0;
    // }  
    int noOfF=0;
    Domain *d = Get_Domain(1);
    
    for(int i=0; i<noOfWalls; i++){
        Thread *tf;
        face_t f;
        Node *node;
        tf = Lookup_Thread(d, walls[i]);
        int n;
        begin_f_loop(f, tf)//Loop over faces
        {
            real f_cent[ND_ND];
            F_CENTROID(f_cent,f,tf);
            //printf("CENTROID %lf,%lf,%lf\n",f_cent[0],f_cent[1],f_cent[2]);
            int iIndex = ceil((f_cent[0]*lengthFactor-xmin)/domainDx);
            int jIndex = ceil((f_cent[1]*lengthFactor-ymin)/domainDy);
            int kIndex = ceil((f_cent[2]*lengthFactor-zmin)/domainDz);
            int cellIndex = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;
            
            // if(bdBox[cellIndex].noOfFaces > NO_OF_FACES){
            //     printf("CELL INDEX %d\n",cellIndex);
            //     printf("bdBox[cellIndex].noOfFaces > NO_OF_FACES\n");
            // }
            // else{
            //     bdBox[cellIndex].surfaces[bdBox[cellIndex].noOfFaces] = f;
            //     bdBox[cellIndex].surfaceThread = tf;

            //     bdBox[cellIndex].noOfFaces++;
            // }
        //    f_node_loop(f,tf,n)//Loop over nodes
        //     {
        //         node = F_NODE(f,tf,n);

        //         //insert face to cell
        //         bdBox[cellIndex].surfaces[]
        
        //     }
            
            noOfF++;
        }
        end_f_loop(f,tf)
        
    }
    printf("FACES %d\n",noOfF);
}

/*
Insert particle to BdBox
param:
pI - particle index
cI - cell index
*/
void insertToBdBox(int p, int cI){
    bdBox[cI].parts[bdBox[cI].noOfParticles] = p;
    bdBox[cI].noOfParticles++;
    if(bdBox[cI].noOfParticles > NO_OF_PARTICLES_IN_BDCELL){
        printf("bdBox[cI].noOfParticles > NO_OF_PARTICLES_IN_BDCELL\n");
    }
}

/*
Delete particle from bdBox
param:
pI - particle index
cI = cell index
*/
void deleteParticle(int p, int cI){
    for(int i=0; i<bdBox[cI].noOfParticles; i++){
        int np = bdBox[cI].parts[i];
        if(p == np){
            bdBox[cI].parts[i] = bdBox[cI].parts[bdBox[cI].noOfParticles-1];
            bdBox[cI].noOfParticles--;
            break;
        }
    }
}


/*
Initially neighbourlist array is empty. Fill neighbourlist. 
*/
void initialize(real *nbList, int *parIndex, int *cellSE, int np,
    real *pos, real *parDia)
{
 
}

/*
Assign scale factors
*/
void setReduceUnits()
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
	powerFactor = 6.0/(pow(gravity,1.5)*PI*pow((real)refLength,3.5)*refDensity);
	velocityFactor = 1.0/sqrt(refLength*gravity);
	accFactor = 1.0/gravity;
	angVelFactor = sqrt(refLength/gravity);
	angAccFactor = refLength/gravity;
	freqFactor = sqrt(refLength/gravity);
	inertiaFactor = 6.0/(PI*pow(refLength,5)*refDensity);

    cutGap = CUTGAP*conversion*lengthFactor;

    dsmaxCff = sfc*(2.0-pois)/(2.0*(1.0-pois));
    writeLogNum("logfile2.log"," DS MAX",dsmaxCff);
    dti = 0.0;
    dd = 0.0;
    dsmax = 0.0;

    // scale particle properties
    ymod = ymod*pressureFactor;
    cyldia = cyldia*lengthFactor;
	timeStep = timeStep*timeFactor;
	maxTime	= maxTime*timeFactor;

    elasticMod = ymod/(1.0-pow(pois,2));

    for (int i=0; i<np; i++){
        // demPart[i].pos[0] = demPart[i].pos[0]*lengthFactor;
        // demPart[i].pos[1] = demPart[i].pos[1]*lengthFactor;
        // demPart[i].pos[2] = demPart[i].pos[2]*lengthFactor;
        demPart[i].dt = timeStep;
        // demPart[i].dia = demPart[i].dia*lengthFactor;
        // demPart[i].mass = (4.0/3.0)*PI*pow((0.5*demPart[i].dia),3.0)*dens*densityFactor;
        //demPart[i].inert = 2.0*demPart[i].mass*pow(0.5*demPart[i].dia,2)/5.0;        
    }

    //Find allowed displacment for neighbourlist update
    rIn = largestParDia*conversion*lengthFactor;
    rOut = 1.55*rIn; //By definition
    allowedDisp = 0.5*(rOut-rIn);

    //Adjust boundary limits
    xmin = xmin*lengthFactor;
    ymin = ymin*lengthFactor;
    zmin = zmin*lengthFactor;
    xmax = xmax*lengthFactor;
    ymax = ymax*lengthFactor;
    zmax = zmax*lengthFactor;

    domainDx = domainDx*lengthFactor;
    domainDy = domainDy*lengthFactor;
    domainDz = domainDz*lengthFactor;

    cellRadius = 0.5*sqrt(domainDx*domainDx+domainDy*domainDy);
    writeLogNum("logfile2.log","CEll R ",cellRadius/lengthFactor);
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
Allocate a real* type array
return: real* 
*/
real *allocateDoubleArray(int size)
{
    real *val = (real*)malloc(size*sizeof(real));
    memset(val,0.0,size*sizeof(real));
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
struct demParticle *allocatePar(int np)
{
    struct demParticle *par = (struct demParticle*)malloc(np*sizeof(struct demParticle));
    return par;
}

