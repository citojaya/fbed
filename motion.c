#include "common.h"

/* 
DEM simulation
After each fluid timestep this method is called by DEFINE_EXECUTE_AT_END
*/
void demLoop(){
    demTime += timeStep;
    updateForce(); //Calculate contact forces 
    move();// update particle position
    //updateNeighList(); //Update neighbour list
}

/* Find contact force */
void updateForce(){
    for (int i=0; i<np; i++){
        //dragForce(i);
        //neighbourContactForce(i);
        boundaryContactForce(i);
    }
}

/* Find drag force on a particle due to fluid flow */
void dragForce(int pI){
    int iIndex = ceil(particle[pI].posX/lengthFactor/domainDx);
    int jIndex = ceil(particle[pI].posY/lengthFactor/domainDy);
    int kIndex = ceil(particle[pI].posZ/lengthFactor/domainDz);
    int cellIndex = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;

    double parVol = PI*pow(particle[pI].dia,3)/6.0; //particle volume

    double dFX = BETA*parVol*(bdBox[cellIndex].fluidVelX-particle[pI].vX)/(1.0-bdBox[cellIndex].fluidVolF);
    double dFY = BETA*parVol*(bdBox[cellIndex].fluidVelY-particle[pI].vY)/(1.0-bdBox[cellIndex].fluidVolF);
    double dFZ = BETA*parVol*(bdBox[cellIndex].fluidVelZ-particle[pI].vZ)/(1.0-bdBox[cellIndex].fluidVolF);

    particle[pI].forceX += dFX;
    particle[pI].forceY += dFY;
    particle[pI].forceZ += dFZ;
    printf("FLUID VEL %lf,%lf,%lf\n",bdBox[cellIndex].fluidVelX,bdBox[cellIndex].fluidVelY,bdBox[cellIndex].fluidVelZ);

    printf("PAR COORD %lf,%lf,%lf\n",particle[pI].posX/lengthFactor,particle[pI].posY/lengthFactor,particle[pI].posZ/lengthFactor);
    printf("DIV\n");
    printf("PAR CELL INDEX %d,%d,%d\n",iIndex,jIndex,kIndex);
}

/*Find contact forces with all boundaries*/
void boundaryContactForce(int pI){
    
    double rStar = 0.5*particle[pI].dia;
	double pXY = particle[pI].posY;
	//double gap = yMinBottom - pXY - rStar;
    double gap = pXY - rStar;
	double sfCoff = sfc;
    printf("demTime%lf\n",demTime/timeFactor);
    printf("PAR COORD %lf,%lf,%lf\n",particle[pI].posX/lengthFactor,particle[pI].posY/lengthFactor,particle[pI].posZ/lengthFactor);
    if(gap < 0){
        uVec[0] = 0.0;
        uVec[1] = -particle[pI].posY/pXY;
        uVec[2] = 0.0;
        ipRVec[0] = -0.5*particle[pI].dia*uVec[0];
        ipRVec[1] = -0.5*particle[pI].dia*uVec[1];
        ipRVec[2] = -0.5*particle[pI].dia*uVec[2];

        jpRVec[0] = -ipRVec[0];
        jpRVec[1] = -ipRVec[1];
        jpRVec[2] = -ipRVec[2];
        printf("U VEC%lf\n",ipRVec[1]);
    }
    double gForce = particle[pI].mass; //gravitational acceleration
    particle[pI].forceY = -gForce;
}


/*Find contact forces with all neighbour particles*/
void neighbourContactForce(int pI){

}

/* Convert DEM Lagrangian values to Euler values & then makes a copy to Fluent 
secondary phase
*/
void copyDEMInfo(){
    printf("copyDEMInfo %d, %d, %d\n",xDiv,yDiv,zDiv);
}


/* Update particle position*/
void move(){
    for (int i=0; i<np; i++){
        double dxDot = 0.5*particle[i].forceX*timeStep/particle[i].mass;
        double dyDot = 0.5*particle[i].forceY*timeStep/particle[i].mass;
        double dzDot = 0.5*particle[i].forceZ*timeStep/particle[i].mass;

        particle[i].vX += dxDot;
        particle[i].vY += dyDot;
        particle[i].vZ += dzDot;

        double dX = particle[i].vX*timeStep;
        double dY = particle[i].vY*timeStep;
        double dZ = particle[i].vZ*timeStep;

        particle[i].posX += dX;
        particle[i].posY += dY;
        particle[i].posZ += dZ;
    }
}


