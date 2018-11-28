#include "common.h"

/* DEM simulation
*/
void demLoop(){
    printf("demLoop\n");

    updateForce(); //Calculate contact forces 
    move();// update particle position
    //updateNeighList(); //Update neighbour list
}

/* Convert DEM Lagrangian values to Euler values & then makes a copy to Fluent 
secondary phase
*/
void copyDEMInfo(){
    printf("copyDEMInfo\n");
}

/* Find contact force */
void updateForce(){
    for (int i=0; i<np; i++){
        neighbourContactForce(i);
        boundaryContactForce(i);
    }
}

/*Find contact forces with all boundaries*/
void boundaryContactForce(int parI){
    double rStar = 0.5*parDia[parI];
	double pXY = sqrt(pow(parPosX[parI],2)+pow(parPosY[parI],2));
	double gap = 0.5*cyldia - pXY - rStar;
	double sfCoff = sfc;

    if(gap < 0){
        uVec[0] = -parPosX[parI]/pXY;
        uVec[1] = -parPosY[parI]/pXY;
        uVec[2] = 0.0;
        ipRVec[0] = -0.5*parDia[parI]*uVec[0];
        ipRVec[1] = -0.5*parDia[parI]*uVec[1];
        ipRVec[2] = -0.5*parDia[parI]*uVec[2];

        jpRVec[0] = -0.5*cyldia*uVec[0];
        jpRVec[1] = -0.5*cyldia*uVec[1];
        jpRVec[2] = -0.5*cyldia*uVec[2];
    }

}

/*Find contact forces with all neighbour particles*/
void neighbourContactForce(int parI){

}

/* Update particle position*/
void move(){

}

/* Update particle position*/
/*void update(double *pX, int np){
    double tempX = pX[1];
    pX[1] = pX[0]; 
    pX[0] = tempX;
}*/

