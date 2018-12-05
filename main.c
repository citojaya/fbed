#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"

//#include "double_math.h"

/* Program begining*/
int main(void){
    run();
}

/* Run the program*/
void run(){
    //***** Test code ***********
    double xmin = 0.0;
    double xmax = 10.0;
    double ymin = 0.0;
    double ymax = 88.0;
    double zmin = 0.0;
    double zmax = 200.0;

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

    //Initialise DEM parameters, read particle information, assign arrays
    demInit(xD, yD, zD, dx, dy, dz);
    

    // for (int i=0; i<xDiv*yDiv*zDiv; i++){
    //     bdBox[i].fluidVelX = 0.0;
    //     bdBox[i].fluidVelY = 10.0;
    //     bdBox[i].fluidVelZ = 0.0;
    //     bdBox[i].fluidVolF = 0.5;
    // }

    for()
    demLoop();
    demSave();
    


 /******* Testing code for neighbourlist, currently not in use *******/
 /*   
    // Initialize neighbourlist array in X direction for the first time
    initialize(sortedList, sortedParIndex, cellSE, np, parPosX, parDia);

    printf("Before sorting\n");
    for (int i=0; i<np*2; i++){
        printf("%lf, %d, %d\n", sortedList[i], sortedParIndex[i], cellSE[i]);
    }

    //Neighbourlist sorting for the first time
    insertionSort(sortedList, np*2, sortedParIndex, cellSE, 1);

    //Assign particle neighbours for the first time
    assignNeighbours(sortedList, sortedParIndex, cellSE, np*2);
*/
    
/*
    //Sort neighbourlist for each iteration
    insertionSort(sortedList, np*2, sortedParIndex, cellSE, 0);

    /****************************/
/*    printf("After sorting\n");
    for (int i=0; i<np*2; i++){
        printf("%lf, %d, %d\n", sortedList[i], sortedParIndex[i], cellSE[i]);
    }

    printf("Neighbours\n");
    for (int i=0; i<np; i++){
        printf("PARTICLE %d\n",i);
        printf("%d\n",parNoOfNb[i]);
        for (int j=0; j<parNoOfNb[i]; j++){
            printf("%d   ",parNb[i*nbSize+j]);
        }
        printf("\n");
    }
    */
/*********************************************************************************/

    // Delete dynamic memeory
    // free(parPosX);
    // free(parPosY);
    // free(parPosZ);
    free(sortedList);
    free(sortedParIndex);
    free(cellSE);
    //free(parIndex);
//    free(parDia);
    free(parNb);
    free(parNoOfNb);
    free(parCIndex);
    // free(parMass);
    // free(parInert);
    //free(parMass);
    free(uVec);
    free(ipRVec);
    free(jpRVec);
    free(bdBox);
    free(particle);
    printf("All good!\n");
    // free(nebListIndex);
}
