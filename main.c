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
    /*double *pForce = malloc(dim*sizeof(double));
    double *uVec = malloc(dim*sizeof(double));
    double *res = malloc(dim*sizeof(double));
    double *res2 = malloc(dim*sizeof(double));

    uVec[0] = 5.0;
    uVec[1] = 6.7;
    uVec[2] = 2.85;

    pForce[0] = 0.0;
	pForce[1] = 100.0;
	pForce[2] = 0.0;


    unitVec(uVec, res);
    //crossProd(uVec, pForce, res);
    printf("Results 1 %lf, %lf, %lf\n", res[0],res[1],res[2]);
    //projVec(pForce, uVec, res2, 0);
    //printf("Results 2 %lf, %lf, %lf\n", res2[0],res2[1],res2[2]);

    
    free(pForce);
    free(uVec);
    free(res);
    free(res2);*/
    //*************************
    
    //Initialise DEM parameters, read particle information, assign arrays


    demInit();
    //Build cell map for DEM particles and FLuent cells
    //buildDEMCFDCellMap();


    /* DEM iteration*/
    //demLoop();


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
    free(parMass);
    free(parInert);
    printf("All good!\n");
    // free(nebListIndex);
}
