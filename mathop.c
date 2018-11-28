#include "common.h"

double getCenterDist(int ip, int jp){
    return sqrt((parPosX[ip]-parPosX[jp])*(parPosX[ip]-parPosX[jp])+
        (parPosY[ip]-parPosY[jp])*(parPosY[ip]-parPosY[jp])+
        (parPosZ[ip]-parPosZ[jp])*(parPosZ[ip]-parPosZ[jp]));
 }

void vecAdd(double *v1, double *v2, double *vec){
   vec[0] = v1[0]+v2[0];
   vec[1] = v1[1]+v2[1];
   vec[2] = v1[2]+v2[2];
}

void vecSub(double *v1, double *v2, double *vec){
   vec[0] = v1[0]-v2[0];
   vec[1] = v1[1]-v2[1];
   vec[2] = v1[2]-v2[2];
}

void crossProd(double *v1, double *v2, double *vec){
    double temp1 = v1[1]*v2[2];
    double temp2 = v2[1]*v1[2];
    vec[0] = temp1 - temp2;

    temp1 = v2[0]*v1[2];
    temp2 = v1[0]*v2[2];
    vec[1] = temp1 - temp2;

    temp1 = v1[0]*v2[1];
    temp2 = v2[0]*v1[1];
    vec[2] = temp1 - temp2;   
}

void unitVec(double *v, double *vec){
    if(vecMag(v) > 0){
        double temp1 = v[0]/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        double temp2 = v[1]/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        double temp3 = v[2]/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        vec[0] = temp1;
        vec[1] = temp2;
        vec[2] = temp3;
    }
    else{
        vec[0] = 0.0;
        vec[1] = 0.0;
        vec[2] = 0.0;
    }
}

double vecMag(double *vec){
    return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
}

void sclVecMult(double scl, double *vec){
    vec[0] = scl*vec[0];
    vec[1] = scl*vec[1];
    vec[2] = scl*vec[2];
}

void projVec(double *v, double *n, double *vec, int type){
    double *tV1 = malloc(dim*sizeof(double));
    double *tV2 = malloc(dim*sizeof(double));
    if(type == 0){
        crossProd(n,v, tV1);
        tV2[0] = -n[0];
        tV2[1] = -n[1];
        tV2[2] = -n[2];
        //printf("tV1 %lf,%lf,%lf\n ",tV1[0],tV1[1],tV1[2]);
        //printf("tV2 %lf,%lf,%lf\n ",tV2[0],tV2[1],tV2[2]);
        crossProd(tV2, tV1, vec);
        //printf("Type 0\n");
    }
    else{
        double *tV3 = malloc(dim*sizeof(double));
        crossProd(n, v, tV1);
        tV2[0] = -n[0];
        tV2[1] = -n[1];
        tV2[2] = -n[2];
        crossProd(tV2, tV1, tV3);
        double temp = sqrt((tV2[0]*tV2[0]+tV2[1]*tV2[1]+tV2[2]*tV2[2])/
                            (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));  
        vec[0] = tV3[0]*temp;
        vec[1] = tV3[1]*temp;
        vec[2] = tV3[2]*temp;      
        free(tV3);
    }
    free(tV1);
    free(tV2);
}

