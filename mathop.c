#include "common.h"

/* Return center distance of two particles
param:
p1 - particle 1
p2 - particle 2 */
real getCenterDist(Tracked_Particle *ip, Tracked_Particle *jp){
    real ipX = P_POS(ip)[0]*lengthFactor;
    real ipY = P_POS(ip)[1]*lengthFactor;
    real ipZ = P_POS(ip)[2]*lengthFactor;
    real jpX = P_POS(jp)[0]*lengthFactor;
    real jpY = P_POS(jp)[1]*lengthFactor;
    real jpZ = P_POS(jp)[2]*lengthFactor;

    return sqrt((ipX-jpX)*(ipX-jpX) + (ipY-jpY)*(ipY-jpY) + (ipZ-jpZ)*(ipZ-jpZ));
 }

/* Add two vecotrs
param:
v1 - input vector 1
v2 - input vector 2
vec - resultant vector */
void vecAdd(real *v1, real *v2, real *vec){
   vec[0] = v1[0]+v2[0];
   vec[1] = v1[1]+v2[1];
   vec[2] = v1[2]+v2[2];
}

/* Substract two vecotrs
param:
v1 - input vector 1
v2 - input vector 2
vec - resultant vector */
void vecSub(real *v1, real *v2, real *vec){
   vec[0] = v1[0]-v2[0];
   vec[1] = v1[1]-v2[1];
   vec[2] = v1[2]-v2[2];
}

/* Find vector cross product
param:
v1 - input vector 1
v2 - input vector 2
vec - resultant vector */
void crossProd(real *v1, real *v2, real *vec){
    real temp1 = v1[1]*v2[2];
    real temp2 = v2[1]*v1[2];
    vec[0] = temp1 - temp2;

    temp1 = v2[0]*v1[2];
    temp2 = v1[0]*v2[2];
    vec[1] = temp1 - temp2;

    temp1 = v1[0]*v2[1];
    temp2 = v2[0]*v1[1];
    vec[2] = temp1 - temp2;   
}

/*
Dot product of two vectors
param:
v1 - vector 1
v2 - vector 2
*/
real dotProduct(real *v1, real *v2){
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

/*
Find normal unit vector to the surface defined by vector v1v2 and v2v3
*/
void getUnitVector(real *v1, real *v2, real *v3, real *uVec){
    real *v1v2 = allocateDoubleArray(DIM);
    real *v1v3 = allocateDoubleArray(DIM);

    vecSub(v1, v2, v1v2);
    vecSub(v1, v3, v1v3);
    crossProd(v1v2, v1v3, uVec);
    unitVec(uVec, uVec);

    free(v1v2);
    free(v1v3);
}

/* Find unit vector
param: 
v - input vector
vec - unit vector */
void unitVec(real *v, real *vec){
    if(vecMag(v) > 0){
        real temp1 = v[0]/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        real temp2 = v[1]/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        real temp3 = v[2]/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
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

/* Find magnitude of a vector*/
real vecMag(real *vec){
    return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
}

/* Multiply input vector by a scaler
param: 
scl - scaler
vec - vector to be multiplied*/
void sclMult(real scl, real *vec){
    vec[0] = scl*vec[0];
    vec[1] = scl*vec[1];
    vec[2] = scl*vec[2];
}

/* Multiply input vector by a scaler and assign to a vector
param: 
scl - scaler
inVec - input vector to be multiplied
outVec - reusltant output vector*/
void sclVecMult(real scl, real *inVec, real *outVec){
    outVec[0] = scl*inVec[0];
    outVec[1] = scl*inVec[1];
    outVec[2] = scl*inVec[2];
}


/* Returns the project vector on the plane defined by normal unit vector
param:
v - input vector
n - unit vector
vec - resultant project vector
type - either 0 or 1, 0 gives project vector, 1 gives project vector scaled by input vector */
void projVec(real *v, real *n, real *vec, int type){
    real *tV1 = malloc(DIM*sizeof(real));
    real *tV2 = malloc(DIM*sizeof(real));
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
        real *tV3 = malloc(DIM*sizeof(real));
        crossProd(n, v, tV1);
        tV2[0] = -n[0];
        tV2[1] = -n[1];
        tV2[2] = -n[2];
        crossProd(tV2, tV1, tV3);
        real temp = sqrt((tV2[0]*tV2[0]+tV2[1]*tV2[1]+tV2[2]*tV2[2])/
                            (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));  
        vec[0] = tV3[0]*temp;
        vec[1] = tV3[1]*temp;
        vec[2] = tV3[2]*temp;      
        free(tV3);
    }
    free(tV1);
    free(tV2);
}



