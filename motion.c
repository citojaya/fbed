#include "common.h"

/* 
DEM simulation
At every DPM time step this method is called by DEFINE_DPM_SOURCE
*/
void demLoop(Tracked_Particle *p){
    long int ip = p->part_id;

    updateNeighbourList(p);
    demPart[ip].force[0] = 0.0;
    demPart[ip].force[1] = -P_MASS(p)*massFactor; //gravitational force
    demPart[ip].force[2] = 0.0;


    demPart[ip].currentTime += demPart[ip].dt;
    //printf("CURRENT TIME %lf\n",CURRENT_TIME);
    //printf("TIME STEP %lf\n",CURRENT_TIMESTEP);

    // if(time_count > 1){
    //     demSave();    
    //     time_count = 0;
    // }
    // time_count++;
    //addToBdBox();
    //updateNeighbourList(); //Update neighbour list
}

/* Assign graviational force*/
void assignGravity(){
    for(int i=0; i<np; i++){
        real gForce = demPart[i].mass; //gravitational acceleration
        demPart[i].force[1] = -gForce; // vertical force component
    }
}

/* Find contact force */
void updateForce(Tracked_Particle *p){
    demPart[p->part_id].force[0] = 0.0;
    demPart[p->part_id].force[1] = -P_MASS(p)*massFactor; //Gravitional force
    demPart[p->part_id].force[2] = 0.0;
}


/*Calculate overlap with contact surface
param:
parPos - particle center
dia - particle diameter
node1, node2, node3 - contact surface nodes
uVec - unit vector normal to surface
*/
real getOverlap(real *parPos, real dia, real *n1, real *n2, real *n3, real *uVec){ 
    real *v1 = allocateDoubleArray(DIM); //vector running from node1 to particles center
    real *v2 = allocateDoubleArray(DIM); //vector running from node1 to particles center
    real *ppDash = allocateDoubleArray(DIM); //vector running from particle center to projection
    vecSub(parPos,n1,v1);  
    projVec(v1, uVec, v2, 0);
    vecSub(v1,v2,ppDash);
    uVec[0] = ppDash[0];
    uVec[1] = ppDash[1];
    uVec[2] = ppDash[2];

    unitVec(uVec,uVec);
    real overlap = vecMag(ppDash)-0.5*dia;
    free(v1);
    free(v2);
    free(ppDash);
    return overlap;
}


/*Find contact forces with all neighbour particles*/
void neighbourContactForce(Tracked_Particle *pI){
    //Loop through all neighbour particles
    for(int i=0; i<demPart[pI->part_id].noOfNeigh; i++){
        //printf("Neighbours\n");
        Tracked_Particle *jp = demPart[pI->part_id].neigh[i];
        real gap = getCenterDist(pI,jp)-(P_DIAM(pI)+P_DIAM(jp))*lengthFactor*0.5;
        if(gap <0){
            partContactForce(pI,jp, -gap);
        }
    }
}
/*
Calculate particle-wall contact forces
param:
ip - ith particle
nrmDisp - overlap
uVec - unit vector normal to contact surface
*/
void surfaceContactForce(Tracked_Particle *p, real nrmDisp, real *uVec){
    real rStar = 0.5*P_DIAM(p)*lengthFactor;
    
    //sclVecMult(-0.5*demPart[ip].dia,uVec,ipRVec);
    
    //crossProd(demPart[ip].angVel,ipRVec,rotVel);
    //vecAdd(demPart[ip].vel,rotVel,ipCntPntVel);

    real *relVel = allocateDoubleArray(DIM);
    real *pVel = allocateDoubleArray(DIM);
    pVel[0] = P_VEL(p)[0]*velocityFactor;
    pVel[1] = P_VEL(p)[1]*velocityFactor;
    pVel[2] = P_VEL(p)[2]*velocityFactor;
    sclVecMult(1.0,pVel,relVel);
    
    real nrmVel = dotProduct(relVel,uVec);
    free(relVel);
    free(pVel);
    
    //sclVecMult(1.0,ipCntPntVel,cntPntVel);

    real *totalForce = allocateDoubleArray(DIM);
    real *momentum = allocateDoubleArray(DIM);

    real nrmCntForce = elasticMod*sqrt(rStar*nrmDisp)*nrmDisp;
    real nrmDampForce = -dmpn*elasticMod*sqrt(rStar*nrmDisp)*nrmVel;

    real *cntPntDisp = allocateDoubleArray(DIM);
    sclVecMult(timeStep,cntPntVel,cntPntDisp);
    real *tipCntPntDisp = allocateDoubleArray(DIM);
    projVec(cntPntDisp, uVec, tipCntPntDisp, 0);
    real slidingDisp = vecMag(tipCntPntDisp);

    real *disp = allocateDoubleArray(DIM);
    projVec(demPart[p->part_id].hisDisp, uVec, disp, 1);
    vecAdd(disp, tipCntPntDisp, disp);
    dsmax = dsmaxCff*nrmDisp;
    dd = vecMag(disp);
    fdt = sfc*nrmCntForce;

    real *fdtVec = allocateDoubleArray(DIM);
    if(dd < 1e-6){
        dd = 0.0;
    }
    if(dti < 1e-6){
        dti = 0.0;
    }
    if(dd >= dsmax){
        sclMult(dsmax/dd,disp);
        dti = vecMag(tipCntPntDisp);
        if(dti != 0){
            sclVecMult(fdt/dti,tipCntPntDisp,fdtVec);
        }
        else{
            sclMult(0.0,fdtVec);
        }
    }
    else{
        if(dd != 0.0){
            fdt = fdt*(1.0 - pow((1.0 - dd/dsmax),1.5));
            sclVecMult(-fdt/dd,disp,fdtVec);
        }
        else{
            sclMult(0.0,fdtVec);
        }
    }

    sclVecMult(1.0,disp, demPart[p->part_id].hisDisp);

    //sum of forces
    real nrmForce = (nrmCntForce + nrmDampForce);
    sclVecMult(nrmForce, uVec, totalForce);
    //vecAdd(totalForce, fdtVec, totalForce);
    //crossProd(ipRVec, totalForce, momentum);
    //real *rotMom = allocateDoubleArray(DIM);
    //sclVecMult(0.5*rf*demPart[ip].dia*nrmCntForce, demPart[ip].angVel, rotMom);
    //vecAdd(momentum, rotMom, momentum);

    //Update force and momentum on particle
    vecAdd(demPart[p->part_id].force, totalForce, demPart[p->part_id].force);
    //vecAdd(demPart[ip].momentum, momentum, demPart[ip].momentum);

    //free(rotMom);
    
    free(fdtVec);
    free(tipCntPntDisp);
    free(cntPntDisp);
    free(totalForce);
    free(momentum);
    free(disp);
    
}

/*
Calculate interparticle forces
param:
ip - ith particle
jp - neighbour particle
nrmDisp - overlap
*/
void partContactForce(Tracked_Particle *ip, Tracked_Particle *jp, real nrmDisp){
    /*real rStar = 0.5*demPart[ip].dia*demPart[jp].dia/(demPart[ip].dia+demPart[jp].dia);

    vecSub(demPart[ip].pos, demPart[jp].pos, ijVec);
    unitVec(ijVec,uVec);
    sclVecMult(-0.5*demPart[ip].dia,uVec,ipRVec);
    sclVecMult(-0.5*demPart[jp].dia,uVec,jpRVec);

    crossProd(demPart[ip].angVel,ipRVec,rotVel);
    vecAdd(demPart[ip].vel,rotVel,ipCntPntVel);

    crossProd(demPart[jp].angVel,jpRVec,rotVel);
    vecAdd(demPart[jp].vel,rotVel,jpCntPntVel);

    real *relVel = allocateDoubleArray(DIM);
    vecSub(demPart[ip].vel,demPart[jp].vel,relVel);

    real nrmVel = dotProduct(relVel,uVec);
    vecSub(ipCntPntVel,jpCntPntVel,cntPntVel);
    free(relVel);
 
    real *totalForce = allocateDoubleArray(DIM);
    real *momentum = allocateDoubleArray(DIM);

    real nrmCntForce = elasticMod*sqrt(rStar*nrmDisp)*nrmDisp;
    real nrmDampForce = -dmpn*elasticMod*sqrt(rStar*nrmDisp)*nrmVel;

    real *cntPntDisp = allocateDoubleArray(DIM);
    sclVecMult(timeStep,cntPntVel,cntPntDisp);
    real *tipCntPntDisp = allocateDoubleArray(DIM);
    projVec(cntPntDisp, uVec, tipCntPntDisp, 0);
    real slidingDisp = vecMag(tipCntPntDisp);

    real *disp = allocateDoubleArray(DIM);
    projVec(demPart[ip].hisDisp, uVec, disp, 1);
    vecAdd(disp, tipCntPntDisp, disp);
    dsmax = dsmaxCff*nrmDisp;
    dd = vecMag(disp);
    fdt = sfc*nrmCntForce;

    real *fdtVec = allocateDoubleArray(DIM);
    if(dd < 1e-6){
        dd = 0.0;
    }
    if(dti < 1e-6){
        dti = 0.0;
    }
    if(dd >= dsmax){
        sclMult(dsmax/dd,disp);
        dti = vecMag(tipCntPntDisp);
        if(dti != 0){
            sclVecMult(fdt/dti,tipCntPntDisp,fdtVec);
        }
        else{
            sclMult(0.0,fdtVec);
        }
    }
    else{
        if(dd != 0.0){
            fdt = fdt*(1.0 - pow((1.0 - dd/dsmax),1.5));
            sclVecMult(-fdt/dd,disp,fdtVec);
        }
        else{
            sclMult(0.0,fdtVec);
        }
    }

    sclVecMult(1.0,disp, demPart[ip].hisDisp);

    //sum of forces
    real nrmForce = (nrmCntForce + nrmDampForce);
    sclVecMult(nrmForce, uVec, totalForce);
    vecAdd(totalForce, fdtVec, totalForce);
    crossProd(ipRVec, totalForce, momentum);
    real *rotMom = allocateDoubleArray(DIM);
    sclVecMult(0.5*rf*demPart[ip].dia*nrmCntForce, demPart[ip].angVel, rotMom);
    vecAdd(momentum, rotMom, momentum);

    //Update force and momentum on particle
    vecAdd(demPart[ip].force, totalForce, demPart[ip].force);
    vecAdd(demPart[ip].momentum, momentum, demPart[ip].momentum);

    free(rotMom);
    free(fdtVec);
    free(tipCntPntDisp);
    free(cntPntDisp);
    free(totalForce);
    free(momentum);
    free(disp);*/
}


/* Convert DEM Lagrangian values to Euler values & then makes a copy to Fluent 
secondary phase
*/
void copyDEMInfo(){
    printf("copyDEMInfo %d, %d, %d\n",xDiv,yDiv,zDiv);
}


/* Update particle position*/
void move(Tracked_Particle *p)
{
    long int ip = p->part_id;
    real parMass = P_MASS(p)*massFactor;
    // real dxDot = demPart[ip].force[0]*P_DT(p)*timeFactor/parMass;
    // real dyDot = demPart[ip].force[1]*P_DT(p)*timeFactor/parMass;
    // real dzDot = demPart[ip].force[2]*P_DT(p)*timeFactor/parMass;
    real dxDot = demPart[ip].force[0]*demPart[ip].dt/parMass;
    real dyDot = demPart[ip].force[1]*demPart[ip].dt/parMass;
    real dzDot = demPart[ip].force[2]*demPart[ip].dt/parMass;

    //if(demPart[ip].currentTime >= P_TIME(p)){
        P_VEL(p)[0] += dxDot/velocityFactor;
        P_VEL(p)[1] += dyDot/velocityFactor;
        P_VEL(p)[2] += dzDot/velocityFactor;

        real dX = P_VEL(p)[0]*velocityFactor*P_DT(p)*timeFactor;
        real dY = P_VEL(p)[1]*velocityFactor*P_DT(p)*timeFactor;
        real dZ = P_VEL(p)[2]*velocityFactor*P_DT(p)*timeFactor;

        P_POS(p)[0] += dX/lengthFactor;
        P_POS(p)[1] += dY/lengthFactor;
        P_POS(p)[2] += dZ/lengthFactor;
    //}
}


