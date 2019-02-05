#include "common.h"

/* 
DEM simulation
At every DPM time step this method is called by DEFINE_DPM_SOURCE
*/
// void demLoop(Tracked_Particle *p){
//     long int ip = p->part_id;

//     updateNeighbourList(p);
//     demPart[ip].currentTime += demPart[ip].dt;
// }

/* Assign graviational force*/
void assignGravity(){
    for(int i=0; i<np; i++){
        real gForce = demPart[i].mass; //gravitational acceleration
        demPart[i].force[1] = -gForce; // vertical force component
    }
}

/* Find contact force */
// void updateForce(Tracked_Particle *p){
//     demPart[p->part_id].force[0] = 0.0;
//     demPart[p->part_id].force[1] = -P_MASS(p)*massFactor; //Gravitional force
//     demPart[p->part_id].force[2] = 0.0;
// }


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
void neighbourContactForce(int pI){
    //Loop through all neighbour particles
    for(int i=0; i<demPart[pI].noOfNeigh; i++){
        //printf("Neighbours\n");
        int jp = demPart[pI].neigh[i];
        real gap = getCenterDist(pI,jp)-(demPart[pI].dia+demPart[jp].dia)*0.5;
        if(gap < 0.0){
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
void surfaceContactForce(int p, real nrmDisp, real *uVec){
    //real rStar = 0.5*P_DIAM(p)*lengthFactor;
    real rStar = 0.5*demPart[p].dia;
    sclVecMult(-0.5*demPart[p].dia,uVec,ipRVec);
    
    crossProd(demPart[p].angVel,ipRVec,rotVel);
    vecAdd(demPart[p].vel,rotVel,ipCntPntVel);

    real *relVel = allocateDoubleArray(DIM);
    // real *pVel = allocateDoubleArray(DIM);
    // pVel[0] = P_VEL(p)[0]*velocityFactor;
    // pVel[1] = P_VEL(p)[1]*velocityFactor;
    // pVel[2] = P_VEL(p)[2]*velocityFactor;

    
    //sclVecMult(1.0,pVel,relVel);
    sclVecMult(1.0,demPart[p].vel,relVel);
    
    real nrmVel = dotProduct(relVel,uVec);
    free(relVel);
    //free(pVel);
    
    sclVecMult(1.0,ipCntPntVel,cntPntVel);

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
    projVec(demPart[p].hisDisp, uVec, disp, 1);
    vecAdd(disp, tipCntPntDisp, disp);
    dsmax = dsmaxCff*nrmDisp;
    dd = vecMag(disp);
    fdt = sfc*nrmCntForce;

    // writeLogNum("logfile2.log"," disp X ",disp[0]);
    // writeLogNum("logfile2.log"," disp Y ",disp[1]);
    // writeLogNum("logfile2.log"," disp Z ",disp[2]);

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

    sclVecMult(1.0,disp, demPart[p].hisDisp);

    //sum of forces
    real nrmForce = (nrmCntForce + nrmDampForce);
    sclVecMult(nrmForce, uVec, totalForce);

    // writeLogNum("logfile2.log"," FDT X ",fdtVec[0]);
    // writeLogNum("logfile2.log"," FDT Y ",fdtVec[1]);
    // writeLogNum("logfile2.log"," FDT Z ",fdtVec[2]);
    vecAdd(totalForce, fdtVec, totalForce);
    crossProd(ipRVec, totalForce, momentum);
    real *rotMom = allocateDoubleArray(DIM);
    sclVecMult(0.5*rf*demPart[p].dia*nrmCntForce, demPart[p].angVel, rotMom);
    vecAdd(momentum, rotMom, momentum);

    //Update force and momentum on particle
    vecAdd(demPart[p].force, totalForce, demPart[p].force);
    vecAdd(demPart[p].momentum, momentum, demPart[p].momentum);

    free(rotMom);
    
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
void partContactForce(int ip, int jp, real nrmDisp){
    real rStar = 0.5*demPart[ip].dia*demPart[jp].dia/(demPart[ip].dia+demPart[jp].dia);
    vecSub(demPart[ip].pos, demPart[jp].pos, ijVec);
    unitVec(ijVec,uVec);

    //writeLog("logfile2.log","GAP ",nrmDisp);
    //real *relVel = allocateDoubleArray(DIM);
    //real *pVel = allocateDoubleArray(DIM);

    sclVecMult(-0.5*demPart[ip].dia,uVec,ipRVec);
    sclVecMult(0.5*demPart[jp].dia,uVec,jpRVec);

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
    //writeLog("logfile2.log","nrmForce ",nrmForce);
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
    free(disp);
}


/* Convert DEM Lagrangian values to Euler values & then makes a copy to Fluent 
secondary phase
*/
void copyDEMInfo(){
    printf("copyDEMInfo %d, %d, %d\n",xDiv,yDiv,zDiv);
}


/* Update particle position*/
void forceCalculation(Particle *p)
{
    real dt = timeStep;
    //writeLog("logfile2.log","MOVE DEM PART time ",P_TIME(p));

    demPart[p->part_id].force[0] = 0.0;
    demPart[p->part_id].force[1] = -demPart[p->part_id].mass; //gravitational force
    demPart[p->part_id].force[2] = 0.0;

    //Add Bouynacy
    Thread *tc = P_CELL_THREAD(p);
    cell_t c = P_CELL(p);
    //demPart[p->part_id].force[1] += partVol(p->part_id)*C_R(c,tc)*densityFactor;

    //Find contact force with wall

    //Domain *d = Get_Domain(1);
    //Thread *t = Lookup_Thread(d,0);

    //Thread **tph = THREAD_SUB_THREADS(t);
    //writeLog("logfile2.log","CELL ",C_VOF(P_CELL(p),t));
    int i;
    face_t f;
    // c_face_loop(c,tc,i)//Loop over faces
    // {
    //     int t_id = THREAD_ID(C_FACE_THREAD(c,tc,i));
    //     //writeLog("logfile2.log","face thread ID ",t_id);
    //     if(t_id == 5)
    //     {
    //         writeLog("logfile2.log","face thread ID ",t_id);
    //     f = C_FACE(c,tc,i);
    //     Thread *tf = C_FACE_THREAD(c,tc,i);

    //     //int face_t_id = THREAD_ID(C_FACE_THREAD(cc,t,i));
    //     int n;
    //     int j=0;
    //     f_node_loop(f,tf,n)//Loop over face nodes
    //     {
    //       Node *node = F_NODE(f,tf,n);
    //       if(j==0)
    //       {
    //         demPart[p->part_id].faceNode1[0] = NODE_X(node)*lengthFactor;
    //         demPart[p->part_id].faceNode1[1] = NODE_Y(node)*lengthFactor;
    //         demPart[p->part_id].faceNode1[2] = NODE_Z(node)*lengthFactor;
    //       }
    //       if(j==1)
    //       {
    //         demPart[p->part_id].faceNode2[0] = NODE_X(node)*lengthFactor;
    //         demPart[p->part_id].faceNode2[1] = NODE_Y(node)*lengthFactor;
    //         demPart[p->part_id].faceNode2[2] = NODE_Z(node)*lengthFactor;
    //       }
    //       if(j==2)
    //       {
    //         demPart[p->part_id].faceNode3[0] = NODE_X(node)*lengthFactor; 
    //         demPart[p->part_id].faceNode3[1] = NODE_Y(node)*lengthFactor;
    //         demPart[p->part_id].faceNode3[2] = NODE_Z(node)*lengthFactor;
    //       }
    //       j++;
    //     }//end of loop over face nodes
    //     real gap = demPart[p->part_id].pos[1] - demPart[p->part_id].dia*0.5;
    //     real *uVec = allocateDoubleArray(DIM);
    //     getUnitVector(demPart[p->part_id].faceNode1,demPart[p->part_id].faceNode2
    //                 ,demPart[p->part_id].faceNode3,uVec);
    //     //real *parPos = allocateDoubleArray(DIM);
    //     //real gap = getOverlap(demPart[p->part_id].pos,demPart[p->part_id].dia, demPart[p->part_id].faceNode1,  
    //     //                        demPart[p->part_id].faceNode2, demPart[p->part_id].faceNode3, uVec);
     
    //     //writeLog("logfile2.log","PAR TIME ",P_TIME(p));
    //     if(gap < 0) //If contact exists calculate contact force
    //     {
    //         surfaceContactForce(p->part_id, -gap, uVec);
    //     } 
    //     free(uVec);
    //     //free(parPos);
    //     }//end of if(t_id == 5)      
    // }//end of loop over faces


    //Fluid drag force calculation
    real velFX = C_U(c,tc)*velocityFactor;  
    real velFY = C_V(c,tc)*velocityFactor;  
    real velFZ = C_W(c,tc)*velocityFactor;    
    real pGX = C_P_G(c,tc)[0]*pressureFactor/lengthFactor;
    real pGY = C_P_G(c,tc)[1]*pressureFactor/lengthFactor;
    real pGZ = C_P_G(c,tc)[2]*pressureFactor/lengthFactor;
    real velPX = demPart[p->part_id].vel[0];
    real velPY = demPart[p->part_id].vel[1];
    real velPZ = demPart[p->part_id].vel[2];
    real density = C_R(c,tc)*densityFactor;
    real visc = VISCOSITY*massFactor/(lengthFactor*timeFactor);

    real instPor = 1.0-solidFraction(p->part_id);//Solid fraction;
    //writeLogNum("logfile2.log", "PORO  ",instPor);
    //--- Test code ---------------------------------------------------------
    // real Re = instPor*demPart[p->part_id].dia*sqrt((velFX-velPX)*(velFX-velPX)  + (velFY-velPY)*(velFY-velPY) 
    //                 + (velFZ-velPZ)*(velFZ-velPZ))*density/visc;
    
    // real coeff, dCoeff,modPor;

    //-- Standard equations -------------------------------------------------
	// if (Re < 1e-6)//if Reynolds number is too small following parameters cannot be defined
	// {
	// 	coeff = 0.0;
	// 	dCoeff = 0.0;
	// 	modPor = 0.0;
	// }
	// else 
	// {
	// 	coeff = 3.7-0.65*exp(-(1.5-log10(Re))*(1.5-log10(Re))*0.5);
	// 	dCoeff = (0.63+4.8/sqrt(Re))*(0.63+4.8/sqrt(Re));
	// 	modPor = pow(instPor,(-coeff));
	// }

    // real pfFX = modPor*dCoeff*density*instPor*(velFX-velPX)*fabs(velFX-velPX)*PI*(pow(demPart[p->part_id].dia*0.5,2));        
    // real pfFY = modPor*dCoeff*density*instPor*(velFY-velPY)*fabs(velFY-velPY)*PI*(pow(demPart[p->part_id].dia*0.5,2));        
    // real pfFZ = modPor*dCoeff*density*instPor*(velFZ-velPZ)*fabs(velFZ-velPZ)*PI*(pow(demPart[p->part_id].dia*0.5,2));        

	//------------------------------------------------------------------------

    real relVelMag = sqrt((velFX-velPX)*(velFX-velPX)+(velFY-velPY)*(velFY-velPY)+(velFZ-velPZ)*(velFZ-velPZ));
    real Re = instPor*demPart[p->part_id].dia*relVelMag*density/visc;
    
    real dragCoeff, beta;
    if(Re <= 1000){
        dragCoeff = 0.5*(24.0/Re)*(1.0+0.15*pow(Re,0.687));
    }
    else if(Re > 1000){
        dragCoeff = 0.43;
    }
    
    //writeLogNum("logfile2.log", "DRAG ", dragCoeff);
    if(instPor <= 0.8){
        real b1 = (1.0-instPor)/(demPart[p->part_id].dia*pow(instPor, 2));
        real b2 = 150.0*(1.0-instPor)*visc/demPart[p->part_id].dia;
        real b3 = 1.75*density*instPor*relVelMag;
        beta = b1*(b2 + b3);
    }
    else if(instPor > 0.8){
        beta = (3.0/4.0)*dragCoeff*density*relVelMag*(1.0-instPor)*pow(instPor,-2.7)/demPart[p->part_id].dia;
    }

    real pfFX = (velFX-velPX)*partVol(p->part_id)*beta/(1.0-instPor);
    real pfFY = (velFY-velPY)*partVol(p->part_id)*beta/(1.0-instPor);
    real pfFZ = (velFZ-velPZ)*partVol(p->part_id)*beta/(1.0-instPor);

    real pGFX = -pGX*PI*pow(demPart[p->part_id].dia,3)/6.0;
    real pGFY = -pGY*PI*pow(demPart[p->part_id].dia,3)/6.0;
    real pGFZ = -pGZ*PI*pow(demPart[p->part_id].dia,3)/6.0;

    // writeLogNum("logfile2.log","pGFX",pGFX);
    // writeLogNum("logfile2.log","pGFY",pGFY);
    // writeLogNum("logfile2.log","pGFZ",pGFZ);

    // writeLogLine("logfile2.log","--------------------");


    // writeLogNum("logfile2.log","pfFX",pfFX);
    // writeLogNum("logfile2.log","pfFY",pfFY);
    // writeLogNum("logfile2.log","pfFZ",pfFZ);

    // writeLogLine("logfile2.log","--------------------");

    // End of fluid drag force calculation
    


    real xMax = 0.010*lengthFactor;
    real zMax = 0.088*lengthFactor;
    real yMax = 0.2*lengthFactor;

    //Contact with x=0 wall
    real gap = demPart[p->part_id].pos[0] - demPart[p->part_id].dia*0.5;
    uVec[0] = 1.0;
    uVec[1] = 0.0;
    uVec[2] = 0.0;  
    if(gap < 0) //If contact exists calculate contact force
    {
       surfaceContactForce(p->part_id, -gap, uVec);
    }  

    //Contact with x=10mm wall
    gap = xMax - demPart[p->part_id].pos[0] - demPart[p->part_id].dia*0.5;
    uVec[0] = -1.0;
    uVec[1] = 0.0;
    uVec[2] = 0.0;  
    if(gap < 0) //If contact exists calculate contact force
    {
       surfaceContactForce(p->part_id, -gap, uVec);
    }  

    //Contact with z=0mm wall
    gap = demPart[p->part_id].pos[2] - demPart[p->part_id].dia*0.5;
    uVec[0] = 0.0;
    uVec[1] = 0.0;
    uVec[2] = 1.0;  
    if(gap < 0) //If contact exists calculate contact force
    {
       surfaceContactForce(p->part_id, -gap, uVec);
    }  

    //Contact with z=88mm wall
    gap = zMax - demPart[p->part_id].pos[2] - demPart[p->part_id].dia*0.5;
    uVec[0] = 0.0;
    uVec[1] = 0.0;
    uVec[2] = -1.0;  
    if(gap < 0) //If contact exists calculate contact force
    {
       //writeLogNum("logfile2.log","GAP ",gap/lengthFactor);
       surfaceContactForce(p->part_id, -gap, uVec);
    }  

    // Contact with bottom
    gap = demPart[p->part_id].pos[1] - demPart[p->part_id].dia*0.5;
    uVec[0] = 0.0;
    uVec[1] = 1.0;
    uVec[2] = 0.0;  

    if(gap < 0) //If contact exists calculate contact force
    {
        //demPart[p->part_id].force[1] += -demPart[p->part_id].mass;
        surfaceContactForce(p->part_id, -gap, uVec);
    }  

    // Contact with Top
    gap = yMax - (demPart[p->part_id].pos[1] + demPart[p->part_id].dia*0.5);
    uVec[0] = 0.0;
    uVec[1] = -1.0;
    uVec[2] = 0.0;  

    if(gap < 0) //If contact exists calculate contact force
    {
        surfaceContactForce(p->part_id, -gap, uVec);
    }  

    //Find particle-particle contact force
    neighbourContactForce(p->part_id);

    // writeLogNum("logfile2.log","pfFY ",pfFY);
    // writeLogNum("logfile2.log","pGFY ",pGFY);
    demPart[p->part_id].force[0]  += pfFX + pGFX;
    demPart[p->part_id].force[1]  += pfFY + pGFY + partVol(p->part_id)*density;
    demPart[p->part_id].force[2]  += pfFZ + pGFZ;


}

void updatePosition(Particle *p){
    
    real dxDot = demPart[p->part_id].force[0]*timeStep/demPart[p->part_id].mass;
    real dyDot = demPart[p->part_id].force[1]*timeStep/demPart[p->part_id].mass;
    real dzDot = demPart[p->part_id].force[2]*timeStep/demPart[p->part_id].mass;
    
    demPart[p->part_id].vel[0] += dxDot;
    demPart[p->part_id].vel[1] += dyDot;
    demPart[p->part_id].vel[2] += dzDot;

    demPart[p->part_id].pos[0] += demPart[p->part_id].vel[0]*timeStep;
    demPart[p->part_id].pos[1] += demPart[p->part_id].vel[1]*timeStep;
    demPart[p->part_id].pos[2] += demPart[p->part_id].vel[2]*timeStep;
   
    demPart[p->part_id].currentTime += timeStep;

    demPart[p->part_id].displacement += sqrt(pow(dxDot,2)+pow(dyDot,2)+pow(dzDot,2));

}

/* Update particle position*/
void moveDPM(Tracked_Particle *p)
{
    //writeLog("logfile2.log","TIME STEP ",(real)demPart[p].dt/timeFactor);
    //writeLog("logfile2.log","MASS ",(real)demPart[p].mass/massFactor);
    //writeLog("logfile2.log","TIME STEP ",(real)demPart[p].dt/timeFactor);
    //if(P_TIME(p) != demPart[p->part_id].currentTime){
    int dt = P_TIME(p) - demPart[p->part_id].currentTime;
    //dt = P_TIME(p)*timeFactor;
    real dxDot = demPart[p->part_id].force[0]*dt*timeFactor/(P_MASS(p)*massFactor);
    real dyDot = demPart[p->part_id].force[1]*dt*timeFactor/(P_MASS(p)*massFactor);
    real dzDot = demPart[p->part_id].force[2]*dt*timeFactor/(P_MASS(p)*massFactor);
    
    P_VEL(p)[0] += dxDot/velocityFactor;
    P_VEL(p)[1] += dyDot/velocityFactor;
    P_VEL(p)[2] += dzDot/velocityFactor;
    demPart[p->part_id].currentTime = P_TIME(p);
    //}
}



