#include "common.h"



/*
Add a new neighbour to the neighbour list
*/
void addNeighbour(Tracked_Particle  *ip, Tracked_Particle *jp){
    demPart[ip->part_id].neigh[demPart[ip->part_id].noOfNeigh] = jp;
    demPart[ip->part_id].noOfNeigh++;
}

/*
Delete neighbour
*/
void deleteNeighbour(Tracked_Particle  *ip, Tracked_Particle *jp){
    for (int i=0; i<demPart[ip->part_id].noOfNeigh; i++){
        Tracked_Particle *neigh = demPart[ip->part_id].neigh[i];
        if(neigh->part_id == jp->part_id){
            demPart[ip->part_id].neigh[i] = demPart[ip->part_id].neigh[demPart[i].noOfNeigh-1];
            demPart[ip->part_id].noOfNeigh--;
            break;
        }
    }

}

/*
Add particles to bounding box
*/
void addToBdBox(){
    printf("************************\n");
    //Reset to zero
    for(int i=0; i<xDiv*yDiv*zDiv; i++){
        bdBox[i].noOfParticles = 0;
    }

    Injection *I;
    Injection *Ilist = Get_dpm_injections();
    //Reset neighbourlsit to zero
    loop(I,Ilist)
    {
        Particle *p;
        loop(p,I->p)
        {
    //Assign to bDBox cells and update neighbour list
    
        int iIndex = ceil((P_POS(p)[0]*lengthFactor-xmin)/domainDx);
        int jIndex = ceil((P_POS(p)[1]*lengthFactor-ymin)/domainDy);
        int kIndex = ceil((P_POS(p)[2]*lengthFactor-zmin)/domainDz);
        int cellIndex = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;
        
        printf("CELL INDEX %d",cellIndex);
        //writeLog("logfile2.log","ADD TO BD BOX",(real)cellIndex);
        //insert particle to cell
        //insertToBdBox(p,cellIndex);
        }
    }
        //printf("Cell index %d, %d, %d\n",iIndex,jIndex,kIndex);
    
}


/*
Update neighbour list for a given particle
For a given particle scan through all neighbour cells and fetch particles within
the neighbour region
param:
pI - particle index
*/
void updateNeighbourList(Tracked_Particle *ip){
    // Injection *I;
    // Injection *Ilist = Get_dpm_injections();
    // //Reset neighbourlsit to zero
    // loop(I,Ilist)
    // {
    //     Particle *p;
    //     loop(p,I->p)
    //     {
            demPart[ip->part_id].noOfNeigh = 0;
            int iIndex = ceil((P_POS(ip)[0]*lengthFactor-xmin)/domainDx);
            int jIndex = ceil((P_POS(ip)[1]*lengthFactor-ymin)/domainDy);
            int kIndex = ceil((P_POS(ip)[2]*lengthFactor-zmin)/domainDz);
            //printf("i,j,k index %d, %d, %d\n",iIndex, jIndex, kIndex);
            for(int r=kIndex-1; r<kIndex+2; r++){
                for(int q=jIndex-1; q<jIndex+2; q++){
                    for(int p=kIndex-1; p<kIndex+2; p++){
                        int neighCellIndex = p + q*xDiv + r*xDiv*yDiv;
                        for(int j=0; j<bdBox[neighCellIndex].noOfParticles; j++){
                            Tracked_Particle *jp = bdBox[neighCellIndex].parts[j];
                            //printf("PAR DIA %lf,%lf\n",P_POS(demPart[ip].cfdp)[0],P_POS(demPart[jp].cfdp)[0]);
                            // if(getCenterDist(p,jp)-0.5*(P_DIAM(p)+P_DIAM(jp))*lengthFactor < cutGap){
                            //     printf("Add neigh\n");
                            //     addNeighbour(p,jp);
                            // }
                        }
                    }
                }
            }
    //     }
    // }
}



