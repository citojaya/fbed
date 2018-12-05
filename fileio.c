#include "common.h"

/*--- Find the reading position in the file--*/
void findRec(FILE *inFile, char* strDest){
	int nbytes = 256;
	char* strSrc;
	strSrc = (char *)malloc(nbytes+1);

	rewind(inFile);
	int n=strlen(strDest);
	while(!feof(inFile)){
		fgets(strSrc, 256, inFile);
		strSrc[n]='\0';
		if (strcmp(strDest, strSrc) == 0){
			break;
		}
	}

	if(strcmp(strDest, strSrc) != 0){
		//free(strSrc);
		//printf("Unable to find relevant info of: %s \n", strDest);
		exit(1);
	}
	free(strSrc);
}

/*---Read input data from a file ----*/
void readData(char *infile, int *np, double *dens, double *ymod, 
			double *pois, double *sfc, double *rec, double *dmpn, double *cyldia, double *dt){
	// input file reading
	char filename[20];
	strcpy(filename, infile);
	strcat(filename ,".in"); 
	FILE *InFile = fopen(filename, "rt");

	if (InFile == NULL){
		fprintf(stderr, "Can't open the parameter file %s! \n", filename);
		char c = getchar();
		exit(1);
	}
	// Log file
	char LogjName[20];
	strcpy(LogjName, infile);
	strcat(LogjName ,"_Log.dat"); 
	FILE *LogFile = fopen(LogjName, "a");

	// expand size for computing area
	double exComp = 0.0;

	// loading phase
	double DieDepth = 1.0;
	double UnDepth  = 0.0;
	double BdDepth  = 0.0;

	double parDia = 0.0;
	
	findRec(InFile, "PAR_NUMBER");
	fscanf(InFile, "%d",  np);

	findRec(InFile, "MATERIAL");
	fscanf(InFile, "%lf", dens);
	fscanf(InFile, "%lf", ymod);
	fscanf(InFile, "%lf", pois);
	fscanf(InFile, "%lf", sfc);
	fscanf(InFile, "%lf", dmpn);

	findRec(InFile, "CylinderBC");
	fscanf(InFile, "%lf", cyldia);
	fprintf(LogFile,"Domain size\n");
	fprintf(LogFile,"xDiv,yDiv,zDiv: %d,%d,%d\n :",xDiv,yDiv,zDiv);

	findRec(InFile, "SIMULATION");
	fscanf(InFile, "%lf", dt);	
	//fprintf(LogFile,"xmin,xmax %lf,%lf\n :",xmin,xmax);
	//fprintf(LogFile,"ymin,ymax %lf,%lf\n :",ymin,ymax);
	//fprintf(LogFile,"zmin,zmax %lf,%lf\n :",zmin,zmax);

    fclose(InFile);
	fclose(LogFile);
}

void diaInput(char *infile, struct Particle *par, int *np){
	char filename[20];
	strcpy(filename, infile);
	strcat(filename ,".in"); 
	FILE *pDiaFile = fopen(filename, "rt");

	if (pDiaFile == NULL){
		fprintf(stderr, "Can't open the parameter file %s! \n", filename);
		char c = getchar();
		exit(1);
	}
	int num = 0;

	double pDia, pX, pY, pZ;
    //printf("No of par %d\n",np);
	findRec(pDiaFile, "PARTICLE");
	for(int i=0; i<*np; i++){
		fscanf(pDiaFile, "%lf", &pDia);
		fscanf(pDiaFile, "%lf", &pX);
		fscanf(pDiaFile, "%lf", &pY);
		fscanf(pDiaFile, "%lf", &pZ);
		particle[i].dia = pDia;
		particle[i].posX = pX;
		particle[i].posY = pY;
		particle[i].posZ = pZ;
		//printf("D[%d]: %lf, %lf, %lf, %lf\n", i, particle[i].dia, particle[i].posX, particle[i].posY, particle[i].posZ);
	}
	fclose(pDiaFile);
}

void demSave(){
	//FILE *outfile; 
	char filename[20];
	sprintf(filename, "particle.dat");
	FILE *outfile = fopen(filename, "a");
	fprintf(outfile, "TIME = %lf\n",demTime/timeFactor);
	//fprintf(outfile, "VARIABLES = X   Y   Z   R   VX   VY   VZ\n");
	for(int i=0; i<np; i++){
		fprintf(outfile, "%11.5lf   %11.5lf   %11.5lf   %11.5f  %11.5lf   %11.5lf  %11.5lf\n",
	 		1e3*particle[i].posX/lengthFactor,1e3*particle[i].posY/lengthFactor,
			 1e3*particle[i].posZ/lengthFactor,1e3*particle[i].dia/lengthFactor,
			 particle[i].vX,particle[i].vY/velocityFactor,particle[i].vZ/velocityFactor); 
	}
	fclose(outfile);
	printf("SAVED\n");

}

// void writeTec(double *pPosX, double *parPosY, double *parPosZ){
// 	FILE *outfile; 
// 	char filename[20];
// 	sprintf(filename, "particle_info.dat");

// 	if (h_dem->Outs == 0)
// 	{
// 		outfile = fopen(filename, "wt");

// 		fprintf(outfile, "TITLE = \" PARTICLE INFORMATION \" \n");
// 		fprintf(outfile, "VARIABLES = X   Y   Z   R   VX   VY   VZ   W   F\n");
// 		fclose(outfile);
// 	}

// 	outfile = fopen(filename, "a");
// 	fprintf(outfile, "ZONE T= \" %12.6lf s \" \n", h_dem->ctime * Rdu->rtunit);
// 	for (int ip = 0; ip<TotalParticle; ip++)
// 	{
// 		fprintf(outfile, "%11.5lf   %11.5lf   %11.5lf   %11.5f  %11.5lf   %11.5lf   %11.5lf   %11.5lf   %11.5lf\n", 
// 		                  hPos[ip].x,   hPos[ip].y,   hPos[ip].z,  hRad[ip],
// 		                  hVel[ip].x,   hVel[ip].y,   hVel[ip].z, 
// 		                  length(hAngVel[ip]), length(hForce[ip]));
// 	}

// 	fclose(outfile);
// }

