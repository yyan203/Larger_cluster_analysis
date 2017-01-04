
// COMPUTE Neighbor Debris Atom probibility according to distance cutoff (species dependent) & time duration 
// DO CLUSTER ANALYSIS on Debris Atoms according to distance cutoff (species dependent) & time duration 
//
// Input file format is lata4olivia type:
// ITEM: ATOMS id type x y z 
#include <cmath>
#include <algorithm>
#include <fstream>
#include <utility>  //std::pair, std::make_pair
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <cstdlib>
#include <map>
#include <vector>
#include <set>

#define ONELINEMAX 2000
#define MAXDUMP 1000
#define MaxAtoms 2000
#define MaxClusters 130000
#define EMPTY -100//for linkcell and other purpose
#define VERBOSE 1//for print debug 
#define TRUE 1
#define FALSE 0

#include "atom.h"
#include "event.h"
//using namespace std;

/////////////////////////////////
// Global variables defination //
/////////////////////////////////
// check if a file exist
bool fexists(const char *filename)
{
  std::ifstream ifile(filename);
    return ifile;
}
// Cluster Analysis
int main(int argc, char *argv[])
{
  Atom ***allatom;
  std::vector <Event> clusterevent;
  int ntotal;
  char oliviafilename[100], outputfilename[100],out2[100];
  FILE *ifp, *ofp, *ofp2;
  char templine[ONELINEMAX], str1[100], str2[100];
  //int count_frame_Atom[MAXDUMP]={0},count_frame_Correlation[MAXDUMP]={0};
  float xlo, xhi, ylo, yhi, zlo, zhi;
  int current_timestep;
  float dx, dy,dz,distance; 
  float ix, iy, iz;
  float cut;
  int id, type;
  float edgeY,maxF,minF,binl;
  char out[200];
  int i,j,k,nmols=-1;
  int n;
  int flag;
  int totalnmols=0;
  int frame;
  //int countAtom=0, countCorrelation=0;
  int fromframe, Step, toframe, extension;
   if(argc>=5)
    {
      sscanf(argv[1], "%s", oliviafilename);
      sscanf(argv[2], "%s", outputfilename);
      sscanf(argv[3], "%d", &fromframe);
      sscanf(argv[4], "%d", &Step);
      sscanf(argv[5], "%d", &toframe);
      sscanf(argv[6], "%d", &extension);
      sscanf(argv[7], "%f", &cut);
    }
    else
    {
      std::cout<<"Correct syntax: "<< argv[0]<< " lata4olivia outputfilename fromframe Step toframe Extension(consider atoms correlation extended to how many frames after current one) Cutoff(define a debris event)\n";
      std::cout<<"!! Output cluster info only contains clusters containing at least one atoms belonging to frame [fromframe to (toframe-extension)],(inclusive))\n";
      std::cout<<"!! Total Debris Atoms must less than 200000, otherwise please increase it in the C++ code by changing constant: MaxAtoms\n";
      std::cout<<"!! Different Atom Species Pairs have different distance cutoff,they are defined in function.h!!! \n";
      exit(0);
    }
//initialized container
    ntotal=(toframe-fromframe)/Step+1;
    allatom = new Atom **[ntotal];
    for(i=0;i<ntotal;i++) 
    {
    allatom[i] = new Atom *[MaxAtoms];
    for(j=0;j<MaxAtoms;j++) 
    allatom[i][j] = new Atom();
    allatom[i][i]->exist=0;
    }
    printf("Initialization of atom has been finished! \n");
    sprintf(out,"%s",oliviafilename);

     if(fexists(out))
     {
     printf("file:%s exist\n",out);
     flag=0;
     ifp = fopen(out,"r");
     //if(ifp==NULL) { std::cout << oliviafilename << " not exist! \n"; exit(0);} 
     printf("It is done here 0\n");
     //fgets(templine, ONELINEMAX, ifp);flag++;// Olivia First Line 
	 // if(atom==NULL) printf("Error in allocating memory!\n");
	 // if(newatom==NULL) printf("Error in allocating memory!\n");
      n=0;

      std::cout<< "start reading in the input file..."<< '\n'<< std::endl;
      frame = 0;
while(!feof(ifp))
      { 
      fgets(templine, ONELINEMAX, ifp); // TIMESTEP line
      sscanf(templine, "%s %s", str1, str2);
  if(strcmp(str1, "ITEM:")==0) {
      fgets(templine, ONELINEMAX, ifp); //actual timesteps
      sscanf(templine, "%ld", &current_timestep);

      //std::cout<< "reading frame "<<current_timestep<<"\n";
      fgets(templine, ONELINEMAX, ifp); // ITEM: NUMBER OF ATOMS
      fgets(templine, ONELINEMAX, ifp); // nmols
      sscanf(templine, "%d", &nmols);  //number of atoms
      fgets(templine, ONELINEMAX, ifp); // ITEM: BOX BOUND
      fgets(templine, ONELINEMAX, ifp); // xlo xhi
      sscanf(templine, "%f %f", &xlo, &xhi);
      fgets(templine, ONELINEMAX, ifp); // ylo yhi
      sscanf(templine, "%f %f", &ylo, &yhi);
      fgets(templine, ONELINEMAX, ifp); // zlo zhi
      sscanf(templine, "%f %f", &zlo, &zhi);
      fgets(templine, ONELINEMAX, ifp); //empty line
//	printf("Nmols=%d, finaloutputfilename=%s\n",ntotal,out);
       //fprintf(ofp,"# Maximum bin number %d\n ", BinNumber);	
       if(nmols<0) { std::cout <<"Error in ntotal: "<< ntotal << std::endl;  exit(0);}
       //std::cout << "Importing " << nmols << " atoms..." <<std::endl;

        for(i=0;i<nmols;i++) {
        fgets(templine, ONELINEMAX, ifp);
        sscanf(templine, "%d %d %f %f %f",&id, &type, &ix, &iy, &iz);
	allatom[frame][id]->id = id; 
	allatom[frame][id]->exist = 1; 
	allatom[frame][id]->type = type;
	allatom[frame][id]->frame = frame; 
	allatom[frame][id]->x[0] = ix;
	allatom[frame][id]->x[1] = iy;
	allatom[frame][id]->x[2] = iz;
	}
      frame++;
      }
    }
     fclose(ifp);
     }
        totalnmols=nmols;
	printf("TotalAtom #:%d\n",totalnmols);
	printf("It is done here1\n");
  //Cluster ananlysis & output
  

        printf("eventsize:%d\n",clusterevent.size());
	WearEvent(allatom,clusterevent,MaxAtoms,ntotal,cut,extension);

        printf("eventsize:%d\n",clusterevent.size());
int s=mergeevent(clusterevent,allatom,extension);

        printf("eventsize:%d,s=%d\n",clusterevent.size(),s);
	sortevent(clusterevent);

	//writetoxyzc_and_clusterDistribution(clusterevent,allatom,outputfilename);
	writetolata_and_clusterDistribution(clusterevent,allatom,outputfilename);
        printf("Done!!\n");
        printf("%d,%d\n",ntotal,MaxAtoms);
  //Deallocate memory
//  for( i=0;i<ntotal;i++) {
//    for( j=0;j<MaxAtoms;j++) {
//      delete [] allatom[i][j];
//    }
//    delete [] allatom[i];
//  }
//    delete [] allatom;    
}//main	
