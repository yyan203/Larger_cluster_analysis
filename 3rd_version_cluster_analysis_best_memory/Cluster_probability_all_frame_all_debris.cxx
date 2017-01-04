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
//#define MaxAtoms 1600000
#define MaxClusters 1000000
#define EMPTY -100//for linkcell and other purpose
#define VERBOSE 1//for print debug 
#define TRUE 1
#define FALSE 0

#include "atom.h"
#include "box.h"
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
  Box box;
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
  float ivx, ivy, ivz;
  float ifx, ify, ifz;
  float cut;
  int id, type;
  float edgeY,maxF,minF,binl;
  char out[200];
  int i,j,k,nmols=-1;
  int n;
  int flag;
  int totalnmols=0;
  int frame;
  unsigned long int MaxAtoms;
  //int countAtom=0, countCorrelation=0;
  int countAtom=0;
  int fromframe, Step, toframe, frameStep;
   if(argc>=7)
    {
     sscanf(argv[1], "%s", oliviafilename);
     sscanf(argv[2], "%s", outputfilename);
     sscanf(argv[3], "%d", &fromframe);
     sscanf(argv[4], "%d", &Step);
     sscanf(argv[5], "%d", &toframe);
     sscanf(argv[6], "%d", &frameStep);
     sscanf(argv[7], "%f", &cut);
     sscanf(argv[8], "%lu", &MaxAtoms);
    }
    else
    {
      std::cout<<"Correct syntax: "<< argv[0]<< " lata4olivia outputfilename fromframe Step toframe Extension(consider atoms correlation extended to how many frames after current one) Cutoff(define a debris event) MaxAtomNum\n";
      std::cout<<"!! Total Debris Atoms must less than 200000, otherwise please increase it in the C++ code by changing last argument: MaxAtoms\n";
      std::cout<<"!! Different Atom Species Pairs have different distance cutoff,they are defined in function.h!!!\n";
      exit(0);
    }
//initialized container
      ntotal=(toframe-fromframe)/Step/frameStep+1;
      allatom = new Atom **[ntotal];
      for(i=0;i<ntotal;i++) 
      {
      allatom[i] = new Atom *[MaxAtoms];
      for(j=0;j<MaxAtoms;j++) 
      {allatom[i][j] = new Atom();
      allatom[i][j]->exist=0;}
      }
      printf("here! \n");
      int* table= new int[20000000];
      int index;
      printf("Initialization of atom has been finished! \n");
      sprintf(out,"%s",oliviafilename);

      if(fexists(out))
      {printf("file:%s exist\n",out);
      flag=0;
      ifp = fopen(out,"r");
      //if(ifp==NULL) { std::cout << oliviafilename << " not exist! \n"; exit(0);} 
      printf("It is done here 0\n");
      //fgets(templine, ONELINEMAX, ifp);flag++;// Olivia First Line 
      //if(atom==NULL) printf("Error in allocating memory!\n");
      //if(newatom==NULL) printf("Error in allocating memory!\n");
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

      if(frame==0){box.set_x(xlo,xhi);box.set_y(ylo,yhi);box.set_z(zlo,zhi);}

       //printf("Nmols=%d, finaloutputfilename=%s\n",ntotal,out);
       //fprintf(ofp,"# Maximum bin number %d\n ", BinNumber);	
       if(nmols<0) { std::cout <<"Error in ntotal: "<< ntotal << std::endl;  exit(0);}
       //std::cout << "Importing " << nmols << " atoms..." <<std::endl;
       if(current_timestep >= fromframe && (current_timestep-fromframe)%(Step*frameStep)==0)
       {for(i=0;i<nmols;i++){
        fgets(templine, ONELINEMAX, ifp);
        //sscanf(templine, "%d %d %f %f %f",&id, &type, &ix, &iy, &iz,&ivx, &ivy, &ivz, &ifx, &ify, &ifz);
        sscanf(templine, "%d %d %f %f %f",&id, &type, &ix, &iy, &iz);
	if(type<3){
	if(frame==0){
	table[id]=countAtom;
	allatom[frame][countAtom]->id = id; 
	allatom[frame][countAtom]->exist = 1; 
	allatom[frame][countAtom]->type = type;
	allatom[frame][countAtom]->frame = frame; 
	allatom[frame][countAtom]->x[0] = ix;
	allatom[frame][countAtom]->x[1] = iy;
	allatom[frame][countAtom]->x[2] = iz;
	countAtom++;
	}
	else{
	index=table[id];
	allatom[frame][index]->id = id; 
	allatom[frame][index]->exist = 1; 
	allatom[frame][index]->type = type;
	allatom[frame][index]->frame = frame; 
	allatom[frame][index]->x[0] = ix;
	allatom[frame][index]->x[1] = iy;
	allatom[frame][index]->x[2] = iz;
	}
	}
	}
      frame++;
      } else {for(i=0;i<nmols;i++){fgets(templine, ONELINEMAX, ifp);}}
      }
     }
     fclose(ifp);
     }
        totalnmols=nmols;
	printf("TotalAtom #:%d\n",totalnmols);
	printf("It is done here1\n");
        //Cluster ananlysis & output
	int* clusterlist= new int[MaxClusters];
	//loop frame
	for(i=0;i<=frame-2;i+=1) {
	WearEventStep(allatom,clusterevent,countAtom,i,1,cut,box);

        mergeevent(clusterevent,allatom,0,box);

	sortevent_ascend(clusterevent);

	writetolata_step(clusterevent,allatom,outputfilename,i*frameStep,box);

	statistic_event(clusterevent,clusterlist);

        printf("frame:%d,eventsize:%d\n",i,clusterevent.size());

	clusterevent.clear();
	}
	
	write_clusterDistribution(clusterlist,MaxClusters,outputfilename);


////////////////////////////////
//      printf("eventsize:%d\n",clusterevent.size());
//	WearEvent(allatom,clusterevent,MaxAtoms,ntotal,cut,extension);
//
//      printf("eventsize:%d\n",clusterevent.size());
//      int s=mergeevent(clusterevent,allatom,extension);
//
//      printf("eventsize:%d,s=%d\n",clusterevent.size(),s);
//	sortevent_ascend(clusterevent);
//
//	//writetoxyzc_and_clusterDistribution(clusterevent,allatom,outputfilename);
//	writetolata_and_clusterDistribution(clusterevent,allatom,outputfilename);
//
//	sortevent_descend(clusterevent);
//	writetolata_all_in_one(clusterevent,allatom,outputfilename);
//      printf("Done!!\n");
//      printf("%d,%d\n",ntotal,MaxAtoms);
///////////////////////////////////
//
//Deallocate memory
//  for( i=0;i<ntotal;i++) {
//    for( j=0;j<MaxAtoms;j++) {
//      delete [] allatom[i][j];
//    }
//    delete [] allatom[i];
//  }
//    delete [] allatom;    

}
