
// COMPUTE Neighbor Debris Atom probibility according to distance cutoff (species dependent) & time duration 
// DO CLUSTER ANALYSIS on Debris Atoms according to distance cutoff (species dependent) & time duration 
//
//
// Input file format is Olivia type:
// Simulation_X Simulation_Y Simulation_Z  Shift
// ITEM: ATOMS id type x y z vx vy vz fx fy fz

#include <cmath>
#include <fstream>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <cstdlib>

#define ONELINEMAX 2000
#define MAXDUMP 1000
#define MaxAtoms 200000

#define MaxClusters 130000
#define EMPTY -100//for linkcell and other purpose
#define VERBOSE 1//for print debug 

#define AAcut 1.40
#define BBcut 1.20 
#define ABcut 1.30 

#define TRUE 1
#define FALSE 0



//using namespace std;

/////////////////////////////////
// Global variables defination //
/////////////////////////////////

class Atom {  // this is about the information of every single atom
public:
  int id;
  int type;
  int frame; // frame that atom become debris
  float x[3];
  float v[3];
  float f[3];
  int next; // for cluster analysis
  void PrintInfo() {
    printf("Atom: id[%d] type[%d] x[%f,%f,%f]\n", id,type,x[0],x[1],x[2]);
  };
};

// check if a file exist
bool fexists(const char *filename)
{
  std::ifstream ifile(filename);
    return ifile;
}


// Cluster Analysis

void clearcandidate(int *candidate)
{
  for(int i=0;i<MaxClusters;i++)
    candidate[i] = EMPTY;
}


int checkcandidate(int *candidate, int ncandidate, int candidateindex)
  //if candidateindex in candidate, then return TRUE
{
  for(int i=0;i<ncandidate;i++)
    if(candidate[i]==candidateindex) return TRUE;
  
  return FALSE;
}



int addatomtocluster(int *cluster, int ncluster, struct Atom **atom, int nmols, int clusterindex, int atomindex)
{
  int ifnewcluster=0;
  //atom should not be bounded to any cluster now
  if(atom[atomindex]->next!=EMPTY) {printf("Atom(%d) is bounded to other cluster already!\n", atomindex);exit(0);}

  if(cluster[clusterindex]==EMPTY) ifnewcluster=1;

  atom[atomindex]->next = cluster[clusterindex];
  cluster[clusterindex] = atomindex;

#if VERBOSE==1
  //  printf("\t atom(%d) add to cluster(%d)\n", atomindex, clusterindex);
#endif
  return ifnewcluster;
}

void addclustertocluster(int *cluster, int ncluster, struct Atom **atom, int nmols, int target, int source)
  //add source to target
{
  int atomindex = cluster[target];

  if(atomindex==EMPTY)//target is empty
    {
      cluster[target] = cluster[source];
      cluster[source] = EMPTY;
      return;
    }

  //advance to the tail of target cluster
  while(atom[atomindex]->next != EMPTY)
    {
      atomindex = atom[atomindex]->next;
    }
  //now atom[atomindex]->next is EMPTY, means it is the last atom in the target cluster
  atom[atomindex]->next = cluster[source];
  cluster[source] = EMPTY;
}


int validcluster(int *cluster, int ncluster)
{
  //recalculate ncluster
  int clustercount, i;
  clustercount=0;
  for(i=0;i<ncluster;i++)
    if(cluster[i]!=EMPTY) clustercount++;		
  return clustercount;
}



void mergecluster(int *cluster, int ncluster, struct Atom **atom, int nmols,  int *candidate, int ncandidate)
//merge all clusters in candidate
//meanwhile delete all empty clusters
{
  int i, j, temp, count=0;

  //merge to the first one
  for(i=1;i<ncandidate;i++)
    {
      addclustertocluster(cluster, ncluster, atom, nmols, candidate[0], candidate[i]);
    }

  //Sort cluster so that empty clusters at end.
  for(i=0;i<ncluster-1;i++)
    for(j=i+1;j<ncluster;j++)
      if(cluster[i]==EMPTY && cluster[j]!=EMPTY) {temp = cluster[i]; cluster[i] = cluster[j]; cluster[j] = temp;}

  for(i=0;i<ncluster;i++)
    if(cluster[i]!=EMPTY) count++;

#if VERBOSE==1
  //  printf("\t\t merge clusters: %d => %d: ", ncluster, count);
  //  for(i=0;i<ncluster;i++) printf("%d ", cluster[i]);
  //  printf("\n");
#endif

}

// Use the number of shared atoms to determine if they are connected!

int GetClusterSize(Atom **atom, int clusterindex, int *cluster)
{
  int count=0;
  int atomindex=cluster[clusterindex];
  while(atomindex!=EMPTY)
    {
      count++;
      atomindex = atom[atomindex]->next;
    };
  return count;
}

void SortCluster(Atom **atom, int *cluster,int ncluster)
{
  int i, j, temp;
  for(i=0;i<ncluster-1;i++)
    for(j=i+1;j<ncluster;j++)
      {
	if(GetClusterSize(atom, i, cluster) < GetClusterSize(atom, j, cluster))
	  {
	    temp = cluster[i];
	    cluster[i] = cluster[j];
	    cluster[j] = temp;
	  }
      }
}



//output xyz
void writetoxyzc_and_cluster_local(int *cluster, int ncluster, struct Atom **atom, int nmols, char *outputfilename, int noutput, int toframe, int extension)
{
  FILE *fp,*fp2;
  int count, i,j, atomindex;
  char clusterinfo[100];
  int flag=0;

  sprintf(clusterinfo, "%s.cluster", outputfilename);
  fp = fopen(outputfilename, "w");
  fp2 = fopen(clusterinfo, "w");
 
  
  fprintf(fp, "%d\n#\n", noutput);//xyzc format!
  

	j=0;
  for(i=0;i<ncluster;i++)
    if(cluster[i]!=EMPTY)
      {
        count = 0;// By Yongjian
	flag = 0;
	atomindex = cluster[i];
	while(atomindex!=EMPTY)
	  {
	    fprintf(fp, "%d %d %f %f %f %d", atom[atomindex]->id, atom[atomindex]->type, atom[atomindex]->x[0], atom[atomindex]->x[1], atom[atomindex]->x[2], atom[atomindex]->frame);
	    fprintf(fp, " %d\n", i);//cluster name as "color"
	    count++;

	    if (atom[atomindex]->frame<= (toframe-extension)){flag=1;}  // By Yongjian for cutting the tailing cluster 

	    atomindex = atom[atomindex]->next;
	  }

	//fprintf(fp2, "%d %d\n", i+1,count);
	if(flag==1){j++;fprintf(fp2, "%d %d\n", j,count);}// By Yongjian cluster id and cluster atom number
      }
#if VERBOSE==1
  printf("plan to output %d atoms, while outputted %d atoms in total.\n", noutput, count);
#endif
  fclose(fp);
  fclose(fp2);
}








double GetDistance(Atom *atom1,Atom *atom2)
{
	double dx,dy,now_distance;
	double dz;
	dx=atom2->x[0]-atom1->x[0];
	dy=atom2->x[1]-atom1->x[1];
	dz=atom2->x[2]-atom1->x[2];
	//while(dy > Lfield[1]*0.5) dy -= Lfield[1];
	//while(dy <-Lfield[1]*0.5) dy += Lfield[1];
	//while(dx > Lfield[0]*0.5) dx -= Lfield[0];
	//while(dx <-Lfield[0]*0.5) dx += Lfield[0];
	//while(dz > Lfield[2]*0.5) dz -= Lfield[2];
	//while(dz <-Lfield[2]*0.5) dz += Lfield[2];
	now_distance = sqrt(dx*dx + dy*dy + dz*dz);
	return now_distance;
}


double GetCutoff(Atom *atom1,Atom *atom2)
{
	if(atom1->type == 0 && atom2->type ==0 ) {return AAcut;}
        if(atom1->type == 1 && atom2->type ==1 ) {return BBcut;}
        if(atom1->type == 0 && atom2->type ==1 ) {return ABcut;}
        if(atom1->type == 1 && atom2->type ==0 ) {return ABcut;}
	
}







void IDClusters(struct Atom **atom, char *outputfilename, int total_nmols, int fromframe, int toframe, int extension)
  //Vmin is the min shear will be considered to be clustered
  //Bondmax is the max bond length to be considered to be adjacent
  //cluster here, the first ncluster is non-emtpy!, need to move them in order to have that.
{
  int i, j, count=0, *candidate, ncandidate, clustercount=0;
    int *cluster;//store first atom index if not empty
    int ncluster=0;//the number of cluster now
  char tempstring[100];
  int nmols=total_nmols;

#if VERBOSE==1
  printf("Start Cluster Analysis\n");
#endif
  //initiate all variables
  cluster = new int[MaxClusters];
  //  clusterpe = new double[MaxClusters];
  candidate = new int[MaxClusters];
  for(i=0;i<MaxClusters;i++) {cluster[i]=EMPTY;candidate[i]=EMPTY;}
  ncluster=0;ncandidate=0;
  for(i=0;i<nmols;i++) atom[i]->next=EMPTY;

  //first identify all the sheared atom ( larger than Vmin )
/*  for(i=0;i<nmols;i++) 
    if((STZType==1&&atom[i]->shear > Vmin) || (STZType==2&&atom[i]->D > Vmin)) count++;
*/
#if VERBOSE==1
  printf("There are %d out of %d atoms considered to be sheared atom\n", count, nmols);
#endif

  //then loop over all sheared atom to put them into clusters
  for(i=0;i<nmols;i++)
      {
	clearcandidate(candidate); ncandidate=0;

	//Now loop over existing clusters, and all atoms already in the cluster
	//if atom[i] is within the bondlength of an atom, then add atom[i] to the current atom
	//current solution!!!
	//should be loop over all clusters, then merge clusters (maybe a number of them) can contain atom[i]
	for(j=0;j<ncluster;j++)//loop over non-empty clusters only	 
	  {
	    int pi = cluster[j]; //pi is the atom index in cluster[j]
	    while(pi != EMPTY)
	      {
		//check distance between atom[i] and atom[pi]
		//originally using Bondmax, now using the nearestnbrceil function, march16, 04, for percolation of stable atoms
		//if( GetDistance(atom[i], atom[pi])< Bondmax )
		if( GetDistance(atom[i], atom[pi]) < GetCutoff(atom[i], atom[pi]) && abs(atom[i]->frame - atom[pi]->frame) <= extension)
		//if( GetDistance(atom[i], atom[pi]) < Bondmax * GetNBRCeil(atom[i], atom[pi],this))//bondmax become a scaling factor
		{		  
		  if(checkcandidate(candidate, ncandidate, j)==FALSE)
		    //add to candidate only if j is not in candidate
		    {
		      candidate[ncandidate] = j; ncandidate++;
		      if(ncandidate>= MaxClusters) {printf("Increase MaxClusters!\n");exit(0);}
#if VERBOSE == 1
		      //		      printf("atom i(%d) found cluster(%d)\n", i, j);
#endif
		    }
		}
		//continue on next atom
		pi = atom[pi]->next;
	      }	    
	  }

	//after looping all clusters
	// add atom[i] into the candidate cluster
	// if no proper candidate cluster, add a new cluster
	if(ncandidate>0)
	  {
	    addatomtocluster(cluster, ncluster, atom, nmols, candidate[0], i);
	    if(ncandidate>1)//merge all candidate clusters together
	      {
		mergecluster(cluster, ncluster, atom, nmols, candidate, ncandidate);
		ncluster = validcluster(cluster, ncluster);
	      }
	  }
	else //there is no candidate, then create another cluster
	  {
	    addatomtocluster(cluster, ncluster, atom, nmols, ncluster, i);ncluster++;
	    if(ncluster>= MaxClusters) {printf("Increase MaxClusters!\n");exit(0);}
	  }
      }


#if VERBOSE==1
  printf("There are totally %d clusters found\n", ncluster);
  //  for(i=0;i<ncluster;i++)
  //    printclusterinfo(cluster, ncluster, atom, nmols, i);
#endif

  //delete all clusters have only one atom
  //  DeleteClusters(cluster, ncluster, 1);

  //sort clusters
  SortCluster(atom,cluster,ncluster);

  //output them
  sprintf(tempstring, "%s.xyzc", outputfilename);
  writetoxyzc_and_cluster_local(cluster, ncluster, atom, nmols, tempstring, nmols, toframe, extension);
  //sprintf(tempstring, "%s.cluster", outputfilename);
  //writecluster(atom, tempstring, count, STZType);

  //sprintf(tempstring, "%s.largest.xyzc", outputfilename);
  //writetoxyzc_shift_big(cluster, 1, atom, nmols, tempstring, count);

  //clear up and return
  //  delete [] cluster;
  delete [] candidate;
  return;
}

// Cluster Analysis Ends here!





int main(int argc, char *argv[])
{

  Atom **debriatom;
  char oliviafilename[100], outputfilename[100],out2[100];
  FILE *ifp, *ofp, *ofp2;
  //Frame *myframe;
  char templine[ONELINEMAX], str1[100], str2[100];
  int count_frame_Atom[MAXDUMP]={0},count_frame_Correlation[MAXDUMP]={0};
  float xlo, xhi, ylo, yhi, zlo, zhi;
  float dx, dy,dz,distance; 
  float ix, iy, iz, ivx, ivy, ivz, ifx, ify, ifz;
  int id, type, startf, endf;
  float avgY, edgeY,maxF,minF,binl;
  char out[200];
  int i,j,k,nmols=0;
  int n;
  int flag;
  int totalnmols=0;
  int countAtom=0, countCorrelation=0;
  int fromframe, Step, toframe, extension;
  //float AAcut, BBcut, ABcut;


   if(argc>=3)
    {
      sscanf(argv[1], "%s", oliviafilename);
      sscanf(argv[2], "%s", outputfilename);
      sscanf(argv[3], "%d", &fromframe);
      sscanf(argv[4], "%d", &Step);
      sscanf(argv[5], "%d", &toframe);
      //sscanf(argv[6], "%f", &AAcut);
      //sscanf(argv[7], "%f", &BBcut);
      //sscanf(argv[8], "%f", &ABcut);
      sscanf(argv[6], "%d", &extension);
    }

   else

    {
      std::cout<<"Correct syntax: "<< argv[0]<< " oliviafile(only leading Name without frame number) outputfilename fromframe Step toframe Extension(consider atoms correlation extended to how many frames after current one)\n";
      std::cout<<"Pay Attention to this: Output cluster info only contains clusters containing at least one atoms belonging to frame [fromframe to (toframe-extension)],(inclusive))\n";
      std::cout<<"Pay Attention to this: Total Debris Atoms must less than 200000, otherwise please increase it in the C++ code by changing constant: MaxAtoms\n";
      std::cout<<"Pay Attention to this: Different Atom Species Pairs have different distance cutoff,they are defined at the beginning Constant!!! Change them to suit your use!!\n";
      std::cout<<"Pay Attention to this: Atom type in oliviafile may have been shifted by -1 because the debriatom olivia file is generated by another code which may have changed the type to suit the MapCorrelation function calculation, so you need to change the current code (compare) to make sure the code works right!!\n";
      exit(0);
    }

    
    debriatom = new Atom *[MaxAtoms];
    for(i=0;i<MaxAtoms;i++) 
    {
    debriatom[i] = new Atom();
    }


    printf("Initialization of derbisatom has been finished! \n");
    
    //ofp = fopen(outputfilename,"w");

     // sprintf(out2,"Debri_Possibility_Evolution_vs_frame.dat");

    //ofp2 = fopen(out2,"w");

	
    for(i=fromframe;i<=toframe;i=i+Step)
	
    { 

     sprintf(out,"%s.%d",oliviafilename,i);

     if(fexists(out))
     {
     printf("file:%s exist\n",out);

     flag=0;
     ifp = fopen(out,"r");
     if(ifp==NULL) { std::cout << oliviafilename << " not exist! \n"; exit(0);} 

     printf("It is done here 0\n");

     fgets(templine, ONELINEMAX, ifp);flag++;// Olivia First Line 

      std::cout<< "start reading atom in this frame: "<<i<<"\n";
	 // if(atom==NULL) printf("Error in allocating memory!\n");
	 // if(newatom==NULL) printf("Error in allocating memory!\n");
      n=0;
      while(fgets(templine, ONELINEMAX, ifp) != NULL)
      {sscanf(templine, "%d %d %f %f %f", &id, &type, &ix, &iy, &iz);
        n++;
	debriatom[nmols]->id = id; 
	debriatom[nmols]->type = type;
	debriatom[nmols]->frame = i; 
	debriatom[nmols]->x[0] = ix;
	debriatom[nmols]->x[1] = iy;
	debriatom[nmols]->x[2] = iz;
	debriatom[nmols]->next = EMPTY;

	printf ("%dth Atom is imported!\n",n);

      nmols++;
      }//while

     fclose(ifp);
     }

    }//for loop
	printf("It is done here1\n");

    totalnmols=nmols;

	printf("TotalAtom #:%d\n",totalnmols);


	
	IDClusters(debriatom,outputfilename,totalnmols,fromframe,toframe,extension);


	  //Do Stuff

/*
     for(j=0;j<totalnmols;j++)
     {
	printf("Atom index:%d\t",j);
	if(debriatom[j]->frame >= fromframe+extension && debriatom[j]->frame <= toframe-extension)
	{
	count_frame_Atom[debriatom[j]->frame]++;
	countAtom++;
	for(k=0;k<totalnmols;k++)
	{
	if(k!=j && abs(debriatom[j]->frame - debriatom[k]->frame) <= extension)
	{
	dx=debriatom[j]->x[0] - debriatom[k]->x[0];
	dy=debriatom[j]->x[1] - debriatom[k]->x[1];
	dz=debriatom[j]->x[2] - debriatom[k]->x[2];
	distance = sqrt (dx*dx+dy*dy+dz*dz);
	printf("distance is:%f \t",distance);
	if(debriatom[j]->type == 0 && debriatom[k]->type ==0 && distance <= AAcut) {countCorrelation++;count_frame_Correlation[debriatom[j]->frame]++;break;} 
	if(debriatom[j]->type == 1 && debriatom[k]->type ==1 && distance <= BBcut) {countCorrelation++;count_frame_Correlation[debriatom[j]->frame]++;break;} 
	if(debriatom[j]->type == 0 && debriatom[k]->type ==1 && distance <= ABcut) {countCorrelation++;count_frame_Correlation[debriatom[j]->frame]++;break;} 
	if(debriatom[j]->type == 1 && debriatom[k]->type ==0 && distance <= ABcut) {countCorrelation++;count_frame_Correlation[debriatom[j]->frame]++;break;} 
	}
	}
	
	}
     
     }//for loop



//	printf("It is done here2\n");
	printf("countAtom %d countCorrelation %d Neighbour DebrisEventPossibility %f\n ", countAtom, countCorrelation,float(countCorrelation)/float(countAtom));

        fprintf(ofp,"# Atoms considered from %d to %d, Delta-frame:%d\n ",fromframe,toframe,extension);
        fprintf(ofp,"# Column 1: Atom#\nColumn 2:CountRelated Time\nColumn 3:Possibility to find another debris event\n");
        fprintf(ofp,"%d %d %f\n ", countAtom, countCorrelation,float(countCorrelation)/float(countAtom));
	
        fclose(ofp);

	fprintf(ofp2,"# Column1 frame\n#Column2 DebrisAtom#\n#Column3 CountRelated\n#Column4 CountRelated_Percent(each frame)\n#Column5 Possibility to find another debris event\n");

	for (j=fromframe;j<=toframe;j++) 
	{
	if(count_frame_Atom[j]>0){fprintf(ofp2,"%d %d %d %f %f\n ", j,count_frame_Atom[j], count_frame_Correlation[j],float(count_frame_Correlation[j])/float(countCorrelation),float(count_frame_Correlation[j])/float(count_frame_Atom[j]));}
	}

	fclose(ofp2);

*/

}//main	

	



