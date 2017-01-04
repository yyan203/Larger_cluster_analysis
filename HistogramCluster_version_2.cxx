// HistogramQuantity.cxx
// by yangy14@rpi.edu 14
//
// Main function is to do histogram analysis of a certain quantity 
// input file format (column-wise, e.g, cluster histogram analysis): 

// clusterid x y z numberofatomsinthecluster X X X
// clusterid x y z numberofatomsinthecluster X X X

// (same as output of MapCorrelation3 28)



#include <math.h>
#include "stdio.h"
#include "stdlib.h"



#define MaxQuantity 130001



int main (int argc, char ** argv)
{


if(argc<5) { printf("Correct usage:\t%s ",argv[0]);
   	     printf("input output TotalColumnNumber QuantityColumn (Column that the quantity located) \n");
	     printf("Description: Do Histogram of Some Quantity\t \n");
	     exit(0);
	        }

	FILE *fp, *ofp,*ofp2;
	char inputfilename[100], outputfilename[100],outputfilename2[100];
	float tmpf1, tmpf3, tmpf4,  tmpf6, tmpf7, tmpf8;
	int tmpf2;
	int TotalCount=0; 
	int TotalAtom=0; 
	int i,ColumNbr;    
	int NthColum;    

	int  *bin = new int[MaxQuantity]; 
	int  *bin_atom = new int[MaxQuantity]; 
       	for(i=0;i<MaxQuantity;i++) {bin[i]=0;bin_atom[i]=0;} 

	int  *StepBin = new int[5];  // Larger and unhomogeneous Bin Size to divide Clusters 
       	for(i=0;i<5;i++) {StepBin[i]=0;} 

	int  *ClusterAtomNbr = new int[5]; 
       	for(i=0;i<5;i++) {ClusterAtomNbr[i]=0;} 



	sscanf(argv[1], "%s", inputfilename); 
	sscanf(argv[2], "%s", outputfilename);
	sscanf(argv[3], "%d", &ColumNbr);
	sscanf(argv[4], "%d", &NthColum);

	fp = fopen(inputfilename, "r");
	ofp = fopen(outputfilename, "w");
 	sprintf(outputfilename2, "%s.LargerBin.dat", outputfilename);	
	ofp2 = fopen(outputfilename2, "w");

	printf("\t\t====%s==== \n", inputfilename);
	if(fp==NULL)
	{
	       printf("error reading file %s \n",inputfilename);
	       exit(0);
	}
	//	char lataoneline[1000];
	while(!feof(fp))
	    {
	//    fgets(lataoneline, 1000, fp);
	    fscanf(fp, "%f %d \n", &tmpf1, &tmpf2);
	TotalCount++;
	TotalAtom+=tmpf2;
	    if (tmpf2>30000){bin[30001]++;bin_atom[30001]=bin_atom[30001]+tmpf2;}
	    else {bin[tmpf2]++;bin_atom[tmpf2]=bin_atom[tmpf2]+tmpf2;}
	    
	    if (tmpf2<1.5 && tmpf2>0){StepBin[0]++;}
	    if (tmpf2<2.5 && tmpf2>1.5){StepBin[1]++;}
	    if (tmpf2>2.5 && tmpf2<10.5) {StepBin[2]++;}
	    if (tmpf2>10.5 && tmpf2<100.5) {StepBin[3]++;}
	    if (tmpf2>100.5 ) {StepBin[4]++;}

	    if (tmpf2<1.5 && tmpf2>0){ClusterAtomNbr[0]+=tmpf2;}
	    if (tmpf2<2.5 && tmpf2>1.5){ClusterAtomNbr[1]+=tmpf2;}
	    if (tmpf2>2.5 && tmpf2<10.5) {ClusterAtomNbr[2]+=tmpf2;}
	    if (tmpf2>10.5 && tmpf2<100.5) {ClusterAtomNbr[3]+=tmpf2;}
	    if (tmpf2>100.5 ) {ClusterAtomNbr[4]+=tmpf2;}
            }
	
	printf("\t\t====TotalCount:%d==== \n",  TotalCount);
	// output histogram results
	fprintf(ofp, "#Column1:ClusterSize(number of Atoms contained) #Column2:ClusterNumber #Column3:portion #Column4:ln(portion) #Column5:Total Atoms# in this sized cluster #Column6:atom portion\n");
	for(i=0;i<MaxQuantity;i++) {if (bin[i]>0) fprintf(ofp, "%d %d %f %f %d %f\n", i, bin[i], (float)(bin[i])/(float)(TotalCount),log((float)(bin[i])/(float)(TotalCount)),bin_atom[i],(float)(bin_atom[i])/(float)(TotalAtom));}
	fprintf(ofp2, "#TotalCluster#:%d TotalAtom#:%d\n#c1:Index\n#c2:Cluster#\n#c3:Atom#\n#c4:ClusterProportion\n#c5:AtomProportion\n#c6:Spread_index\n",TotalCount,TotalAtom);
	for(i=0;i<5;i++) {if(StepBin[i]>0) fprintf(ofp2, "%d %d %d %f %f %f\n", i, StepBin[i],ClusterAtomNbr[i], (float)(StepBin[i])/(float)(TotalCount),(float)(ClusterAtomNbr[i])/(float)(TotalAtom),(float)(StepBin[i])/(float)(ClusterAtomNbr[i]));}
	
	fclose(fp);
	fclose(ofp);
	fclose(ofp2);
}  // for main
