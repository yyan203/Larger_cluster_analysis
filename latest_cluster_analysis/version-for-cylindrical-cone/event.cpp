// Event to store wear event
#include "event.h"
#include <algorithm>
#include "function.h"


Event::Event(int fram,int ext){frame=fram;framewidth=ext;}
//Event::Event(int fram,int ext, float x, float y, float z){frame=fram;x[0]=x;x[1]=y;x[2]=z;framewidth=ext;}

bool smaller_cluster(const Event & event1, const Event & event2){
return event1.clusternum()<event2.clusternum();
}

bool larger_cluster(const Event & event1, const Event & event2){
return event1.clusternum()>event2.clusternum();
}

void sortevent_ascend(std::vector <Event> & event)
//sort all cluster event
{
    std::sort(event.begin(),event.end(),smaller_cluster);
}

void sortevent_descend(std::vector <Event> & event)
//sort all cluster event
{
    std::sort(event.begin(),event.end(),larger_cluster);
}

void writetoxyzc_and_clusterDistribution(const std::vector <Event> & event, Atom ***atom, char *outputfilename)
{
  FILE *fp,*fp2;
  char clusterinfo[100];
  int id,frame;
  sprintf(clusterinfo, "%s.cluster", outputfilename);
  fp = fopen(outputfilename, "w");
  fp2 = fopen(clusterinfo, "w");

  for(int i=0;i<event.size();i++)
  {
        fprintf(fp,"%d\n#id type x y z eventframe\n",event[i].clusternum());
	for (std::set<std::pair<int,int> >::iterator it=event[i].getcluster().begin(); it!=event[i].getcluster().end(); ++it){
	frame=(*it).first;id=(*it).second;
	fprintf(fp, "%d %d %f %f %f %d\n", atom[frame][id]->id, atom[frame][id]->type, atom[frame][id]->x[0], atom[frame][id]->x[1], atom[frame][id]->x[2], atom[frame][id]->frame);
  }
	fprintf(fp2, "%d %d\n", i,event[i].clusternum());// By Yongjian cluster id and cluster atom number
  }
  fclose(fp);
  fclose(fp2);
}

void writetolata_and_clusterDistribution(const std::vector <Event> & event, Atom ***atom, char *outputfilename)
{
  FILE *fp,*fp2;
  char clusterinfo[100];
  int id,frame;
  sprintf(clusterinfo, "%s.cluster", outputfilename);
  fp = fopen(outputfilename, "w");
  fp2 = fopen(clusterinfo, "w");

  for(int i=0;i<event.size();i++)
  {
        fprintf(fp,"ITEM: TIMESTEP\n%d\n",event[i].get_frame());
        fprintf(fp,"ITEM: NUMBER OF ATOMS\n%d\n",event[i].clusternum());
        fprintf(fp,"ITEM: BOX BOUNDS pp ff pp\n0 91.8566\n0 55.1678\n0 73.5287\n");
        fprintf(fp,"ITEM: ATOMS id type x y z eventframe\n",event[i].clusternum());
	for (std::set<std::pair<int,int> >::iterator it=event[i].getcluster().begin(); it!=event[i].getcluster().end(); ++it){
	frame=(*it).first;id=(*it).second;
	fprintf(fp, "%d %d %f %f %f %d\n", atom[frame][id]->id, atom[frame][id]->type, atom[frame][id]->x[0], atom[frame][id]->x[1], atom[frame][id]->x[2], atom[frame][id]->frame);
  }
	fprintf(fp2, "%d %d\n", i,event[i].clusternum());// By Yongjian cluster id and cluster atom number
  }
  fclose(fp);
  fclose(fp2);
}

void writetolata_all_in_one(const std::vector <Event> & event, Atom ***atom, char *outputfilename)
{
  FILE *fp;
  char clusterinfo[100];
  int id,frame;
  sprintf(clusterinfo, "%s.all_cluster", outputfilename);
  fp = fopen(clusterinfo, "w");

  int total_atom=0;
  for(int i=0;i<event.size();i++){total_atom+=event[i].clusternum();}
        fprintf(fp,"ITEM: TIMESTEP\n0\n");
        fprintf(fp,"ITEM: NUMBER OF ATOMS\n%d\n",total_atom);
        fprintf(fp,"ITEM: BOX BOUNDS pp ff pp\n0 91.8566\n0 55.1678\n0 73.5287\n");
        fprintf(fp,"ITEM: ATOMS id type x y z eventframe ClusterSize clusterID\n");
  for(int i=0;i<event.size();i++)
  {
	for (std::set<std::pair<int,int> >::iterator it=event[i].getcluster().begin(); it!=event[i].getcluster().end(); ++it){
	frame=(*it).first;id=(*it).second;
	fprintf(fp, "%d %d %f %f %f %d %d %d\n", atom[frame][id]->id, atom[frame][id]->type, atom[frame][id]->x[0], atom[frame][id]->x[1], atom[frame][id]->x[2], atom[frame][id]->frame, event[i].clusternum(),i);
  }
  }
  fclose(fp);
}

bool findcommonatom(const Event & event1, const Event & event2,   Atom ***atom, int ext)
{
	for (std::set<std::pair<int,int> >::iterator it=event1.getcluster().begin(); it!=event1.getcluster().end(); ++it){
	 for (std::set<std::pair<int,int> >::iterator ij=event2.getcluster().begin(); ij!=event2.getcluster().end(); ++ij){
	     if(*it==*ij){ printf("equal\n"); return true;}
	     if((abs(atom[(*it).first][(*it).second]->frame-atom[(*ij).first][(*ij).second]->frame)<=ext)
		 && (GetDistance(atom[(*it).first][(*it).second],atom[(*ij).first][(*ij).second])<=(GetCutoff(atom[(*it).first][(*it).second],atom[(*ij).first][(*ij).second])))) {
		 //printf("close ");
		 return true;}
	 }
	}
	return false;
}

void combineEvent(std::vector <Event> & event,int firstEvent, int secondEvent){
	for (std::set<std::pair<int,int> >::iterator ij=event[secondEvent].getcluster().begin(); ij!=event[secondEvent].getcluster().end(); ++ij){
	event[firstEvent].addcluster(*ij);
	}
	int newframe=int(0.5*(event[firstEvent].get_frame()+event[secondEvent].get_frame()));
	float width=abs(newframe-event[firstEvent].get_frame());
	event[firstEvent].set_frame(newframe);event[firstEvent].set_width(width);
	//printf("combine:%d(%d)&%d(%d) ",firstEvent,event[firstEvent].getcluster().size(),secondEvent,event[secondEvent].getcluster().size(),width);
	//if(width>5)printf("    width%f     ",width);
	event.erase(event.begin()+secondEvent);
}

int mergeevent(std::vector <Event> & event, Atom ***atom, int ext)
//merge all clusters in event
//meanwhile delete all empty clusters
{
    if(event.size()==0) return 0;
    if(event.size()==1) {event[0].assignatomnum();return event.size();}
    else {
	// event double loop
	int flag=false;
	for (int i=0;i<event.size();i++){
	    for (int j=i+1;j<event.size();j++){
	// cluster double loop
	    //if(false) {combineEvent(event,i,j);flag=true;continue;}
	    if(findcommonatom(event[i],event[j],atom,ext)) {combineEvent(event,i,j);flag=true;continue;}
	    if(flag) continue;
	    }
	}
	if(flag){mergeevent(event,atom,ext);}
	else { for(int i=0;i<event.size();i++) event[i].assignatomnum();return event.size();}
    }
}

void WearEvent(Atom ***atom,std::vector <Event> & event, int total_nmols, int totalframe, float cutoff,int extension)
  //Bondmax is the max bond length to be considered to be adjacent
  //cluster here, the first ncluster is non-emtpy!, need to move them in order to have that.
{
  int nmols=total_nmols;
  //float xx[3];
  int frame1;
  int count;
  int totalevent=0;
  int track[total_nmols];
  //printf("eventsize:%d ",event.size());
  for (int i=0;i<total_nmols;i++) track[i]=0;
  for(int f=0;f<totalframe;f++) // loop frame
  {
      count=0;
      for(int j=0;j<total_nmols;j++){
	if(atom[f][j]->exist>0){
	frame1=track[j];
	//printf("dist:%f ",GetDistance(atom[i][j], atom[frame1][j]));
        if(GetDistance(atom[f][j], atom[frame1][j])> cutoff)
	{
	    count++;
	//printf("dist:%f-frame:%d-id:%d ",GetDistance(atom[i][j], atom[frame1][j]),i,j);
	track[j]=f;
	//printf("changeAtom%dFrameTo%d ",j,f);
	atom[f][j]->wear.insert(f);
	Event newEvent(f,extension); 
	newEvent.add_debri(f,j);
	event.push_back(newEvent);
	} 
      }
      }
      //if(count>0)printf("Frame%devent#%d\n",f,count);
      totalevent+=count;
  }
    //printf("Totalevent#%d\n",totalevent);
   //printf("eventsize:%d ",event.size());
}
