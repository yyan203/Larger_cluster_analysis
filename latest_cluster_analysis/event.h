#ifndef __event_h_
#define __event_h_
#include <set>
#include <vector>
#include <utility>  //std::pair, std::make_pair
#include "atom.h"
#include "function.h"

class Event {  // this is about the information of every single atom
  int frame; // frame that atom become debris
  float framewidth; // frame to consider
  int clusteratomnum; // frame that atom become debris
  std::set< std::pair<int,int> > cluster;// atom associated (frame,ID) 
public:
  //float x[3];
  Event(int,int);
  void add_debri(int Atomid,int frame){cluster.insert(std::make_pair(Atomid,frame));}
  int get_frame() const {return frame;}
  float get_framewidth() const {return framewidth;}
  int get_atomnum() const {return cluster.size();}
  void set_frame(int f){frame=f;}
  void set_width(float f){framewidth=f;}
  const std::set<std::pair<int,int> > & getcluster() const{return cluster;}
  const  int & clusternum() const {return clusteratomnum;}
  bool addcluster(std::pair<int,int> newcluster){if(cluster.insert(newcluster).second){clusteratomnum=int(cluster.size());return true;}else{return false;}}

  void assignatomnum(){clusteratomnum=int(cluster.size());}
  void PrintInfo() {
    //printf("Event: frame[%d] type[%d] x[%f,%f,%f]\n", id,type,x[0],x[1],x[2]);
  };
};
///////////// cluster analysis function
bool smaller_cluster(const Event & event1, const Event & event2);
void WearEvent(Atom ***atom,std::vector <Event> & event, int total_nmols, int totalframe, float cutoff,int extension);
void sortevent(std::vector <Event> & event);

bool findcommonatom(const Event & event1, const Event & event2,  Atom ***atom, int ext);
void combineEvent(std::vector <Event> & event,int firstEvent, int secondEvent);
int mergeevent(std::vector <Event> & event, Atom ***atom, int ext);
void writetoxyzc_and_clusterDistribution(const std::vector <Event> & event, Atom ***atom,char *outputfilename);
void writetolata_and_clusterDistribution(const std::vector <Event> & event, Atom ***atom,char *outputfilename);

#endif
