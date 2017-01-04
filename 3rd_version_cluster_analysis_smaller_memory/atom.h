#ifndef __atom_h_
#define __atom_h_
#include <fstream>
#include <set>
class Atom {  // this is about the information of every single atom
public:
  int id;
  int exist;
  int type;
  int frame; // frame that atom become debris
  float x[3];
  std::set<int> wear; // for cluster analysis
  void PrintInfo() {
    printf("Atom: id[%d] type[%d] x[%f,%f,%f]\n", id,type,x[0],x[1],x[2]);
  };
};
#endif
