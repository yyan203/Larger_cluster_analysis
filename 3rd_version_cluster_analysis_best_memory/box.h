#ifndef __box_h_
#define __box_h_
#include <fstream>
#include <set>
#include <cmath>
class Box {  // this is about the information of every single box
  float x[2];
  float y[2];
  float z[2];
  float l[3];
public:
  float get_xl() const {return x[0];}
  float get_xh() const {return x[1];}
  float get_yl() const {return y[0];}
  float get_yh() const {return y[1];}
  float get_zl() const {return z[0];}
  float get_zh() const {return z[1];}

  float get_lx() const {return l[0];}
  float get_ly() const {return l[1];}
  float get_lz() const {return l[2];}

  void set_x(float xlo, float xhi){x[0]=xlo;x[1]=xhi;l[0]=std::abs(xhi-xlo);}
  void set_y(float ylo, float yhi){y[0]=ylo;y[1]=yhi;l[1]=std::abs(yhi-ylo);}
  void set_z(float zlo, float zhi){z[0]=zlo;z[1]=zhi;l[2]=std::abs(zhi-zlo);}

  void PrintInfo() {
    printf("Box: xlo[%f] xhi[%f] ylo[%f] yhi[%f] zlo[%f] zhi[%f]\n", x[0],x[1],y[0],y[1],z[0],z[1]);
  }
};
#endif
