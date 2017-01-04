// Event to store wear event
//#ifndef __event_h_
//#define __event_h_
#include <cmath>
#include "atom.h" 
#include "function.h" 
#define Boxx 91.8756 
#define Boxz 73.5004
#define AAcut 1.40
#define BBcut 1.20 
#define ABcut 1.30 
float GetDistance(const Atom *atom1, const Atom *atom2)
{
	double dx,dy,now_distance;
	double dz;
	double Lfield[3];
	Lfield[0]=Boxx;
	Lfield[2]=Boxz;
	dx=atom2->x[0]-atom1->x[0];
	dy=atom2->x[1]-atom1->x[1];
	dz=atom2->x[2]-atom1->x[2];
	while(dx > Lfield[0]*0.5) dx -= Lfield[0];
	while(dx <-Lfield[0]*0.5) dx += Lfield[0];
	//while(dy > Lfield[1]*0.5) dy -= Lfield[1];
	//while(dy <-Lfield[1]*0.5) dy += Lfield[1];
	while(dz > Lfield[2]*0.5) dz -= Lfield[2];
	while(dz <-Lfield[2]*0.5) dz += Lfield[2];
	now_distance = sqrt(dx*dx + dy*dy + dz*dz);
	return now_distance;
}
/// getcutoff
float GetCutoff(const Atom *atom1,const Atom *atom2)
{
	if(atom1->type == 1 && atom2->type ==1 ) {return AAcut;}
        if(atom1->type == 2 && atom2->type ==2 ) {return BBcut;}
        if(atom1->type == 1 && atom2->type ==2 ) {return ABcut;}
        if(atom1->type == 2 && atom2->type ==1 ) {return ABcut;}
}
//#endif
