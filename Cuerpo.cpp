# include <iostream>
# include <fstream>
# include <cmath>
# include <vector>
# include "Vector.h"
# include "Random64.h"
# include "Automata3D.h"
# include "Cuerpo.h"
# include "LB3D.h"
# include <math.h>
using namespace std ;
//----------------------------------CUERPO------------
Cuerpo::Cuerpo(void){}
Cuerpo::Cuerpo(double x0 ,double y0, double z0 ,double Vx0 ,double Vy0, double Vz0, double m0, double R0){//Constructor
  r.cargue(x0, y0, z0);
  V.cargue(Vx0, Vy0, Vz0);
  m = m0;
  R = R0;
  Existencia = true;
  Fg.cargue(0,-g,0);
}
Cuerpo::~Cuerpo(void){}
void Cuerpo::AgregaFuerzaArrastreyGravedad(LatticeBoltzmann & Agua){
    double aux;
  // --- Buoyancy Force ---
  aux=(4*M_PI*R0*R0*R0*(sand_density-density))/3;
  Fb=aux*Fg;
  F+=Fb;
  
  double x, y, z;
  x = r.x();
  y = r.y();
  z = r.z(); 
  if(0<x && 0 <z && 0<y){
    Vf = Agua.InterpolaVelocidad(x,y,z);
    RestaVelocidades(Vf);
    Fd = -6*M_PI*viscosity*R*V; // Fuerza de Arrastre
    SumaVelocidades(Vf);
    F += Fd;
  }
}
void Cuerpo::Mueva_r1(double dt){
  r += V*(chi*dt);
}
void Cuerpo::Mueva_V(LatticeBoltzmann & Agua,double dt){
  V += F*(dt/(2*m));
}
void Cuerpo::Mueva_r2(double dt){
  r += V*(Um2chi*dt);
}
void Cuerpo::DesplazaGaussiana(Crandom & ran2, double mu, double sigma){
  double ranx = ran2.gauss(mu ,sigma ), rany = ran2.gauss(mu, sigma), ranz = ran2.gauss(mu, sigma);
  Rgauss.cargue(ranx, rany, ranz);
  r += Rgauss;
}
void Cuerpo::RevisaSuelo(Automata & Suelo ,LatticeBoltzmann & Agua){
  int xprime, yprime, zprime;
  xprime = r.x()*nx; // Lleva a las unidades del suelo ( h [ i ])
  yprime = r.y()*ny;
  zprime = r.z()*nz;
  if(yprime <= Suelo.Get_h(xprime, zprime)){//Si está debajo del nivel suelo
    Suelo.AgregaGrano(xprime, zprime);  //Agrega un grano
    EliminaGrano();//Elimina de cuerpo
  }
}
