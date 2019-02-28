/*
Clase que calcula fuerzas para las partículas inmersas en el fluido
 */
# include <iostream>
# include <fstream>
# include <cmath>
# include <vector>
# include "Vector.h"
# include "Random64.h"
# include "Cuerpo.h"
# include "Colisionador.h"
# include "LB3D.h"
# include <math.h>
using namespace std ;
//----------------------------------COLISIONADOR---------
void Colisionador::CalculeFuerzasDelFluido(vector<Cuerpo> &vectorCuerpo ,LatticeBoltzmann & Agua){
  for(Cuerpo & obj: vectorCuerpo)obj.BorraFuerza();
  for(Cuerpo & obj: vectorCuerpo)obj.AgregaFuerzaArrastreyGravedad(Agua);
  for(Cuerpo & obj: vectorCuerpo)AgregaFuerzaObstaculo(obj);
}
void Colisionador::AgregaFuerzaObstaculo(Cuerpo & Grano ){
  vector3D centro , Vc ,n , r21 , Vcn;
  double normaVcn ,d;
  centro.cargue(Centrox ,Centroy ,Grano.Getz()) ; // Centro del obstaculo
  r21 = centro - Grano.r ; // Distancia al centro del obstaculo ( vector )
  d = norma(r21) ; // distancia al centro del obstaculo
  n = r21/d ; // vector normal
  if( d < R+R0 ) {
    Vc = Grano.V ; // Velocidad de contacto
    normaVcn = Vc*n ; // magnitud velocidad normal - proyeccion en el v .
    Vcn = n*normaVcn ; // Velocidad normal
    Grano.V =(Vc-2*Vcn);
  }
}
