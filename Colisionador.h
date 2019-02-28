/*
Archivo de cabecera para clase que calcula fuerzas.
*/
#ifndef COLISIONADOR_H
#define COLISIONADOR_H
# include "Constantes.h"
# include "Vector.h"
using namespace std ;
//------------------COLISIONADOR---------------
class Cuerpo;
class LatticeBoltzmann;
class Colisionador{
private:
public:
  void CalculeFuerzasDelFluido(vector<Cuerpo> &vectorCuerpo ,LatticeBoltzmann & Agua );
  void AgregaFuerzaObstaculo(Cuerpo & Grano);
  friend class Cuerpo;
  friend class LatticeBoltzmann;
};
#endif
