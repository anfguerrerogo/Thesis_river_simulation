/*
Archivo de cabecera para clase que maneja partículas
*/
#ifndef CUERPO_H
#define CUERPO_H
# include "Constantes.h"
# include "Vector.h"
using namespace std;
//-------------------CUERPO--------------------
class Colisionador;
class LatticeBoltzmann;
class Automata;
class Cuerpo {
private:  
  vector3D Rgauss;  
  vector3D Fg, Fb; // Fuerzas de gravedad y empuje
  vector3D Fd, Vf; // Drag Force ;
  vector3D r, V, F; // Propiedades vectoriales del cuerpo
  double m, R; // Propiedades escalares del cuerpo
  bool Existencia;
public:
  Cuerpo();
  Cuerpo(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0, double R0);// overload constructor
  ~ Cuerpo();
  void AgregaFuerzaArrastreyGravedad(LatticeBoltzmann & Agua);
  void Mueva_r1(double dt);
  void Mueva_V(LatticeBoltzmann & Agua, double dt);
  void Mueva_r2(double dt);
  void DesplazaGaussiana(Crandom & ran2, double mu ,double sigma);
  void RevisaSuelo(Automata & Suelo, LatticeBoltzmann & Agua);
  void RestaVelocidades(vector3D V0){V -= V0;};
  void SumaVelocidades(vector3D V0){V += V0;};
  void BorraFuerza(void){F.cargue(0 ,0 ,0);};
  void EliminaGrano(void){Existencia = false;};
  bool Get_Existencia(void){return Existencia;};
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  double Getz(void){return r.z();};
  friend class Colisionador;
  friend class LatticeBoltzmann;
  friend class Automata;};
#endif
