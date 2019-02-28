/*
Archivo de cabecera para la clase que maneja el Lattice Boltzmann D3Q19
*/

#ifndef LB3D_H
#define LB3D_H
#include "Constantes.h"
#include <vector>
//------------LATTICE BOLTZMANN-----------
class Cuerpo;
class Automata;
class LatticeBoltzmann{
private :
//Funciones de propagación
  std::vector<vector<vector<vector<double> > > > f;
  std::vector<vector<vector<vector<double> > > > fnew;
//Vectores en cada dirección
  std::vector<vector<int> > V;
//Pesos
  std::vector<double> w;
//Variable booleano para saber si hace o no parte del fluido
  std::vector<vector<vector<bool> > > EsFluido;
//Vector auxiliar de velocidad en cada punto
  std::vector<vector<vector<vector<double> > > > Vel;
public :
  LatticeBoltzmann(void); // Constructor
  void Inicie(void);//Inicializa f y fnew en f equilibrio
  void IniciedeArchivo(const char * NombreArchivo);//Inicia Lattice a partir de un archivo de funciones externo
  void VerificaCelda(Automata & Suelo);  //Verifica cada celda del dominio si pertenece o no al fluido
  void Colisione(void);  //Paso de Colisión
  void Adveccione(void);  //Paso de Advección
  double rho(int ix ,int iy, int iz);//Devuelve la densidad en cierto punto
  double Ux(int ix ,int iy, int iz); //Devuelve  la velocidad en x en cierto punto
  double Uy(int ix ,int iy, int iz); //Devuelve  la velocidad en y en cierto punto
  double Uz(int ix ,int iy, int iz); //Devuelve  la velocidad en z en cierto punto
  void ImponerCampos(int ix,int iy,int iz, double & rho0,double & Ux0,double & Uy0, double & Uz0); //Impone valores de velocidad en ciertas regiones (fuentes, partes que no son fluido y obstáculos)
  double feq(double rho0, double Ux0, double Uy0,double Uz0, int i); //Calcula función equilibrio
  void DeltaVelocidad(int t, ofstream &stream); //Imprime el cambio dv entre la velocidad en dos pasos conescutivos a archivo externo
  vector3D Vf(int ix ,int iy, int iz); //Devuelve vector velocidad en cierto punto
  vector3D InterpolaVelocidad(double xd,double yd, double zd); //Devuelve vector velocidad para coordenadas en unidades del autómata usando interpolación trilineal
//--------Funciones de imprimir--------
  void Imprimirz10(const char * NombreArchivo); //Imprime para un z específico los valores: x y Vx Vy densidad
  void Imprimirz10velocidadx(const char * NombreArchivo); //Imprime en archivo externo para un z específico la velocidad en x en cada punto
  void Imprimirz10velocidady(const char * NombreArchivo); //Imprime en archivo externo para un z específico la velocidad en y en cada punto	
  void ImprimirEsFluido(const char * NombreArchivo); //Imprime en un archivo externo: x y z 1_si_pertenece_a_fluido,_0_si_no
  void Imprimir3D(const char * NombreArchivo);//Imprime en archivo externo: x z y Vx Vz Vy densidad
  void ImprimirVorticidad(const char * NombreArchivo);//Calcula la vorticidad y la imprime en archivo externo como: x y z Wx Wy Wz
  void Imprimirfunciones(const char * NombreArchivo); //Imprime funciones de propagación a archivo externo
  friend class Cuerpo;
  friend class Automata;};
#endif
