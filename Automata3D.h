/*
Archivo de cabecera para clase que maneja autómata celular
*/
#ifndef AUTOMATA3D_H
#define AUTOMATA3D_H
# include "Constantes.h"
class Cuerpo;
class LatticeBoltzmann;
class Automata {
private :
  int h[Mx][Mz],hnew[Mx][Mz];
public :
  void IniciaArreglo(void); 
  void AgregaGrano(int i , int k);
  void QuitaGrano(int i , int k);
  void Derrumbe(Crandom & ran64);
  void Erosiona(vector<Cuerpo> & vectorCuerpo, LatticeBoltzmann & Agua);
  void Imprimirh(const char * NombreArchivo);
  void Imprimirhz10(const char * NombreArchivo);	
  int Get_h(int i, int k) {return h[i][k];};
  friend class Cuerpo;
  friend class LatticeBoltzmann;
};

#endif
