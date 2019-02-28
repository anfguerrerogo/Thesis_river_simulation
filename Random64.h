/*
Clase generador aleatorio
*/
#ifndef RANDOM
#define RANDOM
#include<iostream>
#include<cmath>
using namespace std;

//Constantes del generador aleatorio
class Crandom{
  unsigned long long u,v,w;
public:
  Crandom(unsigned long long j);//semilla que entra a un constructor (mismo nombre que clase)
  unsigned long long int64();
  double r() {return 5.42101086242752217E-20 * int64();}
  unsigned int int32(){return (unsigned int) int64();};
  double exponencial(float tau);
  double gauss(float mu,float sigma);
};
#endif
