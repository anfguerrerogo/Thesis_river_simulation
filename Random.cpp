/*
Clase generador aleatorio
 */
#include<iostream>
#include<cmath>
# include "Random64.h"
using namespace std;
Crandom::Crandom(unsigned long long j){
    v=4101842887655102017LL; w=1;
    u = j ^ v; int64();
    v = u; int64();
    w = v; int64();
  }
unsigned long long Crandom::int64() {
    u = u * 2862933555777941757LL + 7046029254386353087LL;
    v ^= v >> 17; v ^= v << 31; v ^= v >> 8;//ciclos xor corridos 17 31 y 8 veces
    w = 4294957665U*(w & 0xffffffff) + (w >> 32);
    unsigned long long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
    return (x + v) ^ w;
  }

double Crandom::exponencial(float tau){
  return -tau*log(r());
}
double Crandom:: gauss(float mu,float sigma){
  return sigma*sqrt(-2*log(r()))*cos(2*M_PI*r())+mu;
}
