/*
Clase Autómata celular que simula el suelo
*/
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
# include <stdlib.h>
using namespace std ;

//------------Autómata celular------------
void Automata::IniciaArreglo(void){
  // Inicia el suelo plano a altura R*número granos en y
  for(int i=0; i<Mx; i++)
    for(int k=0; k<Mz; k++)
      h[i][k]= hnew[i][k]=( suelo*ny );
}
void Automata::AgregaGrano(int i, int k) {
  h[i][k]+=1;
}
void Automata::QuitaGrano(int i, int k){
  h[i][k] -=1;
}
void Automata::Derrumbe(Crandom & ran64){
  int i, k, dummyindex; int aux[Mx][Mz], tocado[Mx][Mz]={0};
  double z[4];
  double prob[4], dummy;

  for(i=0;i<Mx;i++)
    for(k=0;k<Mz;k++)
      aux[i][k]=h[i][k];//Copia h a aux
  
  for (i=1; i<Mx-1; i++) {
    for (k=1; k<Mz-1; k++) {
      z[0] = (h[i][k] - h[i][k-1]); // pendiente arriba
      z[1] = (h[i][k] - h[i-1][k]); // pendiente izquierda
      z[2] = (h[i][k] - h[i][k+1]); // pendiente abajo
      z[3] = (h[i][k] - h[i+1][k]); // pendiente derecha

      //Solo revisar si yo me voy a derrumbar
      for(int j=0;j<4;j++){
	if (z[j]>z_thr){prob[j]=ran64.r();}//pendiente indebida positiva
	else {prob[j]=-1.0;}
      }

      dummy=0.0;
      dummyindex=-1;//Si nadie se puede derrumbar no hace nada
      for(int j=0;j<4;j++){if(dummy<prob[j]){dummy= prob[j];dummyindex=j;}}//dummyindex tiene el índice del que tenga el número al azar más alto o -1 si nadie se derrumbó

      //Calculo derrumbe
      if(dummyindex !=-1){
	//Se busca la dirección en la que se va a derrumbar de acuerdo a dummyindex
	//Al encontrar se quita en la posición actual, se suma en la dirección encontrada
	if(dummyindex == 0){
	  aux[i][k]--;
	  aux[i][k-1]++;
	  hnew[i][k]= aux[i][k];
	  hnew[i][k-1]=aux[i][k-1];
	  tocado[i][k-1]=1;
	}
	else if(dummyindex == 1){
	  aux[i][k]--;
	  aux[i-1][k]++;
	  hnew[i][k]= aux[i][k];
	  hnew[i-1][k]=aux[i-1][k];
	  tocado[i-1][k]=1;
	}
	else if(dummyindex == 2){
	  aux[i][k]--;
	  aux[i][k+1]++;
	  hnew[i][k]= aux[i][k];
	  hnew[i][k+1]=aux[i][k+1];
	  tocado[i][k+1]=1;
	}
	else if(dummyindex == 3){
	  aux[i][k]--;
	  aux[i+1][k]++;
	  hnew[i][k]= aux[i][k];
	  hnew[i+1][k]=aux[i+1][k];
	  tocado[i+1][k]=1;
	}
	hnew[i][k]= aux[i][k];
      }
      
      else if(tocado[i][k]!=0){}
      else hnew[i][k]=h[i][k];
    }
  }
  
  // Intercambio los arreglos
  for(i=0;i<Mx; i++)
    for(k=0;k<Mz; k++)h[i][k] = hnew[i][k] ;
}
void Automata::Erosiona(vector<Cuerpo> & vectorCuerpo , LatticeBoltzmann & Agua){
  vector3D Vf;
  double norma2Vf;
  int i ,z;
  int altura;
  double x_g , y_g , z_g; // Posicion del grano
  for(i=3*nx; i<(Mx-1); i++){// Para todas las columnas del Suelo ( suelo )
    for(z=nz; z<(Mz-1); z++){
      altura = Get_h(i,z); // altura de la columna
      x_g = (1.0*i)/(1.0*nx);//(0,1)->0, (2,3)->1, (4,5)->2, etc
      z_g = (1.0*z)/(1.0*nz);
      y_g = (1.0*altura + 1.0)/(1.0*ny);
      Vf = Agua.InterpolaVelocidad(x_g,y_g,z_g); // Velocidad
      norma2Vf = norma2(Vf) ;
      if ( (norma2Vf > V_Umbral) && (altura>0)) {//si V es mayor a velocidad umbral
	QuitaGrano(i, z) ;//Quita del suelo
	Cuerpo cuerpo(x_g ,y_g +0.5,z_g, 0.0, VySalida, 0.0, m0,R0) ;
	vectorCuerpo.push_back(cuerpo);
      }
    }
  }
}
void Automata::Imprimirh(const char * NombreArchivo){
  ofstream Archivo(NombreArchivo);
  Archivo<<"# vtk DataFile Version 1.0"<<endl;
  Archivo<<"3D triangulation data"<<endl;
  Archivo<<"ASCII"<<endl;
  Archivo<<"DATASET POLYDATA"<<endl;
  Archivo<<"POINTS "<< (Mx-1)*(Mz-1)<<" float"<<endl;
  for(int i=0;i<Mx;i++)
    for (int k=0; k<Mz; k++)Archivo <<i*1.0/nx<<" "<<k*1.0/nz<<" "<<h[i][k]*1.0/ny<<endl;
  Archivo.close();
}
void Automata::Imprimirhz10(const char * NombreArchivo){
  ofstream Archivo(NombreArchivo);
  for(int i=0;i<(Mx-1);i++){Archivo <<i*1.0/nx<<" "<<h[i][10]*1.0/ny<<endl;}
  Archivo.close();
}
