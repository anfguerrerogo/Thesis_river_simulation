/*
Programa para simular comportamiento tridimensional de un dominio rectangular en el que interactuan un fluido y un suelo de arena de altura cambiante en el tiempo para cada punto
Hecho por Andrés Felipe Guerrero González para tesis de pregrado, basado en el trabajo hecho por Juliana García Sarmiento para su tesis de pregrado.
Dirigido por Jose Daniel Muñoz
 */
# include <iostream>
# include <fstream>
# include <cmath>
# include <vector>
# include "Vector.h"
# include "Random64.h"
# include "Automata3D.h"
# include "Cuerpo.h"
# include "Colisionador.h"
# include "LB3D.h"
# include <math.h>
using namespace std;

//---------------------------FUNCIONES------------------
void Transporta(vector<Cuerpo> &vectorCuerpo,LatticeBoltzmann & Agua,Automata & Suelo,Colisionador & Fuerzas,Crandom & ran2, double Dt);
void ImprimirCuerpos(vector<Cuerpo> &vectorCuerpo, const char * NombreArchivo);
//-----------------------------MAIN---------------------
int main(){
  //Declaraciones iniciales
  int tau=1333;//tiempo característico LB
  LatticeBoltzmann Agua;
  Colisionador Fuerzas;
  Automata Suelo;
  Crandom ran2(8);
  vector <Cuerpo> vectorCuerpo;
  int t;

  //Inicialización de autómata celular y LB D3Q19
  Suelo.IniciaArreglo(); // Inicia el perfil del suelo
  Agua.VerificaCelda(Suelo) ; //Asigna EsFluido
  Agua.Inicie();//Inicia Lattice desde 0
  //Agua.IniciedeArchivo("funciones_Lx_200_tau_1333.dat");//Inicia Lattice desde archivo de funciones
  
  //Ciclo de estabilización
  for(t=0; t<3*tau ; t++){
    
    Agua.VerificaCelda(Suelo) ; //Asigna EsFluido
    Agua.Adveccione();//Advección LB
    Agua.Colisione(); //Colisión LB
    
    //--------------------Imprimir valores---------------
    if(t%1000==0){
      std::string a = std::to_string(t);
      //Imprimir campo de velocidades
      Agua.Imprimir3D((a + "_campo.vtk").c_str());     
      //Imprimir Nivel Suelo en 3D para Paraview
      Suelo.Imprimirh((a + "nivel_arena.vtk").c_str());
      //Imprimir campo de velocidades en x y y para un valor particular de z para hacer Streamlines en Python
      Agua.Imprimirz10velocidadx((a + "velocidadx.dat").c_str());
      Agua.Imprimirz10velocidady((a + "velocidady.dat").c_str());
      //Imprimir campo de velocidades bidimensional en un z específico
      Agua.Imprimirz10((a + "_corte.dat").c_str());
      //Imprimir altura del suelo a un z específico
      Suelo.Imprimirhz10((a+"altura.dat").c_str());
    }
  }
  
  //Imprimir las funciones del LB a un archivo externo para guardar la configuración y no tener que repetir la estabilización
  Agua.Imprimirfunciones("funciones_Lx_200_tau_1333.dat");
  
  //Ciclo de interacción entre LB y autómata
  for(t=0; t<10*tau ; t++){
    
    Agua.VerificaCelda(Suelo) ; //Asigna EsFluido
    Agua.Adveccione();//Advección LB
    Agua.Colisione(); //Colisión LB
    
    //---------------------Algoritmo Verlet--------------
    for(int t_auto=0; t_auto<2; t_auto++){Transporta(vectorCuerpo, Agua, Suelo, Fuerzas, ran2, Deltat);}
    for(int t_auto=0; t_auto<4; t_auto++)Suelo.Derrumbe(ran2);
    
    //--------------------Erosión del suelo--------------
    if(t%20==0)Suelo.Erosiona(vectorCuerpo ,Agua);
    
    //---------------------Imprimir valores--------------
    if(t%200==0){
      std::string s = std::to_string(t);
      //Imprimir campo de velocidades
      Agua.Imprimir3D(("Simulation/" + s + "_campo.vtk").c_str());     
      //Imprimir Nivel Suelo en 3D para Paraview
      Suelo.Imprimirh(("Simulation/" + s + "nivel_arena.vtk").c_str());
       //Imprimir campo de velocidades en x y y para un valor particular de z para hacer Streamlines en Python//Imprimir Vorticidad
      Agua.Imprimirz10velocidadx(("Simulation/" + s + "velocidadx.dat").c_str());
      Agua.Imprimirz10velocidady(("Simulation/" + s + "velocidady.dat").c_str());
      //Imprimir campo de velocidades bidimensional en un z específico
      Agua.Imprimirz10(("Simulation/" + s + "_corte.dat").c_str());
      //Imprimir altura del suelo a un z específico
      Suelo.Imprimirhz10(("Simulation/" + s + "altura.dat").c_str());
      //Imprimir partículas en el fluido (sedimentos)
      ImprimirCuerpos(vectorCuerpo, ("Simulation/" + s + "cuerpos.dat").c_str());
      //Imprimir Vorticidad del fluido
      Agua.ImprimirVorticidad(("Simulation/" + s + "vorticidad.dat").c_str());
      //Imprimir si hace o no parte del fluido cada punto del dominio
      Agua.EsFluido(("Simulation/" + s + "esfluido.dat").c_str());    
    }
  }
  return 0;
}


void ImprimirCuerpos(vector<Cuerpo> &vectorCuerpo, const char * NombreArchivo){
  //Imprime cada una de las partículas en vectorCuerpo en un formato para representar en Paraview
  ofstream Archivo(NombreArchivo);
  Archivo<<"# vtk DataFile Version 1.0"<<endl;
  Archivo<<"3D triangulation data"<<endl;
  Archivo<<"ASCII"<<endl;
  Archivo<<"DATASET POLYDATA"<<endl;
  Archivo<<"POINTS "<< vectorCuerpo.size()<<" float"<<endl;
  for(int i=0; (unsigned)i<vectorCuerpo.size()-1; i++){
    Archivo<<vectorCuerpo[i].Getx()<<" "<<vectorCuerpo[i].Getz()<<" "<<vectorCuerpo[i].Gety()<<endl;
    }
  Archivo.close();
}

void Transporta(vector<Cuerpo> &vectorCuerpo,LatticeBoltzmann & Agua,Automata & Suelo,Colisionador & Fuerzas,Crandom & ran2, double Dt){
  //Función que implementa algoritmo de movimiento de partículas de acuerdo a fuerzas del fluido y gravitacional
  
  //Revisa casillas que son o no fluido
  for(int i=0; (unsigned)i < vectorCuerpo.size(); i++){
    vectorCuerpo[i].RevisaSuelo(Suelo,Agua);
    if(vectorCuerpo[i].Get_Existencia() == false)vectorCuerpo.erase(vectorCuerpo.begin()+i);
  }
  //Mueve partículas
  for(int i=0; (unsigned)i<vectorCuerpo.size(); i++)vectorCuerpo[i].Mueva_r1(Dt);
  //Calcula fuerzas
  Fuerzas.CalculeFuerzasDelFluido(vectorCuerpo, Agua);
  //Mueve partículas
  for(int i=0; (unsigned)i<vectorCuerpo.size() ; i++){
    vectorCuerpo[i].Mueva_V(Agua ,Dt);//cambia V
    vectorCuerpo[i].Mueva_r2(Dt);
  }
  //Calcula Fuerzas
  Fuerzas.CalculeFuerzasDelFluido(vectorCuerpo, Agua);
  //Mueve partículas
  for(int i=0; (unsigned)i<vectorCuerpo.size() ; i++){
    vectorCuerpo[i].Mueva_V(Agua ,Dt);//cambia V
    vectorCuerpo[i].Mueva_r1(Dt);
  }
  //Mueve Gaussiana (movimiento aleatorio)
  for(int i=0; (unsigned)i < vectorCuerpo.size(); i++)vectorCuerpo[i].DesplazaGaussiana(ran2, mu, sigma);
  //Revisas si casillas que son o no fluido han cambiado
  for(int i=0; (unsigned)i < vectorCuerpo.size(); i++){
    vectorCuerpo[i].RevisaSuelo(Suelo,Agua);
    if(vectorCuerpo[i].Get_Existencia() == false)vectorCuerpo.erase(vectorCuerpo.begin()+i);
  }
}
