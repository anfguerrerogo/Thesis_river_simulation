/*
Clase Lattice Boltzmann D3Q19
 */
# include <iostream>
# include <fstream>
# include <cmath>
# include <vector>
# include "Vector.h"
# include "Random64.h"
# include "Automata3D.h"
# include "LB3D.h"
# include <math.h>
using namespace std ;
//---------------LATTICE BOLTZMANN-----------
LatticeBoltzmann::LatticeBoltzmann(void){ // Constructor
  //Pone el tamaño de V[x,y,z][i], w[i], EsFluido[ix][iy][iz], f[ix][iy][iz] y fnew[ix][iy][iz]
  V.resize(3);
  for (int i = 0; i < 3; ++i)V[i].resize(19);
  w.resize(19);
  EsFluido.resize(Lx);
  f.resize(Lx);
  fnew.resize(Lx);
  for (int ix = 0; ix < Lx; ++ix){
    EsFluido[ix].resize(Ly);
    f[ix].resize(Ly);
    fnew[ix].resize(Ly);
    for (int iy = 0; iy < Ly; ++iy){
      EsFluido[ix][iy].resize(Lz);
      f[ix][iy].resize(Lz);
      fnew[ix][iy].resize(Lz);
      for(int iz = 0; iz < Lz; ++iz){
	f[ix][iy][iz].resize(19);
	fnew[ix][iy][iz].resize(19);
      }
    }
  }

  // Carga los pesos D3Q19
  w[0]=1.0/3.0;//centro
  w[1]=w[2]=w[3]=w[4]=w[4]=w[5]=w[6]=1.0/18.0;//primeros vecinos
  w[7]=w[8]=w[9]=w[10]=w[11]=w[12]=w[13]=w[14]=w[15]=w[16]=w[17]=w[18]=1.0/36.0;//segundos vecinos
  
  // Carga los vectores D3Q19
  //Centro
  V[0][0]=0;
  V[1][0]=0;
  V[2][0]=0;
  //Right-------Up--------Left--------Down--------Front------Back----
  V[0][1]=1; V[0][2]=0; V[0][3]=-1; V[0][4]=0;  V[0][5]=0; V[0][6]=0;
  V[1][1]=0; V[1][2]=1; V[1][3]=0;  V[1][4]=-1; V[1][5]=0; V[1][6]=0;
  V[2][1]=0; V[2][2]=0; V[2][3]=0;  V[2][4]=0;  V[2][5]=1; V[2][6]=-1;
  
  //Right,up---Right,down---Left,down------Left,up
  V[0][7]=1;   V[0][8]=1;   V[0][9]=-1;  V[0][10]=-1; 
  V[1][7]=1;   V[1][8]=-1;  V[1][9]=-1;  V[1][10]=1;
  V[2][7]=0;   V[2][8]=0;   V[2][9]=0;   V[2][10]=0;

  //Up,front-----Up,back-----Down,back----Down,front
  V[0][11]=0;  V[0][12]=0;  V[0][13]=0;   V[0][14]=0;
  V[1][11]=1;  V[1][12]=1;  V[1][13]=-1;  V[1][14]=-1;
  V[2][11]=1;  V[2][12]=-1; V[2][13]=-1;  V[2][14]=1;

  //Right,front---Right,back-----Left,front----Left,back
  V[0][15]=1;     V[0][16]=1;   V[0][17]=-1;  V[0][18]=-1;
  V[1][15]=0;     V[1][16]=0;   V[1][17]=0;   V[1][18]=0;
  V[2][15]=1;     V[2][16]=-1;  V[2][17]=-1;  V[2][18]=1;
}
void LatticeBoltzmann::Inicie(void){
  //Inicializa f, fnew y Vel
  int ix, iy, iz, i;
  for(ix=0; ix<Lx; ix++)
    for(iy=0; iy<Ly; iy++)
      for(iz=0; iz<Lz; iz++){
	for (i=0; i<19; i++)fnew[ix][iy][iz][i]=f[ix][iy][iz][i]= feq(RHO0, UX0, UY0, UZ0, i);
	rho0=rho(ix,iy,iz);
	Ux0=Ux(ix,iy,iz);
	Uy0=Uy(ix,iy,iz);
	Uz0=Uz(ix,iy,iz);
	ImponerCampos(ix,iy,iz,rho0,Ux0,Uy0,Uz0);
  	Vel[ix][iy][iz][0]=Ux0;
	Vel[ix][iy][iz][1]=Uy0;
	Vel[ix][iy][iz][2]=Uz0;
      }  
}
void LatticeBoltzmann::IniciedeArchivo(const char * NombreArchivo){
  ifstream Archivo(NombreArchivo);
  int ix, iy, iz, i;
  string number;
  double data;
  for(ix=0; ix<Lx; ix++)
    for(iy=0; iy<Ly; iy++)
      for(iz=0; iz<Lz; iz++)
	for (i=0; i<19; i++){
	  if(!Archivo.eof()){
	    getline(Archivo,number,' '); //read number
	    data = atof(number.c_str()); //convert to double
	    fnew[ix][iy][iz][i]=f[ix][iy][iz][i]=data;}
	  else{cout<<"bad data"<<endl;}
	}
  Archivo.close();
}
void LatticeBoltzmann::VerificaCelda(Automata & Suelo){
  //Verifica cada celda del dominio si pertenece o no al fluido
  int ix, iy, iz, i, z;
  int altura_hip, altura_celda_i, altura_celda_i_right, altura_celda_i_down, altura_celda_i_right_down;
  int NGranosCelda;
  int nlim=0.5*(nx*ny*nz);//Mitad del número máximo posible en celda
  for(ix=0; ix<Lx; ix++){
    for(iy=0; iy<Ly; iy++){
      for(iz=0; iz<Lz; iz++){
	i = ix*nx; // Posicion en el arreglo
	z = iz*nz;
	altura_hip = iy*ny;
	altura_celda_i = Suelo.Get_h(i,z) - altura_hip;
	altura_celda_i_right = Suelo.Get_h(i+1,z) - altura_hip;
	altura_celda_i_down = Suelo.Get_h(i,z+1) - altura_hip;
	altura_celda_i_right_down = Suelo.Get_h(i+1,z+1) - altura_hip;
	NGranosCelda = altura_celda_i + altura_celda_i_right + altura_celda_i_down + altura_celda_i_right_down;//Número granos en celda
        //Si hay mas de la mitad de granos no es fluido
	if(NGranosCelda>=nlim)EsFluido[ix][iy][iz]=false;
	else EsFluido[ix][iy][iz]=true;
      }
    }
  }
}
void LatticeBoltzmann::Colisione(void){
  //Paso de Colisión
  int ix,iy,iz,i;
  double rho0,Ux0,Uy0, Uz0;
  for(ix=0; ix < Lx; ix++){
    for(iy=0; iy < Ly ; iy++){
      for(iz=0; iz < Lz ; iz++){
	rho0= rho(ix, iy, iz);
	Ux0 = Ux(ix , iy, iz);
	Uy0 = Uy(ix , iy, iz);
	Uz0 = Uz(ix , iy, iz);
	ImponerCampos(ix, iy, iz,rho0, Ux0, Uy0, Uz0);
	for(i=0; i<19; i++)
	  fnew[ix][iy][iz][i]=f[ix][iy][iz][i]-1/Tau*(f[ix][iy][iz][i] -feq(rho0, Ux0, Uy0,Uz0, i));
      }
    }
  }
}
void LatticeBoltzmann::Adveccione(void){
  //Paso de Advección
  int ix, iy, iz, i;
  for(ix=0; ix < Lx ; ix++){
    for(iy=0; iy < Ly ; iy++){
      for(iz=0; iz < Lz ; iz++){
	if(iz == 0){//abiertas hacia el fondo
	  f[ix][iy][iz][0]=fnew[ix][iy][iz+1][0];
	  f[ix][iy][iz][1]=fnew[(ix-1+Lx)%Lx][iy][iz+1][1];//r
	  f[ix][iy][iz][2]=fnew[ix][(iy-1+Ly)%Ly][iz+1][2];//u
	  f[ix][iy][iz][3]=fnew[(ix+1)%Lx][iy][iz+1][3];//l
	  f[ix][iy][iz][4]=fnew[ix][(iy+1)%Ly][iz+1][4];//d
	  f[ix][iy][iz][5]=fnew[ix][iy][iz+1][5];//f
	  f[ix][iy][iz][6]=fnew[ix][iy][iz+1][6];//b
	  f[ix][iy][iz][7]=fnew[(ix-1+Lx)%Lx][(iy-1+Ly)%Ly][iz+1][7];//r,u
	  f[ix][iy][iz][8]=fnew[(ix-1+Lx)%Lx][(iy+1)%Ly][iz+1][8];//r,d
	  f[ix][iy][iz][9]=fnew[(ix+1)%Lx][(iy+1)%Ly][iz+1][9];//l,d
	  f[ix][iy][iz][10]=fnew[(ix+1)%Lx][(iy-1+Ly)%Ly][iz+1][10];//l,u
	  f[ix][iy][iz][11]=fnew[ix][(iy-1+Ly)%Ly][iz+1][11];//u,f
	  f[ix][iy][iz][12]=fnew[ix][(iy-1+Ly)%Ly][iz+1][12];//u,b
	  f[ix][iy][iz][13]=fnew[ix][(iy+1)%Ly][iz+1][13];//d,b
	  f[ix][iy][iz][14]=fnew[ix][(iy+1)%Ly][iz+1][14];//d,f
	  f[ix][iy][iz][15]=fnew[(ix-1+Lx)%Lx][iy][iz+1][15];//r,f
	  f[ix][iy][iz][16]=fnew[(ix-1+Lx)%Lx][iy][iz+1][16];//r,b
	  f[ix][iy][iz][17]=fnew[(ix+1)%Lx][iy][iz+1][17];//l,b
	  f[ix][iy][iz][18]=fnew[(ix+1)%Lx][iy][iz+1][18];//l,f*/
	}
	else if(iy == Ly-1){//Condiciones de frontera abiertas arriba
	  f[ix][iy][iz][0]=fnew[ix][iy-1][iz][0];
	  f[ix][iy][iz][1]=fnew[(ix-1+Lx)%Lx][iy-1][iz][1];//r
	  f[ix][iy][iz][2]=fnew[ix][iy-1][iz][2];//u
	  f[ix][iy][iz][3]=fnew[(ix+1)%Lx][iy-1][iz][3];//l  El error está en el ix+1!!!!!!
	  f[ix][iy][iz][4]=fnew[ix][iy-1][iz][4];//d
	  f[ix][iy][iz][5]=fnew[ix][iy-1][(iz-1+Lz)%Lz][5];//f
	  f[ix][iy][iz][6]=fnew[ix][iy-1][(iz+1)%Lz][6];//b
	  f[ix][iy][iz][7]=fnew[(ix-1+Lx)%Lx][iy-1][iz][7];//r,u
	  f[ix][iy][iz][8]=fnew[(ix-1+Lx)%Lx][iy-1][iz][8];//r,d
	  f[ix][iy][iz][9]=fnew[(ix+1)%Lx][iy-1][iz][9];//l,d
	  f[ix][iy][iz][10]=fnew[(ix+1)%Lx][iy-1][iz][10];//l,u
	  f[ix][iy][iz][11]=fnew[ix][iy-1][(iz-1+Lz)%Lz][11];//u,f
	  f[ix][iy][iz][12]=fnew[ix][iy-1][(iz+1)%Lz][12];//u,b
	  f[ix][iy][iz][13]=fnew[ix][iy-1][(iz+1)%Lz][13];//d,b
	  f[ix][iy][iz][14]=fnew[ix][iy-1][(iz-1+Lz)%Lz][14];//d,f
	  f[ix][iy][iz][15]=fnew[(ix-1+Lx)%Lx][iy-1][(iz-1+Lz)%Lz][15];//r,f
	  f[ix][iy][iz][16]=fnew[(ix-1+Lx)%Lx][iy-1][(iz+1)%Lz][16];//r,b
	  f[ix][iy][iz][17]=fnew[(ix+1)%Lx][iy-1][(iz+1)%Lz][17];//l,b
	  f[ix][iy][iz][18]=fnew[(ix+1)%Lx][iy-1][(iz-1+Lz)%Lz][18];//l,f
	  }
	else if(ix == Lx-1){ //abiertas a la derecha
	  f[ix][iy][iz][0]=fnew[ix-1][iy][iz][0];
	  f[ix][iy][iz][1]=fnew[ix-1][iy][iz][1];//r
	  f[ix][iy][iz][2]=fnew[ix-1][(iy-1+Ly)%Ly][iz][2];//u
	  f[ix][iy][iz][3]=fnew[ix-1][iy][iz][3];//l
	  f[ix][iy][iz][4]=fnew[ix-1][(iy+1)%Ly][iz][4];//d
	  f[ix][iy][iz][5]=fnew[ix-1][iy][(iz-1+Lz)%Lz][5];//f
	  f[ix][iy][iz][6]=fnew[ix-1][iy][(iz+1)%Lz][6];//b
	  f[ix][iy][iz][7]=fnew[ix-1][(iy-1+Ly)%Ly][iz][7];//r,u
	  f[ix][iy][iz][8]=fnew[ix-1][(iy+1)%Ly][iz][8];//r,d
	  f[ix][iy][iz][9]=fnew[ix-1][(iy+1)%Ly][iz][9];//l,d
	  f[ix][iy][iz][10]=fnew[ix-1][(iy-1+Ly)%Ly][iz][10];//l,u
	  f[ix][iy][iz][11]=fnew[ix-1][(iy-1+Ly)%Ly][(iz-1+Lz)%Lz][11];//u,f
	  f[ix][iy][iz][12]=fnew[ix-1][(iy-1+Ly)%Ly][(iz+1)%Lz][12];//u,b
	  f[ix][iy][iz][13]=fnew[ix-1][(iy+1)%Ly][(iz+1)%Lz][13];//d,b
	  f[ix][iy][iz][14]=fnew[ix-1][(iy+1)%Ly][(iz-1+Lz)%Lz][14];//d,f
	  f[ix][iy][iz][15]=fnew[ix-1][iy][(iz-1+Lz)%Lz][15];//r,f
	  f[ix][iy][iz][16]=fnew[ix-1][iy][(iz+1)%Lz][16];//r,b
	  f[ix][iy][iz][17]=fnew[ix-1][iy][(iz+1)%Lz][17];//l,b
	  f[ix][iy][iz][18]=fnew[ix-1][iy][(iz-1+Lz)%Lz][18];//l,f
	}
	else if(iz == Lz-1){//abiertas hacia el inicio
	  f[ix][iy][iz][0]=fnew[ix][iy][iz-1][0];
	  f[ix][iy][iz][1]=fnew[(ix-1+Lx)%Lx][iy][iz-1][1];//r
	  f[ix][iy][iz][2]=fnew[ix][(iy-1+Ly)%Ly][iz-1][2];//u
	  f[ix][iy][iz][3]=fnew[(ix+1)%Lx][iy][iz-1][3];//l
	  f[ix][iy][iz][4]=fnew[ix][(iy+1)%Ly][iz-1][4];//d
	  f[ix][iy][iz][5]=fnew[ix][iy][iz-1][5];//f
	  f[ix][iy][iz][6]=fnew[ix][iy][iz-1][6];//b
	  f[ix][iy][iz][7]=fnew[(ix-1+Lx)%Lx][(iy-1+Ly)%Ly][iz-1][7];//r,u
	  f[ix][iy][iz][8]=fnew[(ix-1+Lx)%Lx][(iy+1)%Ly][iz-1][8];//r,d
	  f[ix][iy][iz][9]=fnew[(ix+1)%Lx][(iy+1)%Ly][iz-1][9];//l,d
	  f[ix][iy][iz][10]=fnew[(ix+1)%Lx][(iy-1+Ly)%Ly][iz-1][10];//l,u
	  f[ix][iy][iz][11]=fnew[ix][(iy-1+Ly)%Ly][iz-1][11];//u,f
	  f[ix][iy][iz][12]=fnew[ix][(iy-1+Ly)%Ly][iz-1][12];//u,b
	  f[ix][iy][iz][13]=fnew[ix][(iy+1)%Ly][iz-1][13];//d,b
	  f[ix][iy][iz][14]=fnew[ix][(iy+1)%Ly][iz-1][14];//d,f
	  f[ix][iy][iz][15]=fnew[(ix-1+Lx)%Lx][iy][iz-1][15];//r,f
	  f[ix][iy][iz][16]=fnew[(ix-1+Lx)%Lx][iy][iz-1][16];//r,b
	  f[ix][iy][iz][17]=fnew[(ix+1)%Lx][iy][iz-1][17];//l,b
	  f[ix][iy][iz][18]=fnew[(ix+1)%Lx][iy][iz-1][18];//l,f
}
	else{for(i=0; i<19; i++) f[(ix+V[0][i]+ Lx)%Lx][(iy+V[1][i]+ Ly)%Ly][(iz+V[2][i]+ Lz)%Lz][i]=fnew[ix][iy][iz][i];}
      }
    }
  }
}
double LatticeBoltzmann::rho(int ix, int iy, int iz){
  //Devuelve la densidad en cierto punto
  int i; double suma;
  for(suma =0, i =0; i < 19; i++)
    suma += f[ix][iy][iz][i];
  return suma;
}
double LatticeBoltzmann::Ux(int ix, int iy, int iz){
  //Devuelve  la velocidad en x en cierto punto
  int i; double suma;
  for(suma = 0, i = 0; i < 19; i++)
    suma += V[0][i]*f[ix][iy][iz][i];
  return suma/rho(ix,iy, iz) ;
}
double LatticeBoltzmann::Uy(int ix, int iy, int iz){
  //Devuelve  la velocidad en y en cierto punto
  int i; double suma;
  for(suma=0,i=0;i<19; i++)
    suma += V[1][i]*f[ix][iy][iz][i];
  return suma/rho(ix, iy, iz);
}
double LatticeBoltzmann::Uz(int ix, int iy, int iz){
  //Devuelve  la velocidad en z en cierto punto
  int i; double suma;
  for(suma=0,i=0;i<19; i++)
    suma += V[2][i]*f[ix][iy][iz][i];
  return suma/rho(ix, iy, iz);
}
void LatticeBoltzmann::ImponerCampos(int ix,int iy,int iz, double & rho0,double & Ux0,double & Uy0, double & Uz0){
  //Impone valores de velocidad en ciertas regiones (fuentes, partes que no son fluido y obstáculos)
  if(ix <= 3 ){Ux0=UX1;Uy0=UY1;Uz0=UZ1;}//fuente
  if(EsFluido[ix][iy][iz]==false){Ux0=Uy0=Uz0=0;}//Si no es fluido no haga nada
  if((ix - Centrox )*(ix - Centrox )+(iy - Centroy )*(iy - Centroy ) <= R*R){Ux0=Uy0=Uz0=0;}//obstáculo (cilindro 1)
}
double LatticeBoltzmann::feq(double rho0,double Ux0,double Uy0,double Uz0, int i){
  //Calcula función equilibrio
  double U2=Ux0*Ux0+Uy0*Uy0+Uz0*Uz0;
  double UpVi =Ux0*V[0][i]+Uy0*V[1][i]+Uz0*V[2][i];
  return rho0*w[i]*(1+4.5*UpVi*UpVi+3*UpVi-1.5*U2);
}
void LatticeBoltzmann::DeltaVelocidad(int t, ofstream &stream){
  double rho0,Ux0,Uy0,Uz0, delta, total;
  double SumaDelta=0.0;
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      for(int iz=0;iz<Lz;iz++){
	rho0=rho(ix,iy,iz);
	Ux0=Ux(ix,iy,iz);
	Uy0=Uy(ix,iy,iz);
	Uz0=Uz(ix,iy,iz);
	ImponerCampos(ix,iy,10,rho0,Ux0,Uy0,Uz0);
	delta=sqrt(pow(Vel[ix][iy][iz][0]-Ux0,2)+pow(Vel[ix][iy][iz][1]- Uy0,2)+pow(Vel[ix][iy][iz][2]-Uz0,2));//dU=U(t-1)-U(t)
	SumaDelta+=delta;//Sumatoria de dU
	Vel[ix][iy][iz][0]=Ux0;//Actualiza t
	Vel[ix][iy][iz][1]=Uy0;
	Vel[ix][iy][iz][2]=Uz0;
      }
    }
  }
  total=SumaDelta/(Lx*Ly*Lz);
  stream<<t<<" "<<total<<" "<<log10(total)<<endl;
}

vector3D LatticeBoltzmann::Vf(int ix,int iy, int iz){
  //Devuelve vector velocidad en cierto punto
  vector3D Vf0;
  double Ux0, Uy0, Uz0;
  Ux0= Ux(ix,iy, iz);
  Uy0=Uy(ix,iy, iz);
  Uz0=Uz(ix,iy,iz);
  Vf0.cargue(Ux0,Uy0,Uz0);
  return Vf0;
}
vector3D LatticeBoltzmann::InterpolaVelocidad(double x,double y, double z){
  //Devuelve vector velocidad para coordenadas en unidades del autómata usando interpolación trilineal
  vector3D Vint;
  double xd, yd, zd,xr, yr, zr;
  int x0, y0, z0, x1, y1, z1; // Coordenadas de los vertices del cubo
  double intx, inty, intz;// Parte entera de la coordenada
  double ux00, uy00, uz00, ux10, uy10, uz10, ux01, uy01, uz01, ux11, uy11, uz11, ux0, uy0, uz0, ux1, uy1, uz1;
  double ux, uy, uz;
  
  xr=x,yr=y, zr=z;// Posición donde quiero hallar la velocidad
  xd = modf(xr,&intx) ; yd = modf(yr,&inty); zd = modf(zr,&intz); //Escalar entre 0 y 1 usando parte decimal
  
  //Definir cubo de lado 1
  x0 = xr; y0= yr; z0=zr;
  //Si se va a salir de bounds no lo haga
  if(x0>=Lx-1){
    ux=Ux(Lx-1,y0,z0);
    uy=Uy(Lx-1,y0,z0);
    uz=Uz(Lx-1,y0,z0);
    Vint.cargue(ux,uy,uz);
    return Vint;}
  else if(y0>=Ly-1){
    ux=Ux(x0,Ly-1,z0);
    uy=Uy(x0,Ly-1,z0);
    uz=Uz(x0,Ly-1,z0);
    Vint.cargue(ux,uy,uz);
    return Vint;}
  else if(z0>= Lz-1){
    ux=Ux(x0,y0,Lz-1);
    uy=Uy(x0,y0,Lz-1);
    uz=Uz(x0,y0,Lz-1);
    Vint.cargue(ux,uy,uz);
    return Vint;}
  x1 = x0+1; y1 = y0+1; z1= z0+1;
  
  //Interpolar en x
  ux00=Ux(x0,y0,z0)*(1-xd)+Ux(x1,y0,z0)*xd;
  uy00=Uy(x0,y0,z0)*(1-xd)+Uy(x1,y0,z0)*xd;
  uz00=Uz(x0,y0,z0)*(1-xd)+Uz(x1,y0,z0)*xd;
  ux01=Ux(x0,y0,z1)*(1-xd)+Ux(x1,y0,z1)*xd;
  uy01=Uy(x0,y0,z1)*(1-xd)+Uy(x1,y0,z1)*xd;
  uz01=Uz(x0,y0,z1)*(1-xd)+Uz(x1,y0,z1)*xd;  
  ux10=Ux(x0,y1,z0)*(1-xd)+Ux(x1,y1,z0)*xd;
  uy10=Uy(x0,y1,z0)*(1-xd)+Uy(x1,y1,z0)*xd;
  uz10=Uz(x0,y1,z0)*(1-xd)+Uz(x1,y1,z0)*xd;
  ux11=Ux(x0,y1,z1)*(1-xd)+Ux(x1,y1,z1)*xd;
  uy11=Uy(x0,y1,z1)*(1-xd)+Uy(x1,y1,z1)*xd;
  uz11=Uz(x0,y1,z1)*(1-xd)+Uz(x1,y1,z1)*xd;
  
  //Interpolar en y
  ux0 =ux00*(1-yd)+ux10*(yd);
  uy0 =uy00*(1-yd)+uy10*(yd);
  uz0 =uz00*(1-yd)+uz10*(yd);
  ux1 =ux01*(1-yd)+ux11*(yd);
  uy1 =uy01*(1-yd)+uy11*(yd);
  uz1 =uz01*(1-yd)+uz11*(yd);
  
  //Interpolar en z
  ux = ux0*(1-zd)+ux1*zd;
  uy = uy0*(1-zd)+uy1*zd;
  uz = uz0*(1-zd)+uz1*zd;
  
  //Cargar
  Vint.cargue(ux,uy,uz);
  return Vint;
}
//--------Funciones de imprimir--------
void LatticeBoltzmann::Imprimirz10(const char * NombreArchivo){
  //Imprime en archivo externo para un z específico los valores
  // x y Vx Vy densidad
  ofstream Archivo(NombreArchivo);
  double rho0,Ux0,Uy0,Uz0;
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,10);
      Ux0=Ux(ix,iy,10);
      Uy0=Uy(ix,iy,10);
      Uz0=Uz(ix,iy,10);
      ImponerCampos(ix,iy,10,rho0,Ux0,Uy0,Uz0);
      Archivo<<ix<<" "<<iy<<" "<<Ux0<<" "<<Uy0<<" "<<rho0<<endl;
    }
    Archivo<<endl;
  }
  Archivo.close();
}
void LatticeBoltzmann::Imprimirz10velocidadx(const char * NombreArchivo){
  //Imprime en archivo externo para un z específico la velocidad en x en cada punto
  ofstream Archivo(NombreArchivo);
  double rho0,Ux0,Uy0,Uz0;
  //Orden inverso por cómo funciona streamplot en python
  for(int iy=0;iy<Ly;iy++){
    for(int ix=0;ix<Lx;ix++){
      rho0=rho(ix,iy,10);
      Ux0=Ux(ix,iy,10);
      Uy0=Uy(ix,iy,10);
      Uz0=Uz(ix,iy,10);
      ImponerCampos(ix,iy,10,rho0,Ux0,Uy0,Uz0);
      Archivo<<Ux0<<" ";
    }
    Archivo<<endl;
  }
  Archivo.close();
}
void LatticeBoltzmann::Imprimirz10velocidady(const char * NombreArchivo){
  //Imprime en archivo externo para un z específico la velocidad en y en cada punto
  ofstream Archivo(NombreArchivo);
  double rho0,Ux0,Uy0,Uz0;
  //Orden inverso por cómo funciona streamplot en python
  for(int iy=0;iy<Ly;iy++){
    for(int ix=0;ix<Lx;ix++){
      rho0=rho(ix,iy,10);
      Ux0=Ux(ix,iy,10);
      Uy0=Uy(ix,iy,10);
      Uz0=Uz(ix,iy,10);
      ImponerCampos(ix,iy,10,rho0,Ux0,Uy0,Uz0);
      Archivo<<Uy0<<" ";
    }
    Archivo<<endl;
  }
  Archivo.close();
}
void LatticeBoltzmann::ImprimirEsFluido(const char * NombreArchivo){
  //Imprime en un archivo externo: x y z 1_si_pertenece_a_fluido,_0_si_no
  ofstream Archivo(NombreArchivo);
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      for(int iz=0;iz<Lz;iz++){
	Archivo<<ix<<" "<<iz<<" "<<iy<<" "<<EsFluido[ix][iy][iz]<<endl;
      }
    }
  }
  Archivo.close();
}
void LatticeBoltzmann::Imprimir3D(const char * NombreArchivo){
  //Imprime en archivo externo: x z y Vx Vz Vy densidad
  ofstream Archivo(NombreArchivo);
  double rho0,Ux0,Uy0,Uz0;
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      for(int iz=0;iz<Lz;iz++){
	rho0=rho(ix,iy,iz);
	Ux0=Ux(ix,iy,iz);
	Uy0=Uy(ix,iy,iz);
	Uz0=Uz(ix,iy,iz);
	ImponerCampos(ix,iy,iz,rho0,Ux0,Uy0,Uz0);
	Archivo<<ix<<" "<<iz<<" "<<iy<<" "<<Ux0<<" "<<Uz0<<" "<<Uy0<<" "<<rho0<<endl;
      }
    }
    Archivo<<endl;
  }
  Archivo.close();
}
void LatticeBoltzmann::ImprimirVorticidad(const char * NombreArchivo){
  //Calcula la vorticidad y la imprime en archivo externo como: x y z Wx Wy Wz
  ofstream Archivo(NombreArchivo);
  double rho0,Ux0,Uy0,Uz0;
  double rhox, rhoy, rhoz;
  double Uxx,Uxy,Uxz,Uyx,Uyy,Uyz, Uzx, Uzy, Uzz;
  double wx, wy, wz;
  for(int ix=0;ix<Lx-1;ix++){
    for(int iy=0;iy<Ly-1;iy++){
      for(int iz=0;iz<Lz-1;iz++){
	rho0=rho(ix,iy,iz);
	Ux0=Ux(ix,iy,iz);
	Uy0=Uy(ix,iy,iz);
	Uz0=Uz(ix,iy,iz);
	ImponerCampos(ix,iy,iz,rho0,Ux0,Uy0,Uz0);
	
	//Para derivadas
	rhox=rho(ix+1,iy,iz);
	rhoy=rho(ix,iy+1,iz);
	rhoz=rho(ix,iy,iz+1);
	
	Uxx=Ux(ix+1,iy,iz);
	Uxy=Ux(ix,iy+1,iz);
	Uxz=Ux(ix,iy,iz+1);
	Uyx=Uy(ix+1,iy,iz);
	Uyy=Uy(ix,iy+1,iz);
	Uyz=Uy(ix,iy,iz+1);
	Uzx=Uz(ix+1,iy,iz);
	Uzy=Uz(ix,iy+1,iz);
	Uzz=Uz(ix,iy,iz+1);
	ImponerCampos(ix+1,iy,iz,rhox,Uxx,Uyx,Uzx);
	ImponerCampos(ix,iy+1,iz,rhoy,Uxy,Uyy,Uzy);
	ImponerCampos(ix,iy,iz+1,rhoz,Uxz,Uyz,Uzz);
	//Vorticidad
	wx=Uzy-Uz0-Uyz+Uy0;
	wy=Uxz-Ux0-Uzx+Uz0;
	wz=Uyx-Uy0-Uxy+Ux0;
	Archivo<<ix<<" "<<iz<<" "<<iy<<" "<<wx<<" "<<wz<<" "<<wy<<endl;
      }
    }
    Archivo<<endl;
  }
  Archivo.close();
}
void LatticeBoltzmann::Imprimirfunciones(const char * NombreArchivo){
  ofstream Archivo(NombreArchivo);
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      for(int iz=0;iz<Lz;iz++){
	for(int i=0;i<19;i++){
	  Archivo<<f[ix][iy][iz][i]<<" ";
	}
	Archivo<<endl;
      }
      Archivo<<endl;
    }
    Archivo<<endl;
  }
  Archivo.close();
}
