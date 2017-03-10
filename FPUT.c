#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include "omp.h"

#define N 64
#define B 1.0
#define dt 0.005
#define PI 3.14159265
#define T 5*pow(N,2.2)

double energia(double *x,double *v,int m);
double *F(double  *x);

int main(int argc, char **argv){
  
  int threads=atoi(argv[1]);
  omp_set_num_threads(threads);
  int Tiempos = T/dt;
  double *X = malloc(N*sizeof(double));
  double *V = malloc(N*sizeof(double));
  double *V_half = malloc(N*sizeof(double));
  double *aux_X = malloc(N*sizeof(double));
  double *aux_V = malloc(N*sizeof(double));
  double *E1 = malloc(1000*sizeof(double));
  double *E2 = malloc(1000*sizeof(double));
  double *E3 = malloc(1000*sizeof(double));
  
  int a;
  int i;
  int t;
  int j;
  int k;
  int z;
  int e=1;

  //condiciones iniciales
  X[0]=X[N-1]=0; 
  V[0]=V[N-1]=0;
  V_half[0]=V_half[N-1]=0;
  aux_X[0]=aux_X[N-1]=0;
  aux_V[0]=aux_V[N-1]=0;
  
  double ti= omp_get_wtime();

  FILE* Es;
  Es=fopen("datos.dat","wt");
  

#pragma omp parallel for shared(X,V)
  for(i=1;i<N-1;i++){
    X[i] = sin(PI*i/(N-1));
    V[i] = 0;
    V_half[i] = 0;
    
  }
  
  
  E1[0]=energia(X,V,1);
  E2[0]=energia(X,V,2);
  E3[0]=energia(X,V,3);  

  //imprime primeras energías
  fprintf(Es,"%f %f %f \n",E1[0],E2[0],E3[0]);

  //solución paralelizada
  for(t=1;t<Tiempos;t++){
    double *ac;
    double *ac2;
#pragma omp parallel for shared(X,V,aux_X,aux_V)
    for(z=1;z<N-1;z++){
      aux_V[z]=V[z];
      aux_X[z]=X[z];
    }
    ac=F(aux_X);
#pragma omp parallel for shared(X,V_half,aux_X,aux_V)
    for(j=1;j<N-1;j++){
      V_half[j] = aux_V[j] + 0.5*dt*ac[j];
    }
    
#pragma omp parallel for shared(X,V,V_half,aux_X,aux_V)
    for(a=1;a<N-1;a++){
      X[a]=aux_X[a]+dt*V_half[a];
    }
    ac2=F(X);
#pragma omp parallel for shared(X,V,V_half)
    for(k=1;k<N-1;k++){
      V[k] = V_half[k]+0.5*dt*ac2[k];
    }
    if(t%((int)Tiempos/1000)==0){
	E1[e]=energia(X,V,1);
	E2[e]=energia(X,V,2);
	E3[e]=energia(X,V,3);
	fprintf(Es,"%f %f %f \n",E1[e],E2[e],E3[e]);
	e+=1;
      }
  }

  fclose(Es);
  double tf= omp_get_wtime();
  double tt=tf-ti;
  FILE* TIEMPOS;
  TIEMPOS=fopen("tiempos.dat","wt");
  fprintf(TIEMPOS,"%f %d",tt,threads);
  fclose(TIEMPOS);

  return 0;
}



double energia(double *x,double *v,int m)
{
  double Qp = 0;
  double Qk = 0;
  double cte = pow(2.0/N,0.5);
  double w=4*pow(sin(m*PI/(2*N+2)),2);
  int l;
  int y;
  double E;
  #pragma omp parallel for reduction (+ : Qp)
  for(l=0;l<N;l++){
    Qp+= x[l]*sin(PI*m*l/N);
  }
  #pragma omp parallel for reduction (+ : Qk)
  
  for(y=0;y<N;y++){
    Qk+= v[y]*sin(PI*m*y/N);
  }
  E=0.5*(pow(cte*Qk,2) + w*pow(cte*Qp,2));
  return E;
}




double *F(double  *x)
{
    double *ac= malloc(N*sizeof(double));
    int n;
    ac[0]=0;
    ac[N-1]=0;
    #pragma omp parallel for shared(x,ac)    
    for(n=1;n<N-1;n++)
    {
        ac[n]=(x[n+1]-2*x[n]+x[n-1])+B*(pow(x[n+1]-x[n],2)-pow(x[n]-x[n-1],2));    
    }
    return ac;    
}

