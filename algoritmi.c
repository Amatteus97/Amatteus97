#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

struct vec{
  double x;
  double y;
};
typedef struct vec vec;



void eulero(double, double, double, double, double, double,double,  int);
void eulero_cromer(double, double, double, double, double,double,double, int);
void verlet(double, double, double, double, double,double,double, int);
vec sumV(vec, vec);
double scalarV(vec, vec);
vec prodV(vec , double);
vec f(vec, double);
void RK2(vec, double, double, double, double,int);
void RK4(vec, double, double,double, double, int, double *);
void periodo(double * , int, double);
double * ray( vec , double * , int);

int main(int argc, char *argv[]){
  double K, m, v0, x0, tmax, Dt, omega2, E0, size;
  int  npassi, alg;
  vec y;
  double * ay;
  // scanf("%lf %lf %lf %lf %lf %lf %d", &K, &m, &v0, &x0, &tmax, &Dt, &alg);

  if(argc!=8){
    fprintf(stderr, "\n%s <K> <m> <v0> <x0> <tmax> <Dt> <alg> <size>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  K=atof(argv[1]);
  m=atof(argv[2]);
  v0=atof(argv[3]);
  x0=atof(argv[4]);
  tmax=atof(argv[5]);
  Dt=atof(argv[6]);
  alg=atoi(argv[7]);
  
  y.x=x0;
  y.y=v0;

  npassi=(tmax/Dt)+1;
      
  omega2=K/m;

  E0=(K*x0*x0+m*v0*v0)/2;
  

  
  ay = (double *)calloc(npassi, sizeof(double));	
  if(alg==1){
  eulero(K, m ,v0, x0, E0, omega2, Dt, npassi);
  }
  if(alg==2){
  eulero_cromer(K, m ,v0, x0, E0, omega2, Dt, npassi);
  }
  if(alg==3){
  verlet(K, m ,v0, x0, E0, omega2, Dt, npassi);
  }
  if(alg==4){
    RK2(y,E0,Dt,K,m,npassi) ;
  }
  if(alg==5){
    RK4(y,E0,Dt,K,m,npassi, ay) ;
  }


}


void eulero(double K, double m, double v0, double x0, double E0,double omega2, double Dt, int npassi){
  double tmp;
  int i;
  double Ef, errE;
   for(i=0; i<npassi; i++){
    tmp=x0 + v0*Dt;
    v0= v0-omega2*x0*Dt;
    x0=tmp;
    Ef=(K*x0*x0+m*v0*v0)/2;
    errE=(Ef-E0)/E0;
    printf("\n%.5g    %.5g     %.5g       %.5g    %.5g",i*Dt, x0, v0, Ef, errE);
  }
}


void eulero_cromer(double K, double m, double v0, double x0, double E0,double omega2, double Dt, int npassi){
  int i;
  double Ef, errE;
  for(i=0;i<npassi;i++){

    v0= v0-omega2*x0*Dt;
    x0=x0 + v0*Dt;
    Ef=(K*x0*x0+m*v0*v0)/2;
    errE=(Ef-E0)/E0;
    printf("\n%.5g    %.5g     %.5g       %.5g    %.5g",i*Dt, x0, v0, Ef, errE);
  }
}

void verlet(double K, double m, double v0, double x0, double E0,double omega2, double Dt, int npassi){
   int i;
   double tmp;
   double Ef, errE;
  for(i=0;i<npassi;i++){
    tmp=x0;
    x0=x0 + v0*Dt-(0.5*omega2*x0*Dt*Dt);
    v0= v0-omega2*(x0+tmp)*Dt*0.5;
    Ef=(K*x0*x0+m*v0*v0)/2;
    errE=(Ef-E0)/E0;
    printf("\n%.5g    %.5g     %.5g       %.5g    %.5g",i*Dt, x0, v0, Ef, errE);
  }
}

vec sumV(vec x, vec y){
  vec z;
  z.x=x.x+y.x;
  z.y=x.y+y.y;
  return z;
}

double scalarV(vec x, vec y){
  double z;
  z=sqrt(x.x*y.x+x.y*y.y);
  return z;
}

vec prodV(vec x, double a){
  vec z;
  z.x=a*x.x;
  z.y=a*x.y;
  return z;
}

vec f(vec y, double o2){
  vec z={y.y, -o2*y.x};
  return z;
}

void RK2(vec y, double E0, double Dt, double K, double m,int npassi){
  int i;
  double Ef, errE;
  for(i=0; i<npassi; i++){
    y=sumV(y, prodV(f(sumV(y, prodV(f(y,K/m),Dt*0.5)), K/m),Dt));
    Ef=(K*y.x*y.x+m*y.y*y.y)*0.5;
    errE=(Ef-E0)/E0;
    printf("\n%.5g    %.5g     %.5g       %.5g    %.5g",i*Dt, y.x, y.y, Ef, errE);
  }
}
  
void RK4(vec y, double E0, double Dt, double K, double m,int npassi, double * ay){
  int i;
  double Ef, errE;
  vec k1, k2, k3, k4;
  for(i=0; i<npassi; i++){
    k1=f(y,K/m);
    k2=f(sumV(y,prodV(k1,Dt*0.5)),K/m);
    k3=f(sumV(y,prodV(k2,Dt*0.5)),K/m);
    k4=f(sumV(y,prodV(k3,Dt)),K/m);
    y=sumV(y,prodV(sumV(sumV(k1,k4),prodV(sumV(k2,k3),2)),Dt/6));
    Ef=(K*y.x*y.x+m*y.y*y.y)*0.5;
    errE=(Ef-E0)/E0;
    ray(y, ay, i); 
    
    //printf("\n%.5g    %.5g     %.5g       %.5g    %.5g",i*Dt, y.x, y.y, Ef, errE);
  }
  periodo(ay, npassi, Dt);	
} 

double * ray( vec y, double * x, int i){
    x[i]=y.x;
    return x;	
}    			

void periodo(double * ay, int npassi, double Dt){
  int i, n=0;
  double t0;
    			
  for(i=0; i<npassi; i++){	
    if(ay[i] * ay[i+1]<=0.0){
      t0=i*Dt-(((i+1)*Dt-i*Dt)*ay[i]/(ay[i+1]-ay[i]));
      n++;	
      printf("\n%d  %.5g",n, t0);
      if(ay[i+1]==0){
      	i++;
      	}	
      }
  }
}        						
