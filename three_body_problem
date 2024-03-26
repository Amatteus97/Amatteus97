#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// ./pianeti.exe  3.0028*10^-6 1 0 0 4 3.226*10^-7 1.5 0 0 4 50 0.0001 5 >data.dat
// Valori in UA e anni

struct vec{
  double x;
  double y;
};

typedef struct vec vec;

struct gen{
  vec p;
  vec a;

};

typedef struct gen gen;



vec sumV(vec, vec);
double scalarV(vec, vec);
vec prodV(vec , double);
gen sum (gen, gen);
gen prod (gen, double);
double dist (vec, vec);
gen f(gen, gen, double, double, double);
void eulero(gen, gen, double, double, double, int, double);
void eulero_cromer(gen, gen, double, double, double, int, double);
void verlet(gen, gen, double, double, double, int, double);
void RK2(gen, gen, double, double, double, int, double);
void RK4(gen, gen, double, double, double, int, double);

int main(int argc, char *argv[]){
  double accen,muA, xa, ya, vxa, vya, muB, xb, yb, vxb, vyb, tmax, Dt, gamma;
  int  npassi, alg;
  vec posA, velA, posB, velB;
  gen corpoA, corpoB;

  if(argc!=14){
    fprintf(stderr, "\n%s <ua> <xa> <ya> <vxa> <vya> <ub> <xb> <yb> <vxb> <vyb> <tmax> <Dt> <alg>", argv[0]);
    exit(EXIT_FAILURE);
  }
  muA=atof(argv[1]);
  xa=atof(argv[2]);
  ya=atof(argv[3]);
  vxa=atof(argv[4]);
  vya=atof(argv[5]);
  muB=atof(argv[6]);
  xb=atof(argv[7]);
  yb=atof(argv[8]);
  vxb=atof(argv[9]);
  vyb=atof(argv[10]);
  tmax=atof(argv[11]);
  Dt=atof(argv[12]);
  alg=atoi(argv[13]);
  
  gamma=4*M_PI*M_PI;

  posA.x=xa;
  posA.y=ya;
  velA.x=vxa;
  velA.y=vya;

  posB.x=xb;
  posB.y=yb;
  velB.x=vxb;
  velB.y=vyb;

  corpoA.p=posA;
  corpoA.a=velA;
  
  corpoB.p=posB;
  corpoB.a=velB;

  accen=0;
  npassi=(tmax/Dt)+1;

  if(alg==1){
    eulero(corpoA,corpoB,gamma,muA,muB,npassi,Dt);
 }
  if(alg==2){
    eulero_cromer(corpoA,corpoB,gamma,muA,muB,npassi,Dt);
  }
  if(alg==3){
    verlet(corpoA,corpoB,gamma,muA,muB,npassi,Dt);
  }
  if(alg==4){
    RK2(corpoA,corpoB,gamma,muA,muB,npassi,Dt) ;
  }
  if(alg==5){
    RK4(corpoA,corpoB,gamma,muA,muB,npassi,Dt) ;
  }

}


vec sumV(vec x, vec y){
  vec z;
  z.x=x.x+y.x;
  z.y=x.y+y.y;
  return z;
}

gen sum(gen p, gen a){
  gen res;
  res.p=sumV(p.p,a.p);
  res.a=sumV(p.a,a.a);
  
  return res;

}

double scalarV(vec x, vec y){
  double z;
  z=sqrt(fabs(x.x*y.x+x.y*y.y));
  return z;
}

double dist(vec x, vec y){
  double z;
  z=sqrt((x.x-y.x)*(x.x-y.x)+(x.y-y.y)*(x.y-y.y));
  return z;
}


vec prodV(vec x, double a){
  vec z;
  z.x=a*x.x;
  z.y=a*x.y;
  return z;
}

gen prod(gen p, double a ) {
  gen res;
  res.p=prodV(p.p,a);
  res.a=prodV(p.a,a);
  return res;

}

gen f(gen corpo, gen corpo1, double gamma, double mu, double accen){
  gen res;
  res.p=corpo.a;
  res.a=sumV(prodV(corpo.p,-gamma/((scalarV(corpo.p,corpo.p))*(scalarV(corpo.p,corpo.p))*(scalarV(corpo.p,corpo.p)))),prodV(sumV(corpo1.p,prodV(corpo.p,-1.0)),(gamma*accen*mu)/((dist(corpo1.p,corpo.p))*(dist(corpo1.p,corpo.p))*(dist(corpo1.p,corpo.p)))));
  return res;
}

void eulero (gen corpoA, gen corpoB, double gamma, double muA, double muB, int npassi, double Dt){
  double accen=0;
  int i;
  gen tmp;
  for(i=0; i<npassi; i++){
    tmp=corpoA;
    //accen=((double)i/(double)npassi);
    corpoA.a=sumV(corpoA.a,prodV(f(corpoA,corpoB,gamma,muB,accen).a,Dt));
    corpoA.p=sumV(corpoA.p,prodV(f(corpoA,corpoB,gamma,muB,accen).p,Dt));
    corpoB.a=sumV(corpoB.a,prodV(f(corpoB,tmp,gamma,-muA,accen).a,Dt));
    corpoB.p=sumV(corpoB.p,prodV(f(corpoB,tmp,gamma,-muA,accen).p,Dt));
    
    printf("%lf       %lf         %lf      %lf      %lf\n",i*Dt, (corpoA.p).x, (corpoA.p).y, (corpoB.p).x, (corpoB.p).y);
  }
}



void eulero_cromer(gen corpoA, gen corpoB, double gamma, double muA, double muB, int npassi, double Dt){
  int i;
  gen tmp;
  double accen=0;
  for(i=0;i<npassi;i++){
     //accen=((double)i/(double)npassi);
    tmp=corpoA;
    corpoA.a=sumV(corpoA.a,prodV(f(corpoA,corpoB,gamma,muB,accen).a,Dt));
    corpoA.p=sumV(corpoA.p,prodV(corpoA.a,Dt));
    corpoB.a=sumV(corpoB.a,prodV(f(corpoB,tmp,gamma,-muA,accen).a,Dt));
    corpoB.p=sumV(corpoB.p,prodV(corpoB.a,Dt));
   
    printf("%lf     %lf     %lf     %lf     %lf\n",i*Dt, (corpoA.p).x, (corpoA.p).y, (corpoB.p).x, (corpoB.p).y);
  }
}



void verlet(gen corpoA, gen corpoB, double gamma, double muA, double muB, int npassi, double Dt){
   int i;
   gen tmp, tmp1;
   double accen=0;
   for(i=0;i<npassi;i++){
     // accen=((double)i/(double)npassi);
     tmp=corpoA;
     corpoA.p=sumV(sumV(corpoA.p,prodV(corpoA.a,Dt)),prodV(f(corpoA,corpoB,gamma,muB,accen).a,0.5*Dt*Dt));
     corpoA.a=sumV(corpoA.a,prodV(sumV(f(corpoA,corpoB,gamma,muB,accen).a,f(tmp,corpoB,gamma,muB,accen).a),Dt*0.5));
     tmp1=corpoB;
     corpoB.p=sumV(sumV(corpoB.p,prodV(corpoB.a,Dt)),prodV(f(corpoB,tmp,gamma,-muA,accen).a,0.5*Dt*Dt));
     corpoB.a=sumV(corpoB.a,prodV(sumV(f(corpoB,tmp,gamma,-muA,accen).a,f(tmp1,tmp,gamma,-muA,accen).a),Dt*0.5));
	  
     printf("%lf     %lf     %lf     %lf     %lf\n",i*Dt, (corpoA.p).x, (corpoA.p).y, (corpoB.p).x, (corpoB.p).y);
   }
}


void RK2(gen corpoA, gen corpoB, double gamma, double muA, double muB, int npassi, double Dt){
  int i;
  gen tmp;
  double accen=0;
   for(i=0; i<npassi; i++){
      //accen=((double)i/(double)npassi);
     tmp=corpoA;
     corpoA=sum(corpoA, prod(f(sum(corpoA, prod(f(corpoA,corpoB,gamma,muB,accen),Dt*0.5)),corpoB,gamma,muB,accen),Dt));
     corpoB=sum(corpoB, prod(f(sum(corpoB, prod(f(corpoB,tmp,gamma,-muA,accen),Dt*0.5)),tmp,gamma,-muA,accen),Dt));
     printf("%lf     %lf     %lf     %lf     %lf\n",i*Dt, (corpoA.p).x, (corpoA.p).y, (corpoB.p).x, (corpoB.p).y);
   }
}



void RK4(gen corpoA, gen corpoB, double gamma, double muA, double muB, int npassi, double Dt) {
  int j,i;
  double accen=0;
  gen k1, k2, k3, k4,k1b, k2b, k3b, k4b;

  for(i=0; i<npassi; i++){
      k1=f(corpoA, corpoB, gamma, -muB, accen);
      k1b=f(corpoB, corpoA, gamma, muA, accen);
      k2=f(sum(corpoA,prod(k1,Dt*0.5)),corpoB, gamma, -muB, accen);
      k2b=f(sum(corpoB,prod(k1b,Dt*0.5)),corpoA, gamma,muA, accen);
      k3=f(sum(corpoA,prod(k2,Dt*0.5)),corpoB, gamma, -muB, accen);
      k3b=f(sum(corpoB,prod(k2b,Dt*0.5)),corpoA, gamma, muA, accen);
      k4=f(sum(corpoA,prod(k3,Dt)),corpoB, gamma, -muB, accen);
      k4b=f(sum(corpoB,prod(k3b,Dt)),corpoA, gamma, muA, accen);
      corpoA=sum(corpoA,prod(sum(sum(k1,k4),prod(sum(k2,k3),2)),Dt/6));
      corpoB=sum(corpoB,prod(sum(sum(k1b,k4b),prod(sum(k2b,k3b),2)),Dt/6));
      printf("%lf     %lf     %lf     %lf     %lf\n",accen, (corpoA.p).x, (corpoA.p).y, (corpoB.p).x, (corpoB.p).y);
      accen+=(double)1/(double)(100*npassi);
    }
  do{
    k1=f(corpoA, corpoB, gamma, -muB, accen);
    k1b=f(corpoB, corpoA, gamma, muA, accen);
    k2=f(sum(corpoA,prod(k1,Dt*0.5)),corpoB, gamma, -muB, accen);
    k2b=f(sum(corpoB,prod(k1b,Dt*0.5)),corpoA, gamma,muA, accen);
    k3=f(sum(corpoA,prod(k2,Dt*0.5)),corpoB, gamma, -muB, accen);
    k3b=f(sum(corpoB,prod(k2b,Dt*0.5)),corpoA, gamma, muA, accen);
    k4=f(sum(corpoA,prod(k3,Dt)),corpoB, gamma, -muB, accen);
    k4b=f(sum(corpoB,prod(k3b,Dt)),corpoA, gamma, muA, accen);
    corpoA=sum(corpoA,prod(sum(sum(k1,k4),prod(sum(k2,k3),2)),Dt/6));
    corpoB=sum(corpoB,prod(sum(sum(k1b,k4b),prod(sum(k2b,k3b),2)),Dt/6));
    printf("%lf     %lf     %lf     %lf     %lf\n",accen, (corpoA.p).x, (corpoA.p).y, (corpoB.p).x, (corpoB.p).y);
    accen+=(double)1/(double)(100);
  }while(accen<1);
}
