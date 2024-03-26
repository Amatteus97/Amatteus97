#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct vec{
  int x, y;
};

typedef struct vec vec;

void siti(int *, int, double, vec *);
void mov(int *, int, vec *, int);
int bordo(int , int);

int main(int argc, char *argv[]){
  
  int L, i, t, M, seed, N, *occ, j;
  double p;
  vec *avec;
  
  if(argc!=4){
    fprintf(stderr, "\n%s <L> <p> <tmax>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  L=atoi(argv[1]);
  p=atof(argv[2]);
  t=atoi(argv[3]);
  
  
  
  FILE *devram=fopen("/dev/random", "r");
  fread(&seed, sizeof(unsigned int), 1, devram);
  fclose(devram);
  srand48(seed);

  N=L*L;
  occ=(int*)calloc(N, sizeof(int));

  for(i=0; i<N; i++){
    occ[i]=-1;
  }
  M=(int)(p*N);
    
  avec=(vec *)calloc(M, sizeof(vec));

  siti(occ, L, p, avec);
  
  /* for(i=0; i<N; i++){
    printf("%d %d ", occ[i], M);
    printf("%d %d \n", avec[i%M].x, avec[i%M].y);
    }*/




  for(i=0; i<t; i++){
    for(j=0; j<N; j++){
      occ[j]=0;
    }
    mov(occ, L, avec, M);
    for(j=0; j<M; j++){
      //printf("%d\n", occ[i]);
      //printf("%d %d\n", avec[j].x, avec[j].y);
    }
    //printf("\n\n");
  }

  for(i=0; i<M; i++){
    //printf("%d\n", occ[i]);
    //printf("%d %d\n", avec[i].x, avec[i].y);
  }
  
}


void siti(int *occ, int L, double p, vec *v){
  int i=0, M, k=0, N;
  N=L*L;
  M=(int)(p*N);

  do{
    if(drand48()<p && occ[i]<0){
      occ[i]=k;
      v[k].x=(i%L);
      v[k].y=(int)(i/L);
      //printf("%d %d %d\n", occ[i], v[k].x, v[k].y);
      k++;
    }
    i=(i+1)%N;
  }while(k<M);
}


void mov(int *occ, int L, vec *v, int M){
  int i, ind, r, m, n;

  for(i=0; i<M; i++) {
    m=v[i].x;
    n=v[i].y;
    r=drand48();
    if(r<0.25){
      v[i].y++;
    }else{
      if(r<0.50){
	v[i].x++;
      }else{
	if(r<0.75){
	  v[i].y--;
	}else{
	  v[i].x--;
	}
      }
    }
    ind=bordo(L*v[i].y+v[i].x, L);
    v[i].x=(ind%L);
    v[i].y=((int)(ind/L));
    //printf("%d\n", ind);
    //printf("%d %d %d %d \n",(ind%L) , ((int)(ind/L)), v[i].x, v[i].y);
    occ[ind]++;
    //printf("%d\n", occ[ind]);
    if(occ[ind] > 1){
      if(occ[bordo(L*(v[i].y+1)+(v[i].x), L)] <1 || occ[bordo(L*(v[i].y-1)+(v[i].x), L)] <1 || occ[bordo(L*(v[i].y)+(v[i].x+1), L)] <1|| occ[bordo(L*(v[i].y)+(v[i].x-1), L)] <1){
	do{
	  v[i].x=m;
	  v[i].y=n;
	  //printf("%d %d\n", v[i].x, v[i].y);
	  r=drand48();
	  if(r<0.25){
	    v[i].y++;
	  }else{
	    if(r<0.50){
	      v[i].x++;
	    }else{
	      if(r<0.75){
		v[i].y--;
	      }else{
		v[i].x--;
	      }
	    }
	  }
	  //printf("%d %d\n", v[i].x, v[i].y);
	  ind=bordo(L*v[i].y+v[i].x, L);
	  occ[ind]++;
	  //printf("%d %d\n", ind, occ[ind]);
	}while(occ[ind]>1);
      }else{
	v[i].x=m;
	v[i].y=n;
      }
    }
  }
}    

int bordo(int ind, int L){
  if((ind%L) >= L){
    ind=L*((int)(ind/L));
  }else{
    if((ind%L) < 0){
      ind=L*((int)(ind/L))+(L-1);
    }else{
      if(((int)(ind/L)) >= L){
	ind=(ind%L);
      }else{
	if((int)(ind/L) < 0){
	  ind=(L*(L-1)+(ind%L));
	}
      }
    }
  }
  return ind;
}
  
