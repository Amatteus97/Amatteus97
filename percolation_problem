#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int * net;


struct root{
  int data;
  int h;
};
 
typedef struct root root;


int buco(double, int);
void bordo(int);
void capo(int , int, int * );
int perc(int, root *);
void taglia(int *, int);
int label(int );
int maxsize(int *, int);  

  
int main(int argc, char *argv[]){


  unsigned int seed;
  int L, N, i ,nconf, q, j, v, * size, sumc=0;
  double Pmin, Pmax, DP;
  root *vec;
  
  
  if(argc!=6){
    fprintf(stderr, "\n%s <L> <Pmin> <Pmax> <DP> <nconf>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  L=atoi(argv[1]);
  Pmin=atof(argv[2]);
  Pmax=atof(argv[3]);
  DP=atof(argv[4]);
  nconf=atoi(argv[5]);
  

  FILE *devram=fopen("/dev/random", "r");
  fread(&seed, sizeof(unsigned int), 1, devram);
  fclose(devram);
  srand48(seed);

 

  N=L*L;
  net=(int *)calloc(N, sizeof(int));
  q=(int)(Pmax-Pmin)/DP;
  size=(int *)calloc(N, sizeof(int));
  vec=(root *)calloc(N, 2*sizeof(int));
  
  for(j=0; j<=q; j++){
    for(v=0; v<nconf; v++){
      for(i=0; i<N; i++){
	net[i]=buco(j*DP+Pmin,i);
      }
      bordo(L);
      for(i=0; i<N; i++){
	size[i]=0;
      }
      taglia(size, L);
      sumc+=maxsize(size, N);
     
      /*for(i=0; i<N; i++){
	if((i+1)<N && net[i+1]>=0 && net[i+1]!=net[i]){
	  capo(i, i+1, size);
	}
	if((i+L)<N && net[i+L]>=0 && net[i+L]!=net[i] && net[label(i+1)] != net[label(i+L)]){
	  capo(i, i+L, size);
	}
	//printf("%d\n", size[i]);
	}*/
     
      for(i=0; i<N; i++){
	vec[i].data=net[i];
	vec[i].h=(i/L);
      }
      for(i=0; i<L; i++){
	vec[label(i)].h=0;
      }
      if(perc(L, vec) == 1){
	//printf("La configurazione %d con probabilità %lf ha percolato\n", v+1, j*DP+Pmin);
      }
      
    }printf("%lf %lf\n",j*DP+Pmin, (double)sumc/nconf);
      sumc=0;
    
    
  }

}

int buco(double p, int i){
  int x;
  if(((double)lrand48()/RAND_MAX)<p){
    x=i;
  }else{
    x=-1;
  }
  return x;
}



void bordo(int L){
  int nchanges, i, N;
 
  N=L*L;

  do{
    nchanges=0;
    for(i=0; i<N; i++){
      if(net[i]>=0){
	if((i+1)<N && net[i+1]>=0 && net[i+1]!=net[i]){
	  nchanges++;
	  if(net[i+1]<net[i]){
	    net[i]=net[i+1];
	  }else{
	    net[i+1]=net[i];
	  }
	}
	if((i+L)<N && net[i+L]>=0 && net[i+L]!=net[i]){
	  nchanges++;
	  if(net[i+L]<net[i]){
	    net[i]=net[i+L];
	  }else{
	    net[i+L]=net[i];
	  }
	}
      }
    }
  }while(nchanges!=0);
}
	  



int label(int i){
  while(net[i]!=i && net[i]>=0)
    i=net[i];
  return i;
}

void capo(int i, int k, int* size){
  if(label(i)!=label(k)){
    if(size[label(i)]<size[label(k)]){
      net[label(i)]=label(k);
      size[label(k)]+=size[label(i)];
    }else{
      net[label(k)]=label(i);
      size[label(i)]+=size[label(k)];
    } 
  } 
}

int perc(int L, root *min){
  int isperc=0, N=L*L, i;
  for(i=L*(L-1); i<N; i++){
    if(min[label(i)].h == 0){
      isperc=1;
    }
  }
  return isperc;
}


void taglia(int *clsize, int L){
  int i, N;
  N=L*L;
  for(i=0; i<N; i++){
    if(net[i]>=0){
      clsize[net[i]]++;
    }
  }
}


int maxsize(int *size, int N){
  int i, max=0;
  for(i=0; i<N; i++){
    if(size[i] >= max && size[i] >= 0) max=size[i];
  }
  return max;
}



    
    
