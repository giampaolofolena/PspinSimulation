#include <random>
#include <iostream>
#include <time.h>
#include <algorithm>
#include <vector>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <sys/stat.h>

#define LAPACKEE
//gpf$ g++ -o HeToEi /usr/local/opt/lapack/lib/liblapacke.dylib HessToEigen.cpp -I/usr/local/opt/lapack/include -I/usr/local/opt/openblas/include -llapack -lcblas

using namespace std;

#include "declarations.h"

#include "functions.h"

#include "parameters.h"
#include "variables.h"

#include "print.h"
#include "inout.h"

#include "observables.h"

#include "initialization.h"
#include "dynamics.h"


//#define HESSIAN

//#define SADDLES

//#define STDOUT

//#define GD

//#define CHECK



int main(int argc, char *argv[]){

	int N=-1;
	char finitial[200];

	vector<double> H;
	vector<double> G;
	vector<double> S;
	vector<double> E;

	vector<double> Ei;	
	vector<double> gEi;

	char *cvalue = NULL;
    int index;
    int c;
    opterr = 0;

	while ((c = getopt (argc, argv, "N:i:")) != -1) //COLON MEANS THAT A VALUE MUST BE SPECIFIED
    switch (c)
      {
      case 'N': //NUMBER OF SPINS
        cvalue = optarg;
        N=atoi(cvalue); 	if(N>pow(2,16)) { printf("ACTUNG: maximal size of the system exceeded! look at initialization of J"); getchar(); }      
        break;
      case 'i': //initial configuration and parameters from file
        cvalue = optarg;
        sprintf(finitial,"%s",cvalue);
        break;
      case '?':
        fprintf (stderr,"Unknown option character %c.\n",optopt);
        getchar();
        break;
        default:
        abort ();
      }

    if(N==-1 || *finitial=='0') { printf("must specify option N and i!\n"); }

    printf("N = %d; INPUT FILE: %s\n",N,finitial);

	char fname[300];
	sprintf(fname,"%s.hessian",finitial);  H=OpenVector(N*N,fname);
	sprintf(fname,"%s.gradient",finitial); G=OpenVector(N,fname);
	sprintf(fname,"%s.config",finitial);   S=OpenVector(N,fname);
	sprintf(fname,"%s.energy",finitial);   E=OpenVector(1,fname);
	Evaluate_eigenvalues(H, G, S, N, Ei, gEi);
	sprintf(fname,"%s.eigen.print",finitial); 	FILE *fn = fopen(fname,"w"); Print_Eigenvalues(Ei, gEi, fn);
	sprintf(fname,"%s.eigen",finitial); 		SnapVector(Ei,fname);
	sprintf(fname,"%s.projeigen",finitial); 	SnapVector(gEi,fname);
		
}
