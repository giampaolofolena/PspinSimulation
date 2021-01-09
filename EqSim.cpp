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

//#define LAPACKEE
//gpf$ g++ -o Ei PsSim.cpp -llapacke -lcblas

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

#include "methods.h"


//#define HESSIAN

//#define SADDLES

//#define STDOUT

//#define GD

//#define CHECK



int main(int argc, char *argv[]){

    struct parameters p;
    struct variables v;
    Input_Parameters(argc, argv, &v, &p);

    //p.Tin=p.Temp;
    //p.Beta=1./p.Temp;

	#ifdef STDOUT
    p.fout=stdout;
	#endif

	fprintf(p.fout,"# a2 = %f a3 = %f D2 = %f D3 = %f seedJ = %d seedS = %d seedX = %d\n", p.a2,p.a3,p.D2,p.D3,p.seedJ,p.seedS,p.seedX); 

	p.Beta = 1./p.Temp;
	v.S1 = Initialize_System(&v,&p);

	v.E1 = p.TotEN(v.S1,&v,&p);
	v.PG1 = p.TotProjGR(v.S1,&v,&p);

	Initialize_GD(&v, &p);
	v.S3=v.S1;

	v.T = p.Temp;

	//vector <double> S_IS0 = Find_Inherent_Structure(v.S1,&v,&p);


	double DT=0;
	v.Time=0;
	while(v.Time<p.TTime) {

		//REPLICA 1

		v.PX = Projected_Noise(v.S1,&v,&p,v.T);  
		v.PG1 = Sum(v.PG1,v.PX);

		v.alpha = v.DTime/p.sqrtN*sqrt(Norm2(v.PG1));

		v.S2 = Langevin_step(v.S1,v.PG1,&v,&p);

		v.E2 = p.TotEN(v.S2,&v,&p);
		v.PG2 = p.TotProjGR(v.S2,&v,&p);

		v.Time += v.DTime;
		v.Time_prev=v.Time;


		Evaluate_Observables(&v, &p);
		Print_Observables(&v, &p);

		v.E0=v.E1;
		v.S1=v.S2; v.E1=v.E2; v.PG1=v.PG2; v.NPG1=v.NPG2;
		v.CG0=v.CG1; v.H1=v.H2; v.RPG0=v.RPG1; v.RCG0=v.RCG1;

		//if(v.Time<p.TTime/2. && v.Time>p.TTime*(1/2.-0.01)) { v.S3=v.S1; }

		//REPLICA 3
		v.PX = Projected_Noise(v.S3,&v,&p,v.T);  
		v.PG2=Sum(p.TotProjGR(v.S3,&v,&p),v.PX);
		v.S3 = RotateVector(v.S3,v.PG2,v.DTime/p.sqrtN*sqrt(Norm2(v.PG2)));

	}

	/*for(int i=0; i<v.IS.size(); i++){
		for(int j=0; j<i; j++) {
			printf("%1.1f ",Prod(v.IS[i],v.IS[j])/p.N); 
		}
		printf("\n");
	}*/

	fprintf(p.fout, "\n");

	SaveSystem(&v,&p); 
	SaveVector(&v,&p,v.S1);	

    return 0; 
}
