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


//#define HESSIAN

//#define SADDLES

//#define STDOUT

//#define GD

//#define CHECK



int main(int argc, char *argv[]){

    struct parameters p;
    struct variables v;
    Input_Parameters(argc, argv, &v, &p);

	char magn[400];
	sprintf(magn,"%s/magnetization.dat",p.dir); 
	FILE *fmagn = fopen(magn,"w");
	char coupl[400];
	sprintf(coupl,"%s/couplings.dat",p.dir); 
	FILE *fcoupl = fopen(coupl,"w");
	char tap[400];
	sprintf(tap,"%s/gradTAP.dat",p.dir); 
	FILE *ftap = fopen(tap,"w");
	char gradHm[400];
	sprintf(gradHm,"%s/gradHm.dat",p.dir); 
	FILE *fgradHm= fopen(gradHm,"w");

	#ifdef STDOUT
    p.fout=stdout;
    p.fout2=stdout;
	#endif

	fprintf(p.fout,"# a2 = %f a3 = %f D2 = %f D3 = %f seedJ = %d seedS = %d seedX = %d\n", p.a2,p.a3,p.D2,p.D3,p.seedJ,p.seedS,p.seedX); 

	v.S1 = Initialize_System(&v,&p);


	v.E1 = p.TotEN(v.S1,&v,&p);
	v.PG1 = p.TotProjGR(v.S1,&v,&p);

	if(!strcmp(p.method,"GD")) {
		Initialize_GD(&v, &p);
	} 
	else if(!strcmp(p.method,"CG")) {

		Initialize_CG(&v, &p);
	}

	//v.E2=v.E1; SaveSystem(&v,&p); getchar();


	v.T = p.Tin;
	int k=1;
	v.S3=v.S1;

	for(int t=p.ITime;v.Time<(p.ITime+p.TTime);t++) {

		if(!strcmp(p.method,"GD")) {

			if(v.T>0) { 	
						v.PX = Projected_Noise(v.S1,&v,&p,v.T);  
						v.PG1 = Sum(v.PG1,v.PX);
						v.alpha = v.DTime/p.sqrtN*sqrt(Norm2(v.PG1));
			}

			v.S2 = Langevin_step(v.S1,v.PG1,&v,&p);

		} else if(!strcmp(p.method,"CG")) {

			v.S2 = CG_Step(v.S1,v.PG1,&v,&p);
		}

		//if(t==4001) { v.S2 = v.M; v.T = 0; p.Tfin = 0; p.Temp = 0; cout << "normS1 " << Norm2(v.S2) << endl;  }

		v.E2 = p.TotEN(v.S2,&v,&p);
		v.PG2 = p.TotProjGR(v.S2,&v,&p);


		v.Time += v.DTime;
        //v.NormPG = sqrt(Norm2(v.PG1));
		
		if(fabs(v.E2-v.E1)<0.00000001 && fabs(v.E1-v.E0)<0.00000001) { v.Time=p.TTime; } 

		v.Time_prev=v.Time;

		Evaluate_Observables(&v, &p);
		Print_Observables(&v, &p);

		v.E0=v.E1;

		v.S1=v.S2; v.E1=v.E2; v.PG1=v.PG2; v.NPG1=v.NPG2;

		v.CG0=v.CG1; v.H1=v.H2; v.RPG0=v.RPG1; v.RCG0=v.RCG1;

		if(t==0) { v.M = v.S1; }

		v.M = Magnetization(v.S1, v.M, &v, &p); 
		//printf("Time: %f Norm2_M = %f Temp = %f\n",v.Time,Norm2(v.M),p.Temp);

		//if(t%10==0 && t<4001) { v.M = Magnetization(v.S1, v.M, &v, &p); printf("Norm2_M = %f\n", Norm2(v.M));/*Print_Vector_On_File(v.M,fmagn);*/ }
		//else if (t>4000) { v.M = Magnetization(v.S1, v.M, &v, &p); printf("after Norm2_M = %f\n", Norm2(v.M));/*Print_Vector_On_File(v.M,fmagn);*/ }

		//if(t%50 ==0) { char fname[400]; sprintf(fname,"%s/time%d_%f.v",p.dir,t,v.Time); SnapVector(v.S1,fname); }/*printf(fname,"Write file time%d_%f.v",time,v.Time);*/

	}
	Print_Vector_On_File(v.M,fmagn);
	//Print_Couplings_On_File(3,v.J3,fcoupl);

	vector <double> GTAP = Gradient_TAP(&v,&p);
	printf("Norm2_GTAP = %f Norm2_M = %f\n", Norm2(GTAP)/p.N,Norm2(v.M)/p.N);
	Print_Vector_On_File(Total_Gradient(v.M,&v,&p),fgradHm);
	Print_Vector_On_File(GTAP,ftap);

	while(Norm2(GTAP)/p.N>p.alpha0*0.001) {

		GTAP = Gradient_TAP(&v,&p);
		v.M = Lin(1,v.M,0.001,GTAP);
	}

	printf("Norm2_GTAP final = %f Norm2_M = %f\n", Norm2(GTAP)/p.N,Norm2(v.M)/p.N);


	v.G1=GTAP;
	v.H1=Total_Hessian(v.M,&v,&p);

	char fname[400];
	sprintf(p.fnameEi,"%s_eigen_alpha%g.dat",p.fname,p.alpha0);
    p.foutEi=fopen(p.fnameEi,"w");
	/*sprintf(fname,"%s.hessian",p.finitial); SnapVector(v.H1,fname);
	sprintf(fname,"%s.gradient",p.finitial); SnapVector(v.G1,fname);
	sprintf(fname,"%s.config",p.finitial); SnapVector(v.M,fname);
	vector<double> V(1); V[0]=v.E1;
	sprintf(fname,"%s.energy",p.finitial); SnapVector(V,fname);*/

	#ifdef LAPACKEE
		Evaluate_eigenvalues(v.H1, v.M, v.M, p.N, v.Ei, v.gEi);
		Print_Eigenvalues(v.Ei, v.gEi, p.foutEi);
		sprintf(fname,"%s.eigen",p.finitial); SnapVector(v.Ei,fname);
		sprintf(fname,"%s.projeigen",p.finitial); SnapVector(v.gEi,fname);
	#else
		printf("LAPACKEE not activated\n");
	#endif
		

	fclose(ftap);
	fclose(fmagn);
	fclose(fcoupl);

	fprintf(p.fout, "\n");

	SaveSystem(&v,&p); 
	SaveVector(&v,&p,v.S1);



    return 0; 
}
