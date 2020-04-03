void Initialize_GD(variables *v, parameters *p) {

	v->alpha = p->alpha0;
	v->Time = 0;
	v->DTime = p->DT0;
	v->PGPG = 1;

	v->S0=v->S1;
	v->S2=v->S1;
	v->NPG1 = Norm2(v->PG1);

	Evaluate_Observables(v, p);
	Print_Initial_Observables(v, p); 
}

vector<double> Langevin_step(vector<double> S1, vector<double> PG1, variables *v, parameters *p) {

	vector<double> S2;
	
	if(fabs(Prod(S1,PG1))>=0.00000001) {
        cout << "LangevinStep: NON orthogonal > SOFT STEP" << endl;
        v->alpha*=0.9;
        if(p->Temp>0) {
        	S2 = SoftLangevinSTEP(S1,Sum(v->X,v->G1),v->alpha);
        } else {
        	S2 = SoftLangevinSTEP(S1,v->G1,v->alpha);
        }
    }
    else {
		S2 = HardLangevinSTEP(S1,PG1,v->alpha);
	}

	if(p->Temp==0) { 
		if(v->PGPG<p->PGPGpar) { v->alpha*=0.97; } else if (v->PGPG>p->PGPGpar) { v->alpha/=0.97; } 
		v->DTime = sqrt(Norm2(Sub(S2,S1))/Norm2(PG1));
	}

	return S2;
}

void Initialize_CG(variables *v, parameters *p) {

	v->alphaCG = p->alpha0_CG;
	v->Time = 0;
	v->DTime = p->DT0;
	v->PGPG = 1;

	v->S0=v->S1;
	v->S2=v->S1;
	v->NPG1 = Norm2(v->PG1);

	Evaluate_Observables(v, p);
	Print_Initial_Observables(v, p); 

	v->RPG0=v->PG1;
	v->CG1=v->PG1;
	v->CG0=v->PG1;
	v->RCG0=v->PG1;
}

vector<double> CG_Step(vector<double> S1, vector<double> PG1, variables *v, parameters *p) {
	
	v->CG1 = NewConjugateGradient(PG1,v->RPG0,v->RCG0);

	if(fabs(Prod(S1,v->CG1))>=0.00000001) {
   		cout << "NON orthogonal" << endl;
        	//v->alpha*=0.9;
	}

	v->alphaCG = Select_CG_angle(S1,v->CG1,v->alphaCG*1.25,v,p);

	v->RPG1 = RotateGeneralVector(PG1,S1,v->CG1,v->alphaCG);
	v->RCG1 = RotateGeneralVector(v->CG1,S1,v->CG1,v->alphaCG);

	v->alpha=v->alphaCG;

	return RotateVector(S1,v->CG1,v->alphaCG);
}

double Select_CG_angle(vector<double> S, vector<double> CG, double alpha_max, variables *v, parameters *p){

    vector<double> S2;
    double normCG = sqrt(Norm2(CG));

    double invphi2 = 0.6180339887498949; //(sqrt(5) - 1) / 2;                                                                                                              
    double invphi = 1-invphi2;         //           0.3819660112501051

    //Given a function f with a single local minimum in                                                                                                       
    //the interval [a,b], gss returns a subset interval                                                                                                       
    //[c,d] that contains the minimum with d-c <= tol.                                                                                               

    double tol = normCG/1000.; //1e-4;
    //printf("tol%e %f ",tol,normCG);

    double a,b,c,d,h,ya,yb,yc,yd;
    int n,k;

    a=0;
    b=alpha_max/invphi2;

    h = b - a;

    c = a + invphi * h;
    d = a + invphi2 * h;

    S2 = RotateVector(S,CG,c);
    yc = p->TotEN(S,v,p);
    
    S2 = RotateVector(S,CG,d);
    yd = p->TotEN(S,v,p);

    //printf("3s %f \n",Norm2(Ss->Rx,Ss->N));

    k=0;

    while(h>tol) {
        if (yc > yd) {
            a = c;
            ya = yc;
            c = d;
            yc = yd;
            h = b - a;
            d = b - invphi * h;

            S2 = RotateVector(S,CG,d);
            yd = p->TotEN(S,v,p);
        }
        else {
            b = d;
            yb = yd;
            d = c;
            yd = yc;
            h = b - a;
            c = a + invphi * h;

            S2 = RotateVector(S,CG,c);
            yc = p->TotEN(S,v,p);
        }

        k++;
    }

    if (yc > yd) { return d; } else { return c; }
}



/*void DAFARE(variables *p, parameters *v) {

		#ifdef HESSIAN
		H2 = Lin(a2,Hessian(2,J2,S2),a3,Hessian(3,J3,S2));		
		P = Id(N);
    	eigen0 = MinEigen(H2,P,mE0,HighestEigen);
		P = Proj(S2);	
    	eigen1 = MinEigen(H2,P,mE1,HighestEigen);


		eigen00 = Sandwich(H2,mE0);
		eigen11 = Sandwich(H2,mE1);
		DirLaloux = Sandwich(H2,PG1);
		DirS = Sandwich(H2S2);
		DirS0 = Sandwich(H2,S0);
		#endif

				#ifdef LAPACKEE
    	H2 = Lin(a2,Hessian(2,J2,S2),a3,Hessian(3,J3,S2));
    	DH = Sub(H2,H1);
    	Evaluate_eigenvalues(H2,DH,S1,PG1,N,mE0,mE1);
		Print_Eigen(&PAR, &VAR);
		#endif
}*/
