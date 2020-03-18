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

	v->alphaCG = Select_CG_angle(S1,v->CG1,p->a2,v->J2,p->a3,v->J3,v->alphaCG*1.25);

	v->RPG1 = RotateGeneralVector(PG1,S1,v->CG1,v->alphaCG);
	v->RCG1 = RotateGeneralVector(v->CG1,S1,v->CG1,v->alphaCG);

	v->alpha=v->alphaCG;

	return RotateVector(S1,v->CG1,v->alphaCG);
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
