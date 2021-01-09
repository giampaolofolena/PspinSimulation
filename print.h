void Print_Initial_Observables(variables *v, parameters *p) {

	#ifdef CHECK
	cout << "#normS" << Norm2(v.S1) << endl;
	cout << "#normG" << Norm2(v.G1) << endl;
	cout << "#ENERGIA_afterPLANT " << v.E1/p.N << endl;
	cout << "#alpha_STEP " << v.alpha << endl;
	#endif

	fprintf(p->fout, "#v->Time v->q12/N v->E2/N v->Mu1/N sqrt(v->NPG/N) v->e2/N v->alpha v->T");
	#ifdef HESSIAN
	fprintf(p->fout, "v->eigen00 v->eigen11 v->DirLaloux v->DirS v->DirS0");
	#else
    fprintf(p->fout, "\n");
	#endif

	int N=p->N;

    double Ene,Mu;
	if(p->potential=='e') { Ene = v->E1/N; Mu = v->Mu1/N; } 
	else { Ene = Total_Energy(v->S1,v,p)/N; v->G1 = Total_Gradient(v->S1,v,p); Mu = -Prod(v->S1,v->G1)/N; }

	fprintf(p->fout, "%12.11f ",v->Time);
	fprintf(p->fout, "%12.11f ",v->q12/N);
	fprintf(p->fout, "%12.11f ",Ene);
	fprintf(p->fout, "%12.11f ",Mu);
	fprintf(p->fout, "%12.11f ",sqrt(v->NPG1/N));
	fprintf(p->fout, "%12.11f ",v->alpha);
	fprintf(p->fout, "%12.11f ",v->T);
	fprintf(p->fout, "%12.11f ",v->q13/N);

	#ifdef HESSIAN
	fprintf(p->fout, "%12.11f %12.11f %12.11f %12.11f %12.11f\n",v->eigen00,v->eigen11,v->DirLaloux,v->DirS,v->DirS0);
	#else
    fprintf(p->fout, "\n");
	#endif
}

void Print_Observables(variables *v, parameters *p) {

	int N=p->N;

    double Ene,Mu;
	if(p->potential=='e') { Ene = v->E2/N; Mu = v->Mu1/N; } 
	else { Ene = Total_Energy(v->S2,v,p)/N; v->G1 = Total_Gradient(v->S1,v,p); Mu = -Prod(v->S1,v->G1)/N; }

	fprintf(p->fout, "%12.11f ",v->Time);
	fprintf(p->fout, "%12.11f ",v->q12/N);
	fprintf(p->fout, "%12.11f ",Ene);
	fprintf(p->fout, "%12.11f ",Mu);
	fprintf(p->fout, "%12.11f ",sqrt(v->NPG1/N));
	fprintf(p->fout, "%12.11f ",v->e2/N);
	fprintf(p->fout, "%12.11f ",v->alpha);
	fprintf(p->fout, "%12.11f ",v->T);
	fprintf(p->fout, "%12.11f ",v->q13/N);

	if(p->potential!='e') {
		fprintf(p->fout, "%12.11f ",v->E2/N);
		fprintf(p->fout, "%12.11f ",v->Mu1/N);
	}

	#ifdef HESSIAN
	fprintf(p->fout, "%12.11f %12.11f %12.11f %12.11f %12.11f\n",v->eigen00,v->eigen11,v->DirLaloux,v->DirS,v->DirS0);
	#else
    fprintf(p->fout, "\n");
	#endif
}

void Print_Eigen(variables *v, parameters *p) {
	int N=p->N;
	fprintf(p->foutEi, "#%12.11f %12.11f %12.11f\n",v->Time,v->E1/N,v->Mu1/N);
    for(int i = 0; i < N; i++) {
      	fprintf(p->foutEi, "%12.11f %12.11f\n",v->Ei[i],fabs(v->gEi[i]));
    }
    fprintf(p->foutEi, "\n");
}

void Print_Eigenvalues(vector<double> Ei, vector<double> gEi, FILE *fname) {
	int N=Ei.size();
    for(int i = 0; i < N; i++) {
      	fprintf(fname, "%12.11f %12.11f\n",Ei[i],fabs(gEi[i]));
    }
    fprintf(fname, "\n");
}

void Print_Vector(vector<double> V) {
	printf("SIZE %d\n",V.size());
    for(int i = 0; i < 30; i++) {
      	printf("%12.11f ",V[i]);
    }
    printf("\n");
}

void Print_Vector_On_File(vector<double> V, FILE *fname) {
    for(int i = 0; i < V.size(); i++) {
      	fprintf(fname,"%12.11f \n",V[i]);
    }
    fprintf(fname, "\n");
}

void Print_Couplings(vector<struct node> &J) {
	printf("SIZE %d\n",J.size());
    for(int i = 0; i < 30; i++) {
      	printf("%12.11f ",J[i].J);
    }
    printf("\n");
}

void Print_Couplings_On_File(int p, vector<struct node> &JJ, FILE *fname) {

	for(size_t i = 0; i < JJ.size(); i++) {
        double J = JJ[i].J;
        long long int SSS = JJ[i].SSS;

        long long int mask = 0b1111111111111111;
        double multiS = 1.;

		fprintf(fname,"%12.11f ",J);
        for(int j=0;j<p;j++) { 
        	fprintf(fname,"%lld ",(SSS>>(16*j) & mask)); 
        }
        fprintf(fname,"\n");
    }
}