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
	fprintf(p->fout, "%12.11f ",v->Time);
	fprintf(p->fout, "%12.11f ",v->q12/N);
	fprintf(p->fout, "%12.11f ",v->E1/N);
	fprintf(p->fout, "%12.11f ",v->Mu1/N);
	fprintf(p->fout, "%12.11f ",sqrt(v->NPG1/N));
	fprintf(p->fout, "%12.11f ",v->alpha);
	fprintf(p->fout, "%12.11f ",v->T);

	#ifdef HESSIAN
	fprintf(p->fout, "%12.11f %12.11f %12.11f %12.11f %12.11f\n",v->eigen00,v->eigen11,v->DirLaloux,v->DirS,v->DirS0);
	#else
    fprintf(p->fout, "\n");
	#endif
}

void Print_Observables(variables *v, parameters *p) {
	int N=p->N;
	fprintf(p->fout, "%12.11f %12.11f %12.11f %12.11f %12.11f %12.11f %12.11f %12.11f ",v->Time,v->q12/N,v->E2/N,v->Mu1/N,sqrt(v->NPG2/N),v->e2/N,v->alpha,v->T);
	#ifdef HESSIAN
	fprintf(p->fout, "%12.11f %12.11f %12.11f %12.11f %12.11f\n",v->eigen00,v->eigen11,v->DirLaloux,v->DirS,v->DirS0);
	#else
    fprintf(p->fout, "\n");
	#endif
}

void Print_Eigen(variables *v, parameters *p) {
	int N=p->N;
	fprintf(p->fout2, "%12.11f %12.11f %12.11f ",v->Time,v->E1/N,v->Mu1/N);
    for(int i = 0; i < N; i++) {
      	fprintf(p->fout2, "%12.11f %12.11f ",v->Ei[i],fabs(v->gEi[i]));
    }
    fprintf(p->fout2, "\n");
}

void Print_Vector(vector<double> V) {
	printf("SIZE %d\n",V.size());
    for(int i = 0; i < 30; i++) {
      	printf("%12.11f ",V[i]);
    }
    printf("\n");
}

void Print_Couplings(vector<struct node> &J) {
	printf("SIZE %d\n",J.size());
    for(int i = 0; i < 30; i++) {
      	printf("%12.11f ",J[i].J);
    }
    printf("\n");
}