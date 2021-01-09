vector<double> Find_Inherent_Structure(vector<double> S1,variables *v, parameters *p) {

	//v.Stemp1 = v.S1;

	struct variables v_temp;

	v_temp.J2=v->J2;
	v_temp.J3=v->J3;

	v_temp.S1 = v->S1;

	v_temp.E1 = p->TotEN(v_temp.S1,&v_temp,p);
	v_temp.PG1 = p->TotProjGR(v_temp.S1,&v_temp,p);

	Initialize_CG(&v_temp, p);


	while(Norm2(v_temp.TPG) > p->N*0.00000001 && fabs(v_temp.alphaCG)>0.000000001) {
		v_temp.S2 = CG_Step(v_temp.S1,v_temp.PG1,&v_temp,p);

		v_temp.E2 = p->TotEN(v_temp.S2,&v_temp,p);
		v_temp.PG2 =p->TotProjGR(v_temp.S2,&v_temp,p);

		Evaluate_Observables(&v_temp,p);
	
		v_temp.E0=v_temp.E1;

		v_temp.S1=v_temp.S2; v_temp.E1=v_temp.E2; v_temp.PG1=v_temp.PG2; v_temp.NPG1=v_temp.NPG2;
		v_temp.CG0=v_temp.CG1; v_temp.H1=v_temp.H2; v_temp.RPG0=v_temp.RPG1; v_temp.RCG0=v_temp.RCG1;

		//cout << "energy" << v_temp.E2 << " angle" << v_temp.alphaCG << endl; getchar();	
	}

	#ifdef LAPACKEE
		v_temp.G1 = Total_Gradient(v_temp.S1,&v_temp,p);
		v_temp.Mu1 = -Prod(v_temp.S1,v_temp.G1);
		v_temp.H1=Total_Hessian(v_temp.S1,&v_temp,p);
		Evaluate_eigenvalues(v_temp.H1, v_temp.PG1, v_temp.S1, p->N, v_temp.Ei, v_temp.gEi);
		vector<double> VC(p->N,v_temp.Mu1/p->N);
		v_temp.Ei=Sum(v_temp.Ei,VC);

		printf("%f %f \n",Norm2(v_temp.gEi),v_temp.Ei[0]);
	#else
		printf("LAPACKEE not activated\n");
	#endif

	cout << v->Time << ' ' << v->E2 
	<< ' ' << Norm2(v_temp.PG1)/p->N
	<< ' ' << Norm2(v_temp.TPG)/p->N
	<< ' ' << Prod(v_temp.S1,v->S0)/p->N
	<< ' ' << Prod(v_temp.S1,v->S1)/p->N;

	//v.S1 = v.Stemp1;

	return v_temp.S1;
}