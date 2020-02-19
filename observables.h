void Evaluate_Observables(variables *v, parameters *p) {

	v->q12=Prod(v->S0,v->S2);
	v->NPG2 = Norm2(v->PG2);
	v->PGPG = Overlap(v->PG2,sqrt(v->NPG2),v->PG1,sqrt(v->NPG1));
}

double Total_Energy(vector<double> S1, variables *v, parameters *p) {
	v->e2=Energy(2,v->J2,S1); v->e3=Energy(3,v->J3,S1);
	return p->a2*v->e2+p->a3*v->e3;
}

vector<double> Total_Gradient(vector<double> S1, variables *v, parameters *p) {

	return Lin(p->a2,Gradient(2,v->J2,S1),p->a3,Gradient(3,v->J3,S1));
}

vector<double> Total_Projected_Gradient(vector<double> S1, variables *v, parameters *p) {

	v->G1 = Lin(p->a2,Gradient(2,v->J2,S1),p->a3,Gradient(3,v->J3,S1));
	v->Mu1 = -Prod(S1,v->G1);

	return Lin(1.,v->G1,v->Mu1/p->N,S1);
}

vector<double> Projected_Noise(vector<double> S1, variables *v, parameters *p, double Temp) {

	static mt19937 gauss(p->seedX);

	v->X = White_Noise(p->N,gauss,sqrt(2*Temp/v->DTime)); 
	v->MuX = -Prod(S1,v->X);

	return Lin(1.,v->X,v->MuX/p->N,S1);
}
