void Evaluate_Observables(variables *v, parameters *p) {

	v->q12=Prod(v->S0,v->S2);
	v->q13=Prod(v->S3,v->S2);
	v->NPG2 = Norm2(v->PG2);
	v->PGPG = Overlap(v->PG2,sqrt(v->NPG2),v->PG1,sqrt(v->NPG1));
}

double Total_Energy(vector<double> S1, variables *v, parameters *p) {
	v->e2=Energy(2,v->J2,S1); v->e3=Energy(3,v->J3,S1);
	return p->a2*v->e2+p->a3*v->e3;
}

double Total_Energy_W(vector<double> S1, variables *v, parameters *p) {

	v->TPG = Total_Projected_Gradient(S1,v,p);
	return Norm2(v->TPG);
}

vector<double> Total_Gradient(vector<double> S1, variables *v, parameters *p) {

	return Lin(p->a2,Gradient(2,v->J2,S1),p->a3,Gradient(3,v->J3,S1));
}

vector<double> Total_Hessian(vector<double> S1, variables *v, parameters *p) {

	return Lin(p->a2,Hessian(2,v->J2,S1),p->a3,Hessian(3,v->J3,S1));
}

vector<double> Total_Projected_Gradient(vector<double> S1, variables *v, parameters *p) {

	v->G1 = Total_Gradient(S1,v,p);
	v->Mu1 = -Prod(S1,v->G1);
	double NS1 = Norm2(S1);

	return Lin(1.,v->G1,v->Mu1/NS1,S1);
}

vector<double> Total_Projected_Gradient_W(vector<double> S1, variables *v, parameters *p) {

	v->H2 = Total_Hessian(S1,v,p);
	v->G1 = Multiply(v->H2,v->TPG);
	v->Mu1 = -Prod(S1,v->G1);
	
	return Lin(1.,v->G1,v->Mu1/p->N,S1);
}

vector<double> Projected_Noise(vector<double> S1, variables *v, parameters *p, double Temp) {

	static mt19937 gauss(p->seedX);

	v->X = White_Noise(p->N,gauss,sqrt(2*Temp/v->DTime)); 
	v->MuX = -Prod(S1,v->X);

	return Lin(1.,v->X,v->MuX/p->N,S1);
}

vector<double> Magnetization(vector<double> S1, vector<double> M, variables *v, parameters *p) {


	return Lin(1.-1./p->R, M, 1./p->R, S1);
}

vector<double> Gradient_TAP(variables *v, parameters *p) {

	double q = Norm2(v->M)/p->N;
	vector<double> GradM = Total_Gradient(v->M, v, p);
	return Lin(-1,GradM,-p->Temp*R1_TAP(p->a2,p->a3,q,1./p->Temp),v->M);
}

