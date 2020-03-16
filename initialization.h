vector<double> Initialize_System(variables *v, parameters *p) {

	v->S0 = RandomSpin(p->seedS,p->N);

	if(p->a2>0.) {
		/*if(p->dilute=='y') { Initialize_Dilute_J(2,v->J2,p->seedJ,p->N,p->D2,true); }
		else if(p->dilute=='b') { Initialize_Dilute_Boxed_J2(v->J2,p->seedJ,p->N,p->D2,true); }
		else { 		*/Initialize_All_J(2,v->J2,p->seedJ,p->N); 	//}
		cout << "#normJ2" << sqrt(Norm2_J(v->J2)) << endl;
	}
	if(p->a3>0.) {
		/*if(p->dilute=='y') { Initialize_Dilute_J(3,v->J3,p->seedJ,p->N,p->D3,true); }
		else if(p->dilute=='b') {*/ Initialize_Dilute_Boxed_J3(v->J3,p->seedJ,p->N,p->D3,true); //}
		/*else { 		Initialize_All_J(3,v->J3,p->seedJ,p->N); 	}*/
		cout << "#normJ3" << sqrt(Norm2_J(v->J3)) << endl;
	}

	if(p->Beta!=0) { 
		Planting(2,p->Beta,v->J2,v->S0);
		Planting(3,p->Beta,v->J3,v->S0);
	}

	//cout << "ES0 before" << Total_Energy(v->S0,v,p) << endl;  getchar();

	if(strcmp(p->finitial,"none")) { char Vopen[120]; sprintf(Vopen,"%s.v",p->finitial); v->S0=OpenVector(p->N,Vopen); /*Print_Vector(v->S0); Print_Couplings(v->J3);*/ }

	//cout << "ES0 after" << Total_Energy(v->S0,v,p) << endl;  getchar();

	#ifdef HESSIAN
	H2 = Lin(a2,Hessian(2,J2,S2),a3,Hessian(3,J3,S2));		
	P = Id(N);
    eigen0 = MinEigen(H2,P,mE0,HighestEigen,err2);
	P = Proj(S2);	
    eigen1 = MinEigen(H2,P,mE1,HighestEigen,err2);

	eigen00 = Sandwich(H2,mE0);
	eigen11 = Sandwich(H2,mE1);
	DirLaloux = Sandwich(H2,PG1);
	DirS = Sandwich(H2,S2);
	DirS0 = Sandwich(H2,S0);
	#endif

	return v->S0;
}
