struct variables {

	double alpha;
	double alphaCG;

	vector<struct node> J2; 
    vector<struct node> J3; 

    vector<double> S1,S2,S0,Stemp1,S3;
    vector<double> G1;
    vector<double> PG0,PG1,PG2;
    vector<double> TPG;
    vector<double> RPG0,RPG1;
    vector<double> CG0,CG1,RCG0,RCG1;

	/*#ifdef LAPACKEE
	double *Hc;
	Hc = (double *) calloc(N*N, sizeof(double));
    double *wc = (double *) calloc(N,sizeof(double));
    double *S2c = (double *) calloc(N,sizeof(double));
    double *mE0c = (double *) calloc(N,sizeof(double));
	#endif*/

    //default_random_engine (seedX);

	vector<double> EN; //EN(100,0.);

    vector<double> X,PX; //X(N, 0.0);

	double e2,e3;
    double E0,E1,E2;
    double Mu1,MuX;

    double T;

    double Time,Time_prev;
	double DTime;
    double PGPG,CGCG;
    double NPG1,NPG2;

    double q12;

    vector<double> H1,H2,H3; //H1(N),H2(N),DH(N);
    vector<double> P;

    vector<double> Ei;	//mE0(N,1.);
    vector<double> gEi; //mE1(N,1.);

    /*double eigen0,eigen00,eigen1,eigen11;
	double DirLaloux,DirS,DirS0,NPG;
    double HighestEigen = 2*sqrt((a2*a2*2+a3*a3*3*2*1.)/2.);*/
};

