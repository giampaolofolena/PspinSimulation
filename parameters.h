struct parameters {

	double err,err2;

	double a2,a3;

	double DT0;
	double ITime;
	double TTime;

	int N;
	double sqrtN;
	
	double Beta;
	double Temp;

	double Tin;
	double Tfin;
	double TVel;

	double alpha0;
	double alpha0_CG;
	double PGPGpar;

	char potential;

	int seedJ;
	int seedS;
	int seedX;

	bool discreteJ2;
	bool discreteJ3;

	char method[4];
	int dilution;
	char HessianEi;

	FILE *fout,*foutEi;
	char fname[200];
	char fnameEi[200];

	FILE *fin;
	char finitial[200];
	char fpar[200];

	char dir[120];

	double D2,D3,D4;


    double (*TotEN)(vector<double> S1, variables *v, parameters *p);
    vector<double> (*TotProjGR)(vector<double> S1, variables *v, parameters *p);

};

void Input_Parameters(int Ac, char **Av, struct variables *v, struct parameters *p) {

	p->err=0.000001;
	p->err2=pow(p->err,2);

	p->a2 = 1.;
	p->a3 = 1.;

	p->N = 1000;
	p->sqrtN = sqrt(p->N);

	p->DT0 = 0.01;
	p->ITime = 0;
	p->TTime = 10000;
	
	p->Beta = 0.;
	p->Temp = 0.;
	double pTemp = -999;

	p->Tin = 0;
	p->Tfin = 0;
	p->TVel = 0;

	p->PGPGpar = 0.999995;

	p->potential = 'e';

	p->discreteJ2 = false;
	p->discreteJ3 = true;

	p->dilution=3;

	p->HessianEi='n';

	p->alpha0=0.0015;
	p->alpha0_CG = 0.015;

	p->seedJ = 1;
	p->seedS = 3;
	p->seedX = 4;

	sprintf(p->method,"%s","GD");

	sprintf(p->finitial,"%s","none");
	sprintf(p->fpar,"%s","none");

    char *cvalue = NULL;
    int index;
    int c;
    opterr = 0;

    printf("to see input parameters look at input.h!\n");

    while ((c = getopt (Ac, Av, "N:J:S:X:b:T:I:v:F:D:c:A:M:t:2:3:i:s:P:d:H:t")) != -1) //COLON MEANS THAT A VALUE MUST BE SPECIFIED
    switch (c)
      {
      case 'N': //NUMBER OF SPINS
        cvalue = optarg;
        p->N=atoi(cvalue); 	if(p->N>pow(2,16)) { printf("ACTUNG: maximal size of the system exceeded! look at initialization of J"); getchar(); }
        p->sqrtN = sqrt(p->N);      
        break;
      case 'J': //seedJ OF THE COUPLINGS
        cvalue = optarg;
        p->seedJ=atoi(cvalue);
        break;
      case 'S': //seedS OF THE CONFIGURATION
        cvalue = optarg;
        p->seedS=atoi(cvalue);
        break;
      case 'X': //seedX OF THE NOISE
        cvalue = optarg;
        p->seedX=atoi(cvalue);
        break;
      case 'b': //PLANTING TEMPERATURE
        cvalue = optarg;
        p->Beta=atof(cvalue);
        break;
      case 'T': //MONTECARLO TEMPERATURE
        cvalue = optarg;
        p->Temp=atof(cvalue);
        pTemp = p->Temp;
        break;
	  case 'I': //ANNEALING Initial Temperature  
        cvalue = optarg;
        p->Tin = atof(cvalue);
        break;
      case 'v': //ANNEALING VELOCITY IN beta  
        cvalue = optarg;
        p->TVel = atof(cvalue);
        break;
	  case 'F': //ANNEALING Final Temperature  
        cvalue = optarg;
        p->Tfin = atof(cvalue);
        break;
      case 'D': //Diluition if p->dilution=88888888 >> p->dilute='b'
        cvalue = optarg;
        p->dilution = atoi(cvalue);
        break;
      case 'M': //METHOD OF MINIMIZATION
        cvalue = optarg;
        sprintf(p->method,"%s",cvalue);
        break;
      case 't': //Total Time
        cvalue = optarg;
        p->ITime=atof(cvalue);
        break;
      case '2': //weight ot the 2-spin coplings
        cvalue = optarg;
        p->a2 = atof(cvalue);
        break;
      case '3': //weight ot the 3-spin coplings
        cvalue = optarg;
        p->a3 = atof(cvalue);
        break;
      case 'i': //initial configuration and parameters from file
        cvalue = optarg;
        sprintf(p->finitial,"%s",cvalue);
        break;
	  case 's': //system parameters from file
        cvalue = optarg;
        sprintf(p->fpar,"%s",cvalue);
        break;
	  case 'P': //system parameters from file
        cvalue = optarg;
        p->potential = *cvalue;
        break;
      case 'd': //discrete interactions
        cvalue = optarg;
        p->discreteJ3 = (bool)*cvalue;
        break;
      case 'H': //Hessian
        cvalue = optarg;
        p->HessianEi = *cvalue;
        break;
      case '?':
        fprintf (stderr,"Unknown option character %c.\n",optopt);
        getchar();
        break;
        default:
        abort ();
      }

	if(p->Tin==0 && p->Temp!=0) {
		p->Tin = p->Temp;
		p->Tfin = p->Temp;
		p->TVel = 0.;
	}

	if(p->potential=='e') { p->TotEN = Total_Energy; p->TotProjGR = Total_Projected_Gradient; cout << "HAMILTONIAN dynamics" << endl; } 
	else if (p->potential=='w') { p->TotEN = Total_Energy_W; p->TotProjGR = Total_Projected_Gradient_W; cout << "SADDLE-POTENTIAL dynamics" << endl; }
	else { cout << "NO POTENTIAL DEFINED" << endl; }

	// OPEN CONFIGURATION FROM FILE AND IF IT EXISTS OPEN FILE SYSTEM PARAMETERS
	if(strcmp(p->finitial,"none")) { 
		char Sopen[120]; 
		int l=strlen(p->finitial); p->finitial[l-2]=0; sprintf(Sopen,"%s.sy",p->finitial); 

		if(access( Sopen, F_OK ) != -1) { OpenSystem(v,p,Sopen); }
	}

	// OPEN FILE SYSTEM PARAMETERS
	if(strcmp(p->fpar,"none")) {
		char Sopen[120]; 
		sprintf(Sopen,"%s",p->fpar); 

		if(access( Sopen, F_OK ) != -1) { OpenSystem(v,p,Sopen); }
	}

	//STANDARD DILUTION
	p->D3=p->dilution/(double)p->N;
	p->D2=1; //ATTENZIONE: se scalo D2 come D3 va a finire ad energie molto sotto soglia? perché? meglio fissare D2=1


	sprintf(p->dir,"P%g,2+%g,3_N%d_sJ%d_D%d",p->a2,p->a3,p->N,p->seedJ,p->dilution);
	Create_Directory(p->dir);

	
	if(!strcmp(p->finitial,"none")) {
		sprintf(p->fname,"%s/%s_T%g_Tp%g_S%d_X%d_%c.dat",p->dir,p->method,p->Temp,1./p->Beta,p->seedS,p->seedX,p->potential);
		p->fout=fopen(p->fname,"a");
	} else {
		sprintf(p->fname,"%s_%s%c.dat",p->finitial,p->method,p->potential);
		p->fout=fopen(p->fname,"w");
	}

	if (p->fout == NULL) {
    	perror("Failed");
	}

	if(p->HessianEi!='n') {
	#ifdef LAPACKEE
	sprintf(p->fnameEi,"%s.eigen",p->fname);  
    p->foutEi=fopen(p->fnameEi,"w");
    #else
    printf("LAPACKEE not activated\n");
    #endif
	}
}
