struct parameters {

	double err,err2;

	double a2,a3,a4;

	double TTime = 1000;

	int N;
	
	double Beta;
	double Temp;

	double alpha;
	double PGPGpar;

	double alpha0;

	int seedJ = 1;
	int seedS = 3;
	int seedX = 4;

	char method[4];
	char dilute;

	FILE *fout,*fout2;
	char fname[120];
	char fname2[120];

	char dir[120];
}

void Input_Parameters(int Ac, char **Av, struct parameters *PAR) {

	PAR->err=0.000001;
    PAR->err2=pow(err,2);

	PAR->a2 = 1.;
    PAR->a3 = 1.;
    PAR->a4 = 1.;

	PAR->N = 1000;
	PAR->TTime = 1000;
	
	PAR->Beta = 0.;
	PAR->Temp = 0.;

	PAR->alpha=0.0015;
	PAR->PGPGpar=0.999995;

	PAR->alpha0=0.1;


	PAR->seedJ = 1;
	PAR->seedS = 3;
	PAR->seedX = 4;


    char *cvalue = NULL;
    int index;
    int c;
    opterr = 0;

    printf("to see input parameters look at input.h!\n");

    while ((c = getopt (Ac, Av, "N:J:S:b:T:v:D:c:A:M:2:3:4:")) != -1) //COLON MEANS THAT A VALUE MUST BE SPECIFIED
    switch (c)
      {
      case 'N': //NUMBER OF SPINS
        cvalue = optarg;
        PAR->N=atoi(cvalue);        
        break;
      case 'J': //seedJ OF THE COUPLINGS
        cvalue = optarg;
        PAR->seedJ=atoi(cvalue);
        break;
      case 'S': //seedX OF THE CONFIGURATION
        cvalue = optarg;
        PAR->seedX=atoi(cvalue);
        break;
      case 'b': //PLANTING TEMPERATURE
        cvalue = optarg;
        PAR->Beta=atof(cvalue);
        break;
      case 'T': //MONTECARLO TEMPERATURE
        cvalue = optarg;
        PAR->Temp=atof(cvalue);
        break;
      case 'v': //ANNEALING VELOCITY IN beta  
        cvalue = optarg;
        PAR->bVel = atof(cvalue);
        break;
      case 'D': //Diluition  
        cvalue = optarg;
        PAR->D = *cvalue;
        break;
      case 'M': //METHOD OF MINIMIZATION
        cvalue = optarg;
        sprintf(PAR->method,"%s",cvalue);
        break;
      case '2': //weight ot the 2-spin coplings
        cvalue = optarg;
        PAR->a2 = atof(cvalue);
        break;
      case '3': //weight ot the 3-spin coplings
        cvalue = optarg;
        PAR->a3 = atof(cvalue);
        break;
      case '4': //weight ot the 4-spin coplings
        cvalue = optarg;
        PAR->a4 = atof(cvalue);
        break;
      case '?':
        fprintf (stderr,"Unknown option character %c.\n",optopt);
        getchar();
        break;
        default:
        abort ();
      }

    if(dilute=='y') {
		D2=3./N;
		D3=3./N;
		//D2=2000./N;
		//D3=2000./N/N;
	} else { 
		D2=1;
		D3=1;
		dilute='n'; 
	}

	sprintf(fname,"%s_dil%c_%g.2+%g.3-spin_T%g_Tp%g_N%d_sJ%dsS%dsX%d_alpha%g.dat",method,dilute,a2,a3,Temp,1./Beta,N,seedJ,seedS,seedX,alpha);
    sprintf(fname2,"%seigen_dil%c_%g.2+%g.3-spin_T%g_Tp%g_N%d_sJ%dsS%dsX%d_alpha%g.dat",method,dilute,a2,a3,Temp,1./Beta,N,seedJ,seedS,seedX,alpha);
   

  //Default output
  /*PAR->file=NULL;
  int i;
  for(i=0;i<5;i++) { sprintf(PAR->dir[i],"./"); }*/
}
