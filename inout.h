void Create_Directory(char *dir) {
    struct stat st = {0};
    if (stat(dir, &st) == -1) {
        mkdir(dir, 0700);
    }
}

void SnapVector(vector<double> V, char *fname) {
    FILE *pFile;
    pFile = fopen (fname, "wb");
    fwrite (&V[0], sizeof(double), V.size(), pFile);
    fclose (pFile);
}


void SaveVector(variables *v, parameters *p, vector<double> V) {

    char Vsave[200];
    sprintf(Vsave,"%s/En%1.12f.v",p->dir,v->E2/p->N);

    FILE *pFile;
    pFile = fopen (Vsave, "wb");
    fwrite (&V[0], sizeof(double), V.size(), pFile);
    fclose (pFile);
}


vector<double> OpenVector(int N,char *fname) {

    vector<double> X(N);

    if( access( fname, F_OK ) != -1 ) {
        printf("OPENING CONFIG %s\n",fname);

        FILE *fo;

        fo = fopen(fname,"rb");
    

        int result = fread (&X[0],sizeof(double),N,fo);
        if (result != N) {fputs ("Reading error",stderr); exit (3);}

        fclose(fo);

    }
    else {
        printf("FILE %s DOES not EXIST\n",fname);
    }

    return X;
}

void SaveSystem(variables *v, parameters *p) {

    char Ssave[200];
    sprintf(Ssave,"%s/En%1.12f.sy",p->dir,v->E2/p->N);

    FILE *pFile;
    pFile = fopen (Ssave, "w");

    fprintf(pFile,"a2 %lf\n",p->a2);
    fprintf(pFile,"a3 %lf\n",p->a3);

    fprintf(pFile,"N %d\n",p->N);

    fprintf(pFile,"seedJ %d\n",p->seedJ);
    fprintf(pFile,"seedS %d\n",p->seedS);
    fprintf(pFile,"seedX %d\n",p->seedX);
    
    fprintf(pFile,"Temp %lf\n",p->Temp);
    fprintf(pFile,"T %lf\n",v->T);
    fprintf(pFile,"Beta %lf\n",p->Beta);

    fprintf(pFile,"DT0 %lf\n",p->DT0);
    fprintf(pFile,"ITime %lf\n",v->Time);
    fprintf(pFile,"TTime %lf\n",p->TTime);

    fprintf(pFile,"TVel %lf\n",p->TVel);

    fprintf(pFile,"discreteJ2 %d\n",p->discreteJ2);
    fprintf(pFile,"discreteJ3 %d\n",p->discreteJ3);

    fprintf(pFile,"dilution %d\n",p->dilution);

    fprintf(pFile,"alpha0 %lf\n",p->alpha0);
    fprintf(pFile,"alpha0_CG %lf\n",p->alpha0_CG);

    fprintf(pFile,"PGPGpar %lf\n",p->PGPGpar);

   fclose (pFile);
}

void OpenSystem(variables *v, parameters *p, char *fname) {
    FILE *pFile;
    if((pFile = fopen (fname, "r"))==NULL) {
        printf("Errore nell'apertura del file'");
        exit(1);
    }

    printf("OPENING SYSTEM %s\n",fname);

    char op[120];
    int b;

    fscanf(pFile,"%s %lf\n",op,&p->a2);   if(strcmp(op,"a2")) { cout << "ERROR: not a2 "; }
    fscanf(pFile,"%s %lf\n",op,&p->a3);   if(strcmp(op,"a3")) { cout << "ERROR: not a3 "; }

    fscanf(pFile,"%s %d\n",op,&p->N);    if(strcmp(op,"N")) { cout << "ERROR: not N "; p->sqrtN = p->N; }

    fscanf(pFile,"%s %d\n",op,&p->seedJ);    if(strcmp(op,"seedJ")) { cout << "ERROR: not seedJ "; }
    fscanf(pFile,"%s %d\n",op,&p->seedS);    if(strcmp(op,"seedS")) { cout << "ERROR: not seedS "; }
    fscanf(pFile,"%s %d\n",op,&p->seedX);    if(strcmp(op,"seedX")) { cout << "ERROR: not seedX "; }
    
    fscanf(pFile,"%s %lf\n",op,&p->Temp);     if(strcmp(op,"Temp")) { cout << "ERROR: not Temp "; }
    fscanf(pFile,"%s %lf\n",op,&v->T);        if(strcmp(op,"T")) { cout << "ERROR: not T "; }
    fscanf(pFile,"%s %lf\n",op,&p->Beta);     if(strcmp(op,"Beta")) { cout << "ERROR: not Beta "; }

    fscanf(pFile,"%s %lf\n",op,&p->DT0);      if(strcmp(op,"DT0")) { cout << "ERROR: not DT0 "; }
    fscanf(pFile,"%s %lf\n",op,&p->ITime);    if(strcmp(op,"ITime")) { cout << "ERROR: not ITime "; }
    fscanf(pFile,"%s %lf\n",op,&p->TTime);    if(strcmp(op,"TTime")) { cout << "ERROR: not TTime "; }

    fscanf(pFile,"%s %lf\n",op,&p->TVel);     if(strcmp(op,"TVel")) { cout << "ERROR: not TVel "; }

    fscanf(pFile,"%s %d\n",op,&b);  p->discreteJ2=b;  if(strcmp(op,"discreteJ2")) { cout << "ERROR: not discreteJ2 "; }
    fscanf(pFile,"%s %d\n",op,&b);  p->discreteJ3=b;  if(strcmp(op,"discreteJ3")) { cout << "ERROR: not discreteJ3 "; }

    fscanf(pFile,"%s %d\n",op,&p->dilution);  if(strcmp(op,"dilution")) { cout << "ERROR: not dilution "; }

    fscanf(pFile,"%s %lf\n",op,&p->alpha0);    if(strcmp(op,"alpha0")) { cout << "ERROR: not alpha0 "; }
    fscanf(pFile,"%s %lf\n",op,&p->alpha0_CG); if(strcmp(op,"alpha0_CG")) { cout << "ERROR: not alpha0_CG "; }

    fscanf(pFile,"%s %lf\n",op,&p->PGPGpar);   if(strcmp(op,"PGPGpar")) { cout << "ERROR: not PGPGpar "; }

    fclose (pFile);
}