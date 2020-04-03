//MODIFICHE DELL''ULTIMO MINUTO
//#include <Eigen/Dense>

double Sandwich(vector<double> &M, vector<double> &A);

double DoubleSandwich(vector<double> &A, vector<double> &M, vector<double> &B);

#ifdef LAPACKEE
//#include "/usr/include/gsl/gsl_cblas.h"  /// NEL CASO DI MAC metti solo "cblas.h" requires the link flag -lcblas
#include "cblas.h"
#include "lapacke.h" ///it requires the link flag -llapack

/*vector<double> FirstOrderEigenPert(vector<double> H, vector<double> L) {
    
    int N = l.size();
    vector<double> Dl(N);

    for(int i=0; i<N ; i++) {
        cblas_dtrmv(CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, const gsl_matrix * A, gsl_vector * x)
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, L, N, H, N, 0., HP, N);
        cblas_ddot(N,x,1,y,1);
    }
}*/




void Evaluate_eigenvalues(vector<double> &H, vector<double> G, vector<double> S, int N, vector<double> &W, vector<double> &VW) {
    
    double n2S=0;
    double n2G=0;
    //for(int i; i<N; i++) { n2X+=X[i]*X[i]; }

    double *P,*HH,*HP,*w,*x,*g,*L,*dH,*dHL;
    P = (double *) calloc(N*N, sizeof(double));
    HH = (double *) calloc(N*N, sizeof(double));
    HP= (double *) calloc(N*N, sizeof(double));
    /*L = (double *) calloc(N*N, sizeof(double));
    dH= (double *) calloc(N*N, sizeof(double));
    dHL= (double *) calloc(N*N, sizeof(double));*/
    w= (double *) calloc(N, sizeof(double));
    x= (double *) calloc(N, sizeof(double));
    g= (double *) calloc(N, sizeof(double));


    for(int i = 0; i < N*N; i++) {
        HH[i]=H[i];
    }

    /*or(int i=0;i<N;i++){
        P[N*i+i] = 1-X[i]*X[i]/n2X;
        for(int j=i+1;j<N;j++){
            P[N*i+j] = -X[i]*X[j]/n2X;
            P[N*j+i] = P[N*i+j];
        }
    }*/

    for(int i=0;i<N;i++){
        P[N*i+i]=1;
        x[i]=S[i];
        n2S+=S[i]*S[i];
        g[i]=G[i];
        n2G+=G[i]*G[i];
    }

    cblas_dsyr(CblasRowMajor,CblasUpper, N, -1./n2S, x, 1, P, N);

//( __Order, __TransA, __TransB, __M, __N, __K, __alpha, *__A, __lda, *__B, __ldb, __beta, *__C, __ldc)
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, HH, N, P, N, 0., HP, N);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, P, N, HP, N, 0., HH, N);



    int info=0;
    int lda = N;

    // Solve eigenproblem 
    info = LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V', 'U', N, HH, lda, w );
    // Check for convergence 
    if( info > 0 ) {
        printf( "The algorithm failed to compute eigenvalues.\n" );
        exit( 1 );
    }

    cblas_dgemv(CblasRowMajor,CblasNoTrans,N,N,1.,HH,N,g,1,0.0000000000000001,x,1);


    n2G=sqrt(n2G);

    if(W.size()<N) { W.resize(N,0); }
    if(VW.size()<N) { VW.resize(N,0); }


    for(int i = 0; i < N; i++) {
        
        /*double d = 0;
        for(int j = 0; j < N; j++) { d+=HH[N*i+j]*HH[N*i+j]; }
        cout << d << ' '; getchar();*/

        W[i] = w[i];
        VW[i] = x[i]/n2G;
    }

    /*for(int i = 0; i < N*N; i++) {
        L[i]=HH[i];
        dH[i]=DH[i];
    }

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1.0, dH, N, L, N, 0., dHL, N);

    for(int i=0; i<N ; i++) {
        //cblas_dtrmv(CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, const gsl_matrix * A, gsl_vector * x)
        VW[i] = 0;
        for(int j=0; j<N ; j++) {
            VW[i] += dHL[j*N+i]*L[i*N+j];
        }
    }

    vector<double> l(N);
    vector<double> l0(N);


    for(int i=0; i<N ; i++) {
        l0[i] = L[i];
    }

    VW[0] = DoubleSandwich(l0,DH,l0);


    VW[1]=0;
    for(int i=1; i<N ; i++) {
        for(int j = 0; j < N; j++) {
            l[j]=L[N*i+j];
        }
        VW[1] += pow(DoubleSandwich(l0,DH,l),2)/(w[0]-w[i]);
    }*/
}

double Norm22(const double * X, const int n){
    double norm2=0.;
    int i;
    for(i=0;i<n;i++) norm2 += X[i]*X[i];
    return norm2;
}
void Projj(double * H, const double * X, const int n) {

    double n2X = Norm22(X,n);
    int i,j,k;

    double *P;
    P = (double *) calloc(n*n, sizeof(double));

    for(i=0;i<n;i++){
        P[n*i+i] = 1-X[i]*X[i]/n2X;
        for(j=i+1;j<n;j++){
            P[n*i+j] = -X[i]*X[j]/n2X;
            P[n*j+i] = P[n*i+j];
        }
    }

    double *HP;
    HP = (double *) calloc(n*n, sizeof(double));

    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            HP[n*i+j]=0.;
            for(k=0;k<n;k++){
                HP[n*i+j] += H[n*i+k]*P[n*k+j];
            }
        }
    }

    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            H[n*i+j]=0.;
            for(k=0;k<n;k++){
                H[n*i+j] += P[n*i+k]*HP[n*k+j];
            }
        }
    }
    free(HP);
    free(P);
}

//FINE MODIFICHE DELL?ULTIMO MINUTO
#endif

#include <random>  
#define U_pot 10000000


struct node {
    double J;
    long long int SSS;
};

double Prod(vector<double> A, vector<double> B) {
    double C=0.;
    int N = A.size();

    for(int i=0; i<N; i++) {  C+=A[i]*B[i]; }
    /*double *a = (double *) calloc(N, sizeof(double));
    double *b = (double *) calloc(N, sizeof(double));
    for(int i=0; i<N; i++) {  a[i]=A[i];
                                b[i]=B[i]; }
    C=cblas_ddot(N,a,1,b,1);*/

    return C;
}

double Mean(vector<double> A) {
    double C=0.;
    int N = A.size();
    for(int i=0; i<N; i++) {  C+=A[i]; }
    return C/N;
}

vector<double> Scal(vector<double> A, double c) {
    int N = A.size();
    vector<double> C(N);
    for(int i=0; i<N; i++) { C[i]=c*A[i]; }
    return C;
}

vector<double> Sum(vector<double> A, vector<double> B) {
    int N = A.size();
    vector<double> C(N);    
    for(int i=0; i<N; i++) { C[i]=A[i]+B[i]; }
    return C;
}

vector<double> Sub(vector<double> A, vector<double> B) {
    int N = A.size();
    vector<double> C(N);    
    for(int i=0; i<N; i++) { C[i]=A[i]-B[i]; }
    return C;
}

vector<double> Lin(double a, vector<double> A, double b, vector<double> B) {
    int N = A.size();
    vector<double> C(N);
    for(int i=0; i<N; i++) { C[i]=a*A[i]+b*B[i]; }
    return C;
}

vector<double> ProdMat(vector<double> &M, vector<double> &A) {
    int N = A.size();

    vector<double> B(N, 0.0);

    for(int i=0; i<N; i++) {
        for(int k=0; k<N; k++) {
            B[i] += M[i*N+k]*A[k];
        }
    }
    return B;
}

double Sandwich(vector<double> &M, vector<double> &A) {
    int N = A.size();

    double B=0;
    double NormA=0;

    for(int i=0; i<N; i++) {
        NormA+=A[i]*A[i];
        for(int k=0; k<N; k++) {
            B += M[i*N+k]*A[k]*A[i];
        }
    }
    return B/NormA;
}

double DoubleSandwich(vector<double> &A, vector<double> &M, vector<double> &B) {
    int N = A.size();

    double C=0;

    for(int i=0; i<N; i++) {
        for(int k=0; k<N; k++) {
            C += M[i*N+k]*B[k]*A[i];
        }
    }
    return C;
}

double Norm2_J(vector<struct node> &J) {
    double norm2=0.;
    for(size_t i = 0, size = J.size(); i < size; i++) {
        norm2 += J[i].J*J[i].J;
    }
    return norm2;
}

double Norm2(vector<double> V) {
    int N = V.size();

    double norm2=0.;
    for(size_t i=0; i<N; i++) 
    { 
        norm2 += V[i]*V[i];
    }
    return norm2;
}

vector<double> RotateVector(vector<double> S, vector<double> PG, double angle) {


    double SPG = Prod(S,PG);
    if(SPG>0.0000000001) { cout << "ERR_not_ort " << SPG << endl; }

    int N = S.size();
    vector<double> RS(N);

    double cs = cos(angle);
    double sn = sin(angle);
    double normS = sqrt(Norm2(S));
    double normPG = sqrt(Norm2(PG));
    double cs_x=cs;//normS;
    double sn_x=sn*normS/normPG;
    
    for(int i=0; i<N; i++) {
        RS[i] = S[i]*cs_x-PG[i]*sn_x;//(S[i]*cs_x-PG[i]*sn_x)*normS;
    }

    return RS;


}

double Overlap(vector<double> &V, double normV, vector<double> &W, double normW) {
    double sc=0.;
    int N = V.size();
    int ni;
    for (ni=0; ni<N; ni++) {
        sc += V[ni]/normV*W[ni]/normW;
    }
    return sc;
}

vector<double> RotateGeneralVector(vector<double> X, vector<double> A, vector<double> B, double angle) {
    int i;

    int N = X.size();

    double cs = cos(angle);
    double sn = sin(angle);
    double normA = sqrt(Norm2(A));
    double normB = sqrt(Norm2(B));

    double XA=0.,XB=0.,YA=0.,YB=0.;
    vector<double> a(N),b(N);

    for(i=0;i<N;i++){
        a[i]=A[i]/normA;
        b[i]=B[i]/normB;
    }

    double ab = Prod(a,b);
    if(fabs(ab)>0.000001) { printf("a and b not orthogonal!! ab = %f \n",ab); }

    for(i=0;i<N;i++){
        XA += X[i]*a[i];
        XB += X[i]*b[i];
    }

    YA = XA*cs+XB*sn;
    YB = -XA*sn+XB*cs;
    
    vector<double> Y(N);

    for(i=0;i<N;i++){
        Y[i] = X[i]+(YA-XA)*a[i]+(YB-XB)*b[i];
    }

    return Y;
}

void RenormalizeNodes(vector<struct node> &J, double oldVar, double newVar) {
        
    double rj = sqrt(oldVar/newVar);
    for(size_t i = 0, size = J.size(); i < size; i++) {
        J[i].J/=rj;
    }
}

void RenormalizeVector(vector<double> &V, double oldVar, double newVar) {
        
    int N = V.size();

    double rv = sqrt(oldVar/newVar);

    for(size_t i=0; i<N; i++) { V[i]/=rv; }
}

void Initialize_Dilute_J(int p, vector<struct node> &J, int seedJ, int N, double D, bool discrete) {

    mt19937 generator(seedJ);
    //uniform_real_distribution<double> unif(0,1);
    //bernoulli_distribution d(0.5);

    normal_distribution<double> distribution(0,1.0);

    //default_random_engine re;

    uniform_int_distribution<int> int_dist(0,N-1);

    double dNJ = D*N/p;
    for(int i=1; i<p; i++) { dNJ *= (N-i)/i; }
    int NJ = (int)dNJ;


    int nj;
    double norm2=0.;
    int rep=0;

    for(size_t nj=0;nj<NJ;nj++){
        struct node newJ;
        vector<long long int> newS;

        for(int i=0; i<p; i++) { newS.push_back(int_dist(generator)); }

        sort(newS.begin(), newS.end());
        auto it = std::unique(newS.begin(),newS.end());
        bool wasUnique = (  it == newS.end());

        /*long long int checkSSS=0;
        for(int i=0; i<p; i++) { checkSSS += newS[i]<<(16*i); }

        for(size_t ni=0;ni<nj;ni++){
            if(J[ni].SSS==checkSSS) { wasUnique=false; rep++; break; }
        }*/

        if(wasUnique) {

            //cout << newJ.S[0] << ' ' <<  newJ.S[1] << ' ' << newJ.S[2] << ' ' << endl;
            double JJ = distribution(generator);

            if(discrete==true) {
                if(JJ>0) {newJ.J=1;}
                if(JJ<0) {newJ.J=-1;}
            }
            else {
                newJ.J = JJ;
            }

            long long int newSSS=0;
            for(int i=0; i<p; i++) { newSSS += newS[i]<<(16*i);  }
           
            newJ.SSS = newSSS;

            J.push_back(newJ);
 
            norm2 += newJ.J*newJ.J;
        }
        else
        { nj--; }
    }

    //cout << "REPETITIONS: " << rep << endl;
    cout << "#" << p << "-Interactions: " << nj << ' ' << NJ << endl;

    RenormalizeNodes(J,norm2,N/2.);
}

void Initialize_Dilute_Boxed_J2(vector<struct node> &J, int seedJ, int N, double D, bool discrete) {

    mt19937 generator(seedJ);

    int p=2;

    normal_distribution<double> distribution(0,1.0);

    double dNJ = D*N/p;
    for(int i=1; i<p; i++) { dNJ *= (N-i)/i; }
    int NJ = (int)dNJ;

    int boxL = (int) pow(D,-1/2.);

    int bowL2=boxL*boxL;

    int boxN = N/boxL;

    int boxN2 = boxN*boxN;

    uniform_int_distribution<int> int_dist(0,boxL-1);


    cout << "boxL" << boxL << ' ' << pow(D,-1/2.) << ' ' << D << endl;
    cout << "boxN" << boxN << endl;


    int nj = 0;
    double norm2=0;

    int rep=0;

    for(int i=0; i<boxN; i++) {
        for(int j=0; j<i; j++) {

            struct node newJ;
                
            double JJ = distribution(generator);
            if(discrete==true) {
            if(JJ>0) {newJ.J=1;}
            if(JJ<0) {newJ.J=-1;}
            }
            else {
                newJ.J = JJ;
            }

            long long int newSSS=0;

            long long int X=i*boxL + int_dist(generator); newSSS += X<<(16*0);
            long long int Y=j*boxL + int_dist(generator); newSSS += Y<<(16*1);

                //cout << X << ' '<< Y << ' '<< Z << ' ' << endl; getchar(); 
           
            newJ.SSS = newSSS;

            J.push_back(newJ);  
 
            norm2 += newJ.J*newJ.J;

                /*for(size_t ni=0;ni<nj;ni++){
                    if(J[ni].SSS==newSSS) { rep++; cout << ni << endl;  break; }
                }*/
            nj++;
        }
    }


    cout << "REPETITIONS: " << rep << endl;
    cout << "Interactions: " << nj << ' ' << NJ << endl;

    RenormalizeNodes(J,norm2,N/2.);
}

void Initialize_Dilute_Boxed_J3(vector<struct node> &J, int seedJ, int N, double D, bool discrete) {

    mt19937 generator(seedJ);

    int p=3;

    normal_distribution<double> distribution(0,1.0);

    double dNJ = D*N/p;
    for(int i=1; i<p; i++) { dNJ *= (N-i)/i; }
    int NJ = (int)dNJ;

    int boxL = (int) pow(D,-1/3.);

    int bowL2=boxL*boxL;

    int boxN = N/boxL;

    int boxN2 = boxN*boxN;

    uniform_int_distribution<int> int_dist(0,boxL-1);


    cout << "boxL" << boxL << ' ' << pow(D,-1/3.) << ' ' << D << endl;
    cout << "boxN" << boxN << endl;


    int nj = 0;
    double norm2=0;

    int rep=0;

    for(int i=0; i<boxN; i++) {
        for(int j=0; j<i; j++) {
            for(int k=0; k<j; k++) {

                struct node newJ;
                
                double JJ = distribution(generator);
                if(discrete==true) {
                if(JJ>0) {newJ.J=1;}
                if(JJ<0) {newJ.J=-1;}
                }
                else {
                    newJ.J = JJ;
                }

                long long int newSSS=0;

                long long int X=i*boxL + int_dist(generator); newSSS += X<<(16*0);
                long long int Y=j*boxL + int_dist(generator); newSSS += Y<<(16*1);
                long long int Z=k*boxL + int_dist(generator); newSSS += Z<<(16*2);

                //cout << X << ' '<< Y << ' '<< Z << ' ' << endl; getchar(); 
           
                newJ.SSS = newSSS;

                J.push_back(newJ);

                
 
                norm2 += newJ.J*newJ.J;

                /*for(size_t ni=0;ni<nj;ni++){
                    if(J[ni].SSS==newSSS) { rep++; cout << ni << endl;  break; }
                }*/
                nj++;

            }
        }
    }


    cout << "REPETITIONS: " << rep << endl;
    cout << "Interactions: " << nj << ' ' << NJ << endl;

    RenormalizeNodes(J,norm2,N/2.);
}


void Initialize_All_J(int p, vector<struct node> &J, int seedJ, int N) {

    mt19937 generator(seedJ);
    normal_distribution<double> distribution(0,1.0);
    //default_random_engine re;

    long long int w[4];

    int NJ = N/p;
    for(int i=1; i<p; i++) { NJ *= (N-i)/i; }

    double norm2=0.;

    for(int i=0;i<N;i++){
        w[0]=i;
        if(p==1) { 
            struct node newJ;
            newJ.J = distribution(generator);
            newJ.SSS = 0.;
            for(int s=0; s<p; s++) { newJ.SSS += w[s]<<(16*s); }
            J.push_back(newJ);
            norm2 += newJ.J*newJ.J;
        } else {

            for(int j=i+1;j<N;j++){
                w[1]=j;
                if(p==2) { 
                    struct node newJ;
                    newJ.J = distribution(generator);
                    newJ.SSS = 0.;
                    for(int s=0; s<p; s++) { newJ.SSS += w[s]<<(16*s); }
                    J.push_back(newJ);
                    norm2 += newJ.J*newJ.J;

                } else {

                    for(int k=j+1;k<N;k++){
                        w[2]=k;
                        if(p==3) { 
                            struct node newJ;
                            newJ.J = distribution(generator);
                            newJ.SSS = 0.;
                            for(int s=0; s<p; s++) { newJ.SSS += w[s]<<(16*s); }
                            J.push_back(newJ);
                            norm2 += newJ.J*newJ.J;

                            /*long long int mask = 0b1111111111111111;
                            for(int s=0; s<p; s++) { cout << (newJ.SSS>>(16*s) & mask) << ' '; }
                            getchar();*/
                        } else {

                            for(int l=k+1;l<N;l++){
                                w[3]=l;
                                if(p==4) { 
                                    struct node newJ;
                                    newJ.J = distribution(generator);
                                    newJ.SSS = 0.;
                                    for(int s=0; s<p; s++) { newJ.SSS += w[s]<<(16*i); }
                                    J.push_back(newJ);
                                    norm2 += newJ.J*newJ.J;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    RenormalizeNodes(J,norm2,N/2.);
}


void Planting(int p, double beta, vector<struct node> &JJ, vector<double> &S) {

    int N = S.size();
    int NJ = JJ.size();

    for(size_t i = 0; i < NJ; i++) {

        double plant_tilt = beta * N / NJ / 2.;
        long long int mask = 0b1111111111111111;
        
        for(int k=0; k<p; k++) { plant_tilt *= S[(JJ[i].SSS >> (16*k) & mask)]; }
    
        JJ[i].J -= plant_tilt;
    }
}

vector<double> White_Noise(int N, mt19937 &gauss, double std) {
 
    //cout << "NOT";
    normal_distribution<double> distribution(0,std);

    vector<double> X(N);

    for(int i=0;i<N;i++){    
        X[i]=distribution(gauss);
    }

    return X;
}

vector<double> White_discrete_Noise(int N, mt19937 &gauss, double std) {
 

    //cout << "YES";

    uniform_int_distribution<long long unsigned> distribution(0,0xFFFFFFFFFFFFFFFF);
    
    vector<double> X(N);

    long long unsigned randomNum;

    for(int i=0;i<N;i++){
        if(i%64==0) { 
            randomNum=distribution(gauss);
        }  
        X[i]= std*((randomNum >> i) & 1);
    }

    return X;
}


vector<double> RandomSpin(int seedS, int N) {

    double std = sqrt((double)N);

    vector<double> S(N, 0.0);

    mt19937 generator(seedS);
    normal_distribution<double> distribution(0,std);

    double norm,norm2=0.;
   
    int i;
    for(i=0;i<N;i++){        
        S[i] = distribution(generator);
        norm2 += S[i]*S[i];
    }

    RenormalizeVector(S,norm2,(double)N);

    return S;
}

double Energy(int p, vector<class node> &JJ, vector<double> &S) {
    double En=0.;

    int NJ = JJ.size();

    for(size_t i = 0; i < NJ; i++) {
        double J = JJ[i].J;
        long long int SSS = JJ[i].SSS;

        long long int mask = 0b1111111111111111;
        double multiS = 1.;

        for(int i=0;i<p;i++) { multiS *= S[(SSS>>(16*i) & mask)]; /*cout << i << ' ' << S[i] << ' ' << multiS << endl;*/ }
        En += J*multiS;
    }
    return En;
}

vector<double> Gradient(int p, const vector<struct node> &JJ, vector<double> &S) {
    
    int N = S.size();
    int NJ = JJ.size();

    vector<double> G(N, 0.0);

    for(size_t j = 0; j < NJ; j++) {
            double J = JJ[j].J;

            long long int mask = 0b1111111111111111;
            long long int SSS = JJ[j].SSS;

            int i,i2;
            for(i=0;i<p;i++){
                double g=J;
                for(i2=1;i2<p;i2++){
                    int ii = (i+i2)%p;
                    g *= S[(SSS>>(16*ii) & mask)];
                }
                G[(SSS>>(16*i) & mask)]+=g;
            }
    }
    return G;
}

vector<double> Hessian(int p, vector<struct node> &JJ, vector<double> &S) {

    int N = S.size();
    int NJ = JJ.size();

    vector<double> H(N*N, 0.0);

    for(size_t j = 0; j < NJ; j++) {
            double J = JJ[j].J;

            long long int mask = 0b1111111111111111;
            long long int SSS = JJ[j].SSS;

            for(int i=0;i<p;i++){
                int II = (SSS>>(16*i) & mask);
                for(int k=i+1;k<p;k++){
                    int K = (SSS>>(16*k) & mask);
                    double g=J;
                    for(int a=0;a<p;a++){
                        //if(a!=i && a!=k) { g*=S[(SSS>>(16*a) & mask)]; }
                        (a==i||a==k) ?  g*=1 : g*=S[(SSS>>(16*a) & mask)];
                    }
                    H[II*N + K]+=g;
                    H[K*N + II]+=g;
                }
            }
    }
    return H;
}

vector<double> Proj(vector<double> &X) {


    int N = X.size();

    vector<double> P(N*N,0.);

    double N2X = Norm2(X);
    int i,j,k;

    for(i=0;i<N;i++){
        P[N*i+i] += 1;
        for(j=0;j<N;j++){
            P[N*i+j] += -X[i]*X[j]/N2X;
        }
    }

    return P;
}

vector<double> Id(int N) {

    vector<double> P(N*N,0.);

    int i;

    for(i=0;i<N;i++){
        P[N*i+i] = 1;
    }

    return P;
}

vector<double> HProj(vector<double> &H, vector<double> &P, int N) {

    int i,j,k;

    vector<double> HP(N*N,0.);

    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            HP[N*i+j]=0.;
            for(k=0;k<N;k++){
                HP[N*i+j] += H[N*i+k]*P[N*k+j];
            }
        }
    }

    return HP;
    /*for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            H[N*i+j]=0.;
            for(k=0;k<N;k++){
                H[N*i+j] += P[N*i+k]*HP[N*k+j];
            }
        }
    }*/
}

vector<double> ProjVec(vector<double> &V, vector<double> &S) {
	int N = V.size();
    vector<double> PV(N,0.);
	double mu = -Prod(S,V);

 	return Lin(1.,V,mu/N,S);
}

vector<double> Multiply(vector<double> &H, vector<double> &V) {
	int N = V.size();
    vector<double> HV(N,0.);

    int i,k;

	for(i=0;i<N;i++){
        for(k=0;k<N;k++){
            
            HV[i] += H[N*i+k]*V[k];
        }
    }
    return HV;
}

double MinEigen(vector<double> &H, vector<double> &P, vector<double> &mE, double HighestEigenValue, double err) {

    //using namespace Eigen;

    int N = mE.size();
    
    /*MatrixXd HH(N,N);
    MatrixXd PP(N,N);
    VectorXd mmEE(N);
    VectorXd mmEE_prev(N);*/


    //vector<double> PmE(N,0.);
    vector<double> mE_prev(N,0.);
    //RenormalizeVector(mE,Norm2(mE),1);

    for(int i = 0; i < N; i++) {
        H[N*i+i] -= HighestEigenValue;
    }



    /*for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            HH(i,j) = H[i*N+j];
            PP(i,j) = P[i*N+j];
        }
        mmEE(i) = mE[i];
    }*/

    double eigen,prev_eigen;

    int k=0;
    eigen=1;
    do {

        prev_eigen = eigen;

        //mmEE_prev=mmEE;
        mE_prev=mE;
        //mmEE=PP*mmEE;
        mE=ProdMat(P,mE);
        //mmEE=HH*mmEE;
        mE=ProdMat(H,mE);
        //eigen=mmEE_prev(0)/mmEE(0);
        eigen = (mE_prev[0]/mE[0]+mE_prev[1]/mE[1])/2.;
        //mmEE=mmEE*eigen;
        mE=Scal(mE,eigen);   
        //RenormalizeVector(mE,Norm2(mE),1);
        k++;

        //cout << "ppp " << Norm2(mE) << ' ' << k << endl;

    } while ((eigen-prev_eigen)*(eigen-prev_eigen) > err);// (k<10); //Norm2(Sub(mE,mE_prev))/N/N>err);

    //cout << "ppp " << Norm2(mE) << ' ' << k << endl;

    //mmEE_prev=mmEE;
    mE_prev=mE;
    //mmEE=PP*mmEE;
    mE=ProdMat(P,mE);
    //mmEE=HH*mmEE;
    mE=ProdMat(H,mE);
    //mmEE=PP*mmEE;
    mE=ProdMat(P,mE);

    //double ee = mmEE(0)/mmEE_prev(0)+HighestEigenValue;
    double ee = mE[0]/mE_prev[0]+HighestEigenValue;

    //mmEE=mmEE/sqrt(mmEE.dot(mmEE));
    RenormalizeVector(mE,Norm2(mE),1);

    for(int i = 0; i < N; i++) {
        H[N*i+i] += HighestEigenValue;
    }

    //for(int i = 0; i < N; i++) {
    //    mE[i] = mmEE(i);
    //}

    //printf("e %d %1.12f",k,mE[1]/mE_prev[1]);

    //printf("H %1.10f %1.10f %1.10f %1.10f\n",H[0],H[1],H[2],sqrt(Norm2(mE))-HighestEigenValue);

    return ee;
}


vector<double> SoftLangevinSTEP(vector<double> S1, vector<double> G, double alpha) {


    double N = S1.size();
    vector<double> S2(N);

    //vector<double> sG2 = SoftenGradient(S1);

    for (int i=0; i<N; i++) {
        S2[i] = S1[i] - alpha*G[i];
    }

    RenormalizeVector(S2,Norm2(S2),(double)N);

    return S2;
}

vector<double> HardLangevinSTEP(vector<double> S1, vector<double> PG, double angle) {

    vector<double> S2;

    double N = S1.size();

    S2 = RotateVector(S1,PG,angle);

    return S2;    
}

//(const double *Gp, const double *RCG, double *CG, const double beta, const double n)
vector<double> NewConjugateGradient(vector<double> PG1, vector<double> RPG0, vector<double> RCG0) {

    double BetaPR=0.;
    int N = PG1.size();
    
    for(int i=0; i<N; i++) { BetaPR+=PG1[i]*(PG1[i]-RPG0[i]); }
    BetaPR/=Norm2(RPG0);

    vector<double> CG1(N);

    for(int i=0; i<N; i++) {
        CG1[i] = PG1[i]+BetaPR*RCG0[i];
    }

    return CG1;
}
