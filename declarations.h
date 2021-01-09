
struct parameters;
struct variables;

//dynamics.h
void Initialize_GD(variables *v, parameters *p);
vector<double> Langevin_step(vector<double> S1, vector<double> PG1, variables *v, parameters *p);
void Initialize_CG(variables *v, parameters *p);
vector<double> CG_Step(vector<double> S1, vector<double> PG1, variables *v, parameters *p, double (*TotEN)(vector<double> S1, variables *v, parameters *p));
double Select_CG_angle(vector<double> S, vector<double> CG, double alpha_max, variables *v, parameters *p);

//initialization.h
vector<double> Initialize_System(variables *v, parameters *p);

//inout.h
void Create_Directory(char *dir);
void SnapVector(vector<double> V, char *fname);
void SaveVector(variables *v, parameters *p, vector<double> V);
vector<double> OpenVector(int N,char *fname);
void SaveSystem(variables *v, parameters *p);
void OpenSystem(variables *v, parameters *p, char *fname);

//observables.h
void Evaluate_Observables(variables *v, parameters *p);
double Total_Energy(vector<double> S1, variables *v, parameters *p);
double Total_Energy_W(vector<double> S1, variables *v, parameters *p);
vector<double> Total_Gradient(vector<double> S1, variables *v, parameters *p);
vector<double> Total_Projected_Gradient(vector<double> S1, variables *v, parameters *p);
vector<double> Total_Projected_Gradient_W(vector<double> S1, variables *v, parameters *p);
vector<double> Projected_Noise(vector<double> S1, variables *v, parameters *p, double Temp);
vector<double> Magnetization(vector<double> S1, vector<double> M, variables *v, parameters *p);

//parametr.h
void Input_Parameters(int Ac, char **Av, struct variables *v, struct parameters *p);

//print.h
void Print_Initial_Observables(variables *v, parameters *p);
void Print_Observables(variables *v, parameters *p);
void Print_Eigen(variables *v, parameters *p);

//
vector<double> Find_Inherent_Structure(vector<double> S1,variables *v, parameters *p);

