void auto_display_format(int panel,int control,double min,double max);

void	Refresh_rho(void);
void	Refresh_sigma(void);
void	Refresh_delta(void);
void	Draw_bunch(void);
void	Set_zero_color(void);
void	Set_zero_color_one_control(void);
void	Set_y_z(void);
void	RecallPS (int P, char* fn);
void	SavePS (int P, char* fn);

int PH;

int NUM,yz;

double XMIN,XMAX;
double A_R,B_R,C_R,D_R,O_R,Fi_R,E_R,X0_R,S_R,F_R;
double Normirovka;
double A_S[2],B_S[2],C_S[2],D_S[2],O_S[2],Fi_S[2],E_S[2],X0_S[2],S_S[2],F_S[2];
double A_D[2],B_D[2],C_D[2],D_D[2],O_D[2],Fi_D[2],E_D[2],X0_D[2],S_D[2],F_D[2];
double RHO[10000],SIGMA[2][10000],DELTA[2][10000],X[10000];

double matrix_make_bunch[3][10][2];

int CT[256];
