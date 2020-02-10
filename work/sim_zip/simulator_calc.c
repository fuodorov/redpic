#include <analysis.h>
#include <utility.h>
#include <formatio.h>
#include <userint.h>
#include "toolbox.h"
#include <ansi_c.h>
#include "simulator_vars.h"
#include "simulator_draw.h"
#include "simulator_beam.h"
#include "make_bunch_my.h"
#include "lDSystem.h"
#include "math.h"

int plotHandle=0, N; double *K,*KX;
struct gauss_bunch bunch;

double gauss(double sigma)
{
	double x1,x2;
	do
	{
		x1=Random(0,1);
		x2=Random(0,1);
	}
	while (x1==0);
	return sigma*sqrt(-2*log(x1))*cos(2*pi*x2);
}

double gaussrand(double sigma)
{
	double U, V, Z;
	U = (rand() + 1.) / (RAND_MAX + 2.);
	V = rand() / (RAND_MAX + 1.);
	Z = sqrt(-2 * log(U)) * cos(2 * pi * V);
	return Z*sigma;
}

int sign (double x)
{
	if(x>=0) return  1;
	else 	 return -1;
}

void loadbeam(struct particle* b1,struct parametrs p)
{
	int i;
	switch (p.t_distr)
	{
		case 0: //распределение по сетке квадратное
			for (i=0;i<p.t_number;i++)
			{
				b1->x=(   ( ceil((i+1)/sqrt(p.t_number))-1) /(sqrt(p.t_number)-1)   )*2*p.t_radius - p.t_radius;
				b1->y=(   (((i+1)%(int)sqrt(p.t_number))-0) /(sqrt(p.t_number)-1)   )*2*p.t_radius - p.t_radius;
				b1->z=(p.t_length*((double)(i-p.t_number)/p.t_number));
				b1->px=Random(-sqrt(p.t_temp*((p.t_temp/me)+2)/me),sqrt(p.t_temp*((p.t_temp/me)+2)/me));
				b1->py=Random(-sqrt(p.t_temp*((p.t_temp/me)+2)/me),sqrt(p.t_temp*((p.t_temp/me)+2)/me));
				b1->pz=sqrt(p.t_energy*1e3*((p.t_energy*1e3/me)+2)/me);
				b1->status=in_cathode;
				b1++;
			}
		break;
		case 1: //равномерное распределение квадратное
			for (i=0;i<p.t_number;i++)
			{
				b1->x=Random(-p.t_radius,p.t_radius);
				b1->y=Random(-p.t_radius,p.t_radius);
				b1->z=Random(-p.t_length,0);
				b1->px=Random(-sqrt(p.t_temp*((p.t_temp/me)+2)/me),sqrt(p.t_temp*((p.t_temp/me)+2)/me));
				b1->py=Random(-sqrt(p.t_temp*((p.t_temp/me)+2)/me),sqrt(p.t_temp*((p.t_temp/me)+2)/me));
				b1->pz=sqrt(p.t_energy*1e3*((p.t_energy*1e3/me)+2)/me);
				b1->status=in_cathode;
				b1++;
			}
		break;
		case 2: //равномерное распределение круглое
			for (i=0;i<p.t_number;i++)
			{
				do
				{
					b1->x=Random(-p.t_radius,p.t_radius);
					b1->y=Random(-p.t_radius,p.t_radius);
				}
				while(pow(b1->x,2)+pow(b1->y,2)>pow(p.t_radius,2));
				b1->z=Random(-p.t_length,0);
				do
				{
					b1->px=Random(-sqrt(p.t_temp*((p.t_temp/me)+2)/me),sqrt(p.t_temp*((p.t_temp/me)+2)/me));
					b1->py=Random(-sqrt(p.t_temp*((p.t_temp/me)+2)/me),sqrt(p.t_temp*((p.t_temp/me)+2)/me));
				}
				while(pow(b1->px,2)+pow(b1->py,2)>pow(sqrt(p.t_temp*((p.t_temp/me)+2)/me),2));
				b1->pz=sqrt(p.t_energy*1e3*((p.t_energy*1e3/me)+2)/me)*(1+gauss((p.t_spread/p.t_energy)*(p.t_energy*1e3+me)/(p.t_energy+2*me)));
				//b1->pz=sqrt(p.t_energy*1e3*((p.t_energy*1e3/me)+2)/me);
				b1->status=in_cathode;
				b1++;
			}
		break;
		case 3: //нормальное распределение squared
			for (i=0;i<p.t_number;i++)
			{
				//b1->x=Random(-p.t_radius,p.t_radius);
				//b1->y=Random(-p.t_radius,p.t_radius);
				b1->x=gauss(p.t_radius);
				b1->y=gauss(p.t_radius);
				b1->z=Random(-p.t_length,0);
				b1->px=gauss(sqrt(p.t_temp*((p.t_temp/me)+2)/me));
				b1->py=gauss(sqrt(p.t_temp*((p.t_temp/me)+2)/me));
				b1->pz=sqrt(p.t_energy*1e3*((p.t_energy*1e3/me)+2)/me);
				b1->status=in_cathode;
				b1++;
			}
	}
}

struct gauss_bunch loadbunch(struct parametrs p)
{
	FILE *f;
	int i; double integral=0,m1=0,summa=0,bb;
	struct gauss_bunch b={0};
	f = fopen ("gauss.txt", "r");
	if (!f) { MessagePopup ("Error", "File 'gauss.txt' not found"); return b; }
	fscanf (f, "%d", &b.n);
	b.t = malloc(b.n*sizeof(double));
	b.i = malloc(b.n*sizeof(double));
	b.s = malloc(b.n*sizeof(double));
	b.y = malloc(b.n*sizeof(double));
	b.z = malloc(b.n*sizeof(double));
	for (i=0;i<b.n;i++) fscanf (f, "%lf%lf%lf%lf%lf%lf", &b.t[i],&b.i[i],&b.s[i],&b.y[i],&bb,&b.z[i]);
	fclose (f);
	for (i=0;i<b.n-1;i++) integral+=(b.i[i+1]+b.i[i])*(b.t[i+1]-b.t[i])/2;
	for (i=0;i<b.n;i++) { m1+=b.t[i]*b.i[i]; summa+=b.i[i]; }
	m1/=summa;
	for (i=0;i<b.n;i++)
	{
		b.t[i]=(b.t[i]-m1)/p.step;
		b.i[i]/=p.b_v*integral/p.step;
		b.s[i]*=fl/p.step;
		b.y[i]*=fl/p.step;
		b.z[i]*=fl/p.step;
	}
	return b;
}

struct vector E_gun(struct particle b, struct parametrs p)
{
	struct vector E={0};
	int i; double r;
	r=sqrt(b.x*b.x+b.y*b.y);
	if(  r>(0.03*fl/p.step)  )  r = 0.03*fl/p.step;
	i=99*b.z/p.gun_length;
	if(i==99) i=98;
	E.z=(          Ez_gun[i][0]+(Ez_gun[i+1][0]-Ez_gun[i][0])*(99*b.z/p.gun_length-i)
		      +(  (Ez_gun[i][1]+(Ez_gun[i+1][1]-Ez_gun[i][1])*(99*b.z/p.gun_length-i))-(Ez_gun[i][0]+(Ez_gun[i+1][0]-Ez_gun[i][0])*(99*b.z/p.gun_length-i))   )*r/(0.03*fl/p.step)                  )*100000*fe*p.step/30000;
	E.x=(          Er_gun[i][0]+(Er_gun[i+1][0]-Er_gun[i][0])*(99*b.z/p.gun_length-i)
		      +(  (Er_gun[i][1]+(Er_gun[i+1][1]-Er_gun[i][1])*(99*b.z/p.gun_length-i))-(Er_gun[i][0]+(Er_gun[i+1][0]-Er_gun[i][0])*(99*b.z/p.gun_length-i))   )*r/(0.03*fl/p.step)                  )*100000*fe*p.step/30000
		*b.x/r;
	E.y=(          Er_gun[i][0]+(Er_gun[i+1][0]-Er_gun[i][0])*(99*b.z/p.gun_length-i)
		      +(  (Er_gun[i][1]+(Er_gun[i+1][1]-Er_gun[i][1])*(99*b.z/p.gun_length-i))-(Er_gun[i][0]+(Er_gun[i+1][0]-Er_gun[i][0])*(99*b.z/p.gun_length-i))   )*r/(0.03*fl/p.step)                  )*100000*fe*p.step/30000
		*b.y/r;
	return E;
}

struct vector E_obj(int t, struct particle b, struct parametrs p, struct object ob)
{
	struct vector E = {0}, E1={0}, E2 = {0};
	int i,ind = 0, t_ob = (ob.parametrs[0]+p.t_length/2)/p.t_v; //time of flight until object center
	double bx_1=0, by_1=0, r_1=0, r1, sigma_r_1=0, Er_1=0, Er1=0, nx_1=0, at, Xt, bz=0;
	double bx_2=0, by_2=0, r_2=0, r2, sigma_r_2=0, Er_2=0, Er2=0, nx_2=0, dx;
	switch (ob.name)
	{
		case in_scan:
			//if(((t-t_ob)>=-ob.parametrs[5]/2)&&((t-t_ob)<+ob.parametrs[5]/2))
			if( fabs(t-t_ob)<=ob.parametrs[5]/2 )
			{
				switch ((int)ob.parametrs[6])
				{
					case 0: // linear scan field
						E.x=(ob.parametrs[3]/ob.parametrs[2])*2*(t-t_ob)/ob.parametrs[5];
						E.y=(ob.parametrs[4]/ob.parametrs[2])*2*(t-t_ob)/ob.parametrs[5];
						E.z=0;
					break;
					case 1: // scan field from file
						DSfield(-ob.parametrs[3]*((t-t_ob)/ob.parametrs[5])*300/(fe*fl), b.x*p.step/fl, b.y*p.step/fl, b.z*p.step/fl, &E.y, &E.x, &E.z);
						E.x*=fe*p.step/30000;
						E.y*=fe*p.step/30000;
						E.z*=fe*p.step/30000;
					break;
					case 2:
						E.x=(ob.parametrs[3]/ob.parametrs[2])*(t-t_ob)*pow(     1+2*pow((t-t_ob)/ob.parametrs[5],1)     ,1)/ob.parametrs[5];
						E.y=(ob.parametrs[4]/ob.parametrs[2])*(t-t_ob)/ob.parametrs[5];
						E.z=0;
				}
			}
			return E;
		case in_bunch:
			switch (p.cb) // 1)формируем продольную зависимость в банче 
			{
				case 0: //банч из простых уставок
					by_1 = ob.Y_1; bx_1 = ob.X_1;
					by_2 = ob.Y_2; bx_2 = ob.X_2;
					bz = ob.parametrs[0];
					r_1 = sqrt(  pow( b.y - by_1,2)+pow( b.z - bz - bx_1,2));
					r_2 = sqrt(  pow( b.y - by_2,2)+pow( b.z - bz - bx_2,2));
					sigma_r_1=ob.R_1;
					sigma_r_2=ob.R_2;
					switch ((int)ob.parametrs[8])  //выбираем продольное распределение
					{
						case 0://gauss
					//		if (fabs(b.x-p.b_v*(t-t_ob))<5*ob.parametrs[2])
								nx_1 = exp(    -pow(b.x-ob.Z_1-p.b_v*(t-t_ob),2)/(2*pow(ob.L_1,2))   )/(ob.L_1*sqrt(2*pi));
								nx_2 = exp(    -pow(b.x-ob.Z_2-p.b_v*(t-t_ob),2)/(2*pow(ob.L_2,2))   )/(ob.L_2*sqrt(2*pi));
					//		else
					//			nx_1 = 0; nx_2 = 0;
						break;
						case 1://uniform
							if (fabs(b.x-ob.Z_1-p.b_v*(t-t_ob))> (ob.L_1/2)) nx_1 = 0; 			else nx_1 = 1/ob.L_1;
							if (fabs(b.x-ob.Z_2-p.b_v*(t-t_ob))<=(ob.L_2/2)) nx_2 = 1/ob.L_2;   else nx_2 = 0;
						break;
						case 2://triangle
							if (fabs(b.x-ob.Z_1-p.b_v*(t-t_ob))>(ob.L_1/2)) nx_1 = 0; else nx_1 = (1-2*fabs(b.x-ob.Z_1-p.b_v*(t-t_ob))/ob.L_1)/(ob.L_1/2);
							if (fabs(b.x-ob.Z_2-p.b_v*(t-t_ob))>(ob.L_2/2)) nx_2 = 0; else nx_2 = (1-2*fabs(b.x-ob.Z_2-p.b_v*(t-t_ob))/ob.L_2)/(ob.L_2/2);
						break;
						case 3://user non-simmetric 2 different Gauss
						//	if (b.x-p.b_v*(t-t_ob)>=0) nx_1=exp(    -pow(b.x-p.b_v*(t-t_ob),2)/(2*pow(ob.parametrs[2]/4,2))   )/(ob.parametrs[2]*sqrt(2*pi));
						//	else                       nx_1=exp(    -pow(b.x-p.b_v*(t-t_ob),2)/(2*pow(ob.parametrs[2],2))   )/((ob.parametrs[2]/1)*sqrt(2*pi));
							if (b.x-ob.Z_1-p.b_v*(t-t_ob)>=0) nx_1 = exp(    -pow(b.x-ob.Z_1-p.b_v*(t-t_ob),2)/(2*pow(ob.L_1,2))   )/(ob.L_1*sqrt(2*pi));
							else                       		  nx_1 = exp(    -pow(b.x-ob.Z_1-p.b_v*(t-t_ob),2)/(2*pow(ob.L_2,2))   )/(ob.L_1*sqrt(2*pi));
					}
				break;
				case 1: //банч из файла
					if (  (  ((t-t_ob)*p.b_v-b.x)<bunch.t[0]*p.b_v  ) || (  ((t-t_ob)*p.b_v-b.x)>bunch.t[bunch.n-1]*p.b_v  )  ) return E;
					for (i=0;i<bunch.n;i++)
					{
						if (  bunch.t[i]*p.b_v>((t-t_ob)*p.b_v-b.x)  ) { ind=i-1; break; }
					}
					if ((ind==-1)||(ind==bunch.n-2)) ind=0;
					at=(  (t-t_ob-(b.x/p.b_v)) - bunch.t[ind]  )/(  bunch.t[ind+1] - bunch.t[ind]  );
					by_1 = ob.parametrs[4] + bunch.y[ind] + (  bunch.y[ind+1] - bunch.y[ind]  )*at;
					bz = ob.parametrs[0]   + bunch.z[ind] + (  bunch.z[ind+1] - bunch.z[ind]  )*at;
					r_1 = sqrt(  pow( b.y - by_1 ,2) + pow( b.z - bz ,2)  );
					sigma_r_1 = bunch.s[ind] + (  bunch.s[ind+1] - bunch.s[ind]  )*at;
					nx_1      = bunch.i[ind] + (  bunch.i[ind+1] - bunch.i[ind]  )*at;
				break;
				case 2: //банч из сложных уставок
					Xt=(t-t_ob)*p.step;
					by_1 = (matrix_make_bunch[2][0][0]*Xt*Xt+matrix_make_bunch[2][1][0]*Xt+matrix_make_bunch[2][2][0])*
						 (matrix_make_bunch[2][3][0]*sin(matrix_make_bunch[2][4][0]*Xt+matrix_make_bunch[2][5][0])+matrix_make_bunch[2][6][0]*exp(-(Xt-matrix_make_bunch[2][7][0])*(Xt-matrix_make_bunch[2][7][0])/(2*matrix_make_bunch[2][8][0]*matrix_make_bunch[2][8][0]))+matrix_make_bunch[2][9][0])*
						 fl/p.step;
					bz = ob.parametrs[0] + 
					     (matrix_make_bunch[2][0][1]*Xt*Xt+matrix_make_bunch[2][1][1]*Xt+matrix_make_bunch[2][2][1])*
					     (matrix_make_bunch[2][3][1]*sin(matrix_make_bunch[2][4][1]*Xt+matrix_make_bunch[2][5][1])+matrix_make_bunch[2][6][1]*exp(-(Xt-matrix_make_bunch[2][7][1])*(Xt-matrix_make_bunch[2][7][1])/(2*matrix_make_bunch[2][8][1]*matrix_make_bunch[2][8][1]))+matrix_make_bunch[2][9][1])*
					     fl/p.step;
					r_1 = sqrt(  pow( b.y - by_1,2) + pow( b.z - bz ,2));
					yz=0;
					sigma_r_1 = (matrix_make_bunch[1][0][yz]*Xt*Xt+matrix_make_bunch[1][1][yz]*Xt+matrix_make_bunch[1][2][yz])*
					          (matrix_make_bunch[1][3][yz]*sin(matrix_make_bunch[1][4][yz]*Xt+matrix_make_bunch[1][5][yz])+matrix_make_bunch[1][6][yz]*exp(-(Xt-matrix_make_bunch[1][7][yz])*(Xt-matrix_make_bunch[1][7][yz])/(2*matrix_make_bunch[1][8][yz]*matrix_make_bunch[1][8][yz]))+matrix_make_bunch[1][9][yz])*
					          fl/p.step;
					nx_1 = (Normirovka*p.step)*(matrix_make_bunch[0][0][0]*Xt*Xt+matrix_make_bunch[0][1][0]*Xt+matrix_make_bunch[0][2][0])*
					     (matrix_make_bunch[0][3][0]*sin(matrix_make_bunch[0][4][0]*Xt+matrix_make_bunch[0][5][0])+matrix_make_bunch[0][6][0]*exp(-(Xt-matrix_make_bunch[0][7][0])*(Xt-matrix_make_bunch[0][7][0])/(2*matrix_make_bunch[0][8][0]*matrix_make_bunch[0][8][0]))+matrix_make_bunch[0][9][0]);
			}
// ========== bunch 1 считаем радиальное поле от первого банча ===============
			switch ((int)ob.TRD_1) //выбираем поперечное распределение
			{
				case 0://gauss
					Er_1 =     ob.K_1*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(1 - exp(-r_1*r_1/(2*pow(sigma_r_1,2)))  )*nx_1/r_1;
				break;
				case 1://uniform
					if (r_1>sigma_r_1) Er_1 = ob.K_1*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(  1  							  )*nx_1/r_1;
					else 			   Er_1 = ob.K_1*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(  r_1*r_1/(sigma_r_1*sigma_r_1)  )*nx_1/r_1;
				break;
				case 2://triangle
					if (r_1>sigma_r_1) Er_1 = ob.K_1*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(  1  )*nx_1/r_1;
					else               Er_1 = ob.K_1*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(  3*r_1*r_1/(sigma_r_1*sigma_r_1)-2*r_1*r_1*r_1/(pow(sigma_r_1,3))  )*nx_1/r_1;
				break;
				case 3://user
					if ( r_1> sigma_r_1) 					 Er_1 = ob.K_1*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(  1  )*nx_1/r_1;
					if ((r_1<=sigma_r_1)&&(r_1>sigma_r_1/2)) Er_1 = ob.K_1*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(  (sigma_r_1*sigma_r_1+4*r_1*r_1)/(5*sigma_r_1*sigma_r_1)  )*nx_1/r_1;
					if (r_1<=sigma_r_1/2) 					 Er_1 = ob.K_1*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(  8*r_1*r_1/(5*sigma_r_1*sigma_r_1)  )*nx_1/r_1;
				break;
				case 4: //flat horizontal
					for (i=1;i<=p.N_Flat;i++)
					{
						dx = (2*i-p.N_Flat-1)*2*sigma_r_1/2;
						r1 = sqrt(  pow( b.y - by_1,2)+pow( b.z - bz - bx_1 - dx,2));
						Er1 = ob.K_1*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(1 - exp(-r1*r1/(2*pow(sigma_r_1,2)))  )*nx_1/r1;
						Er1*= exp( -pow((2*i-p.N_Flat-1)*sigma_r_1,2)/(2*pow(sigma_r_1*10,2))); // gauss flat
						E1.y += Er1 * (b.y-by_1)/r1;
						E1.z += Er1 * (b.z-bz-bx_1-dx)/r1;
					}
					E.y += E1.y/p.N_Flat;
					E.z += E1.z/p.N_Flat;
				break;
				case 5: //flat vertical
					for (i=1;i<=p.N_Flat;i++)
					{
						dx = (2*i-p.N_Flat-1)*2*sigma_r_1/2;
						r1 = sqrt(  pow( b.y - by_1 - dx,2)+pow( b.z - bz - bx_1,2));
						Er1 = ob.K_1*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(1 - exp(-r1*r1/(2*pow(sigma_r_1,2)))  )*nx_1/r1;
						Er1*= exp( -pow((2*i-p.N_Flat-1)*sigma_r_1,2)/(2*pow(sigma_r_1*10,2))); // gauss flat
//ob.K_1*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(1 - exp(-r_1*r_1/(2*pow(sigma_r_1,2)))  )*nx_1/r_1;
						E1.y += Er1 * ( b.y - by_1 - dx )/r1;
						E1.z += Er1 * ( b.z - bz   - bx_1 )/r1;
					}
					E.y += E1.y/p.N_Flat;
					E.z += E1.z/p.N_Flat;
				break;
			}
// ============= bunch 2 считаем радиальное поле от второго банча ==============
			switch ((int)ob.TRD_2)
			{
				case 0://gauss
					Er_2 = ob.K_2*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(1 - exp(    -r_2*r_2/(2*pow(sigma_r_2,2))    ))*nx_2/r_2;
				break;
				case 1://uniform
					if (r_2>sigma_r_2) Er_2 = ob.K_2*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(  1  )*nx_2/r_2;
					else               Er_2 = ob.K_2*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(  r_2*r_2/(sigma_r_2*sigma_r_2)  )*nx_2/r_2;
				break;
				case 2://triangle
					if (r_2>sigma_r_2) Er_2 = ob.K_2*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(  1  )*nx_2/r_2;
					else               Er_2 = ob.K_2*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(  3*r_2*r_2/(sigma_r_2*sigma_r_2)-2*pow(r_2,3)/(pow(sigma_r_2,3))  )*nx_2/r_2;
				break;
				case 3://user
					if (r_2>sigma_r_2) 						 Er_2 = ob.K_2*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(  1  )*nx_2/r_2;
					if ((r_2<=sigma_r_2)&&(r_2>sigma_r_2/2)) Er_2 = ob.K_2*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(  (sigma_r_2*sigma_r_2+4*r_2*r_2)/(5*sigma_r_2*sigma_r_2)  )*nx_2/r_2;
					if (r_2<=sigma_r_2/2)					 Er_2 = ob.K_2*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(  8*r_2*r_2/(5*sigma_r_2*sigma_r_2)  )*nx_2/r_2;
				break;
				case 4: //flat horizontal
					for (i=1;i<=p.N_Flat;i++)
					{
						dx = (2*i-p.N_Flat-1)*2*sigma_r_2/2;
						r2 = sqrt(  pow( b.y - by_2,2)+pow( b.z - bz - bx_2 - dx,2));
						Er2 = ob.K_2*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(1 - exp(-r2*r2/(2*pow(sigma_r_2,2)))  )*nx_2/r2;
						Er2*= exp( -pow((2*i-p.N_Flat-1)*sigma_r_2,2)/(2*pow(sigma_r_2*10,2))); // gauss flat
						E2.y += Er2 * (b.y-by_2)/r2;
						E2.z += Er2 * (b.z-bz-bx_2-dx)/r2;
					}
					E.y += E2.y/p.N_Flat;
					E.z += E2.z/p.N_Flat;
				break;
				case 5: //flat vertical
					for (i=1;i<=p.N_Flat;i++)
					{
						dx = (2*i-p.N_Flat-1)*2*sigma_r_2/2;
						r2 = sqrt(  pow( b.y - by_2 - dx,2)+pow( b.z - bz - bx_2,2));
						Er2 = ob.K_2*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(1 - exp(-r2*r2/(2*pow(sigma_r_2,2)))  )*nx_2/r2;
						Er2*= exp( -pow((2*i-p.N_Flat-1)*sigma_r_2,2)/(2*pow(sigma_r_2*10,2))); // gauss flat
						E2.y += Er2 * ( b.y - by_2 - dx )/r2;
						E2.z += Er2 * ( b.z - bz   - bx_2 )/r2;
					}
					E.y += E2.y/p.N_Flat;
					E.z += E2.z/p.N_Flat;
				break;
			}
			E.y += ((b.y-by_1   )/r_1)*Er_1+((b.y-by_2   )/r_2)*Er_2;
			E.z += ((b.z-bz-bx_1)/r_1)*Er_1+((b.z-bz-bx_2)/r_2)*Er_2;
			//E.z = ((z-bz)  /r_1)*Er_1+((z-bz)/  r_2)*Er_2; % original from 2009
			return E;
		case in_RF:
			if (ob.parametrs[5]==0)   //ideal
			{
				E.x=(ob.parametrs[2]/ob.parametrs[1])*cos(2*pi*(t-t_ob)/ob.parametrs[3]);//*cos(2*pi*(z-ob.parametrs[0])/ob.parametrs[4]);
				return E;
			}
			//if (ob.parametrs[5]==1)   //real
			//if (ob.parametrs[5]==2)   //file
		break;
		case in_acc:
			E.z = -(ob.parametrs[2]/ob.parametrs[1])*cos(2*pi*(t-t_ob+(b.z-ob.parametrs[0]))/ob.parametrs[3]);//*cos(2*pi*(z-ob.parametrs[0])/ob.parametrs[4]);
			return E;
	}
	return E;
}

struct vector H_obj(int t, struct particle b, struct parametrs p, struct object ob)
{
	struct vector H={0},H1={0}, H2={0};
	double alpha,betta,delta,d;
	//double An[10]={3.386,3185,-0.5636,0.1438,-0.01614,-0.786,-0.06714,0.01292,-0.06255,-0.21975};
	//double Bn[10]={1.790,30.09,-0.5121,-0.1921,-0.09311,-0.06929,-0.08558,0.1308,-0.06829,-0.09957};
	double An[10]={0.001,1,0.01,0.01,0.01,0.01,0.1,0.01,0.01,0.01};
	double Bn[10]={0.001,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
	double Hr_1=0,Hr_2=0,Hfi=0, fi, r0, A2;
	int i,j,k,k1,k2,ind=0;
	int t_ob=(ob.parametrs[0]+p.t_length/2)/p.t_v;	 //time of flight until object 
	double bx_1=0, by_1=0, r_1=0, r1, sigma_r_1=0, nx_1=0, at, Xt, bz=0;
	double bx_2=0, by_2=0, r_2=0, r2, sigma_r_2=0, nx_2=0, dx, Hr=0;
	double x1,y1,z1,z2,phi;
	double a1,b1,c1,d1,e1,f1,g1,h1;
	switch (ob.name)
	{
		case in_quad: //двойной квадруполь
			phi=ob.parametrs[6]*pi/180; //поворот системы координат на 45 градусов
			//phi=45*pi/180;
			x1 = (b.x - ob.parametrs[7])*cos(phi)-(b.y - ob.parametrs[8])*sin(phi);
			y1 = (b.x - ob.parametrs[7])*sin(phi)+(b.y - ob.parametrs[8])*cos(phi);
			if ((fabs(x1)>quadrupole_parametrs.size_x)||(fabs(y1)>quadrupole_parametrs.size_y))
			{
				H.x = 0;
				H.y = 0;
				H.z = 0;
				return H;
			}
			if (b.z<(ob.parametrs[3]-10*fl/p.step))
			{
				z1 = b.z-ob.parametrs[2];
				if (fabs(z1)>quadrupole_parametrs.size_z)
				{
					H.x = 0;
					H.y = 0;
					H.z = 0;
					return H;
				}
				i = fabs(x1)/(quadrupole_parametrs.size_dx);
				j = fabs(y1)/(quadrupole_parametrs.size_dy);
				k = fabs(z1)/(quadrupole_parametrs.size_dz);
				b.x = (fabs(x1)/quadrupole_parametrs.size_dx)-i;
				b.y = (fabs(y1)/quadrupole_parametrs.size_dy)-j;
				b.z = (fabs(z1)/quadrupole_parametrs.size_dz)-k;
				/*
				H.x = 1*quadrupole_fields1[i][j][k].x*sign(y1);
				H.y = 1*quadrupole_fields1[i][j][k].y*sign(x1);
				H.z = 1*quadrupole_fields1[i][j][k].z*sign(x1*y1*z);
				*/
				a1=quadrupole_fields1[i  ][j  ][k  ].x;
				b1=quadrupole_fields1[i+1][j  ][k  ].x-a1;
				c1=quadrupole_fields1[i  ][j+1][k  ].x-a1;
				d1=quadrupole_fields1[i  ][j  ][k+1].x-a1;
				e1=quadrupole_fields1[i+1][j+1][k  ].x-(a1+b1+c1);
				f1=quadrupole_fields1[i+1][j  ][k+1].x-(a1+b1+d1);
				g1=quadrupole_fields1[i  ][j+1][k+1].x-(a1+c1+d1);
				h1=quadrupole_fields1[i+1][j+1][k+1].x-(a1+b1+c1+d1+e1+f1+g1);
				H1.x=(a1+b1*b.x+c1*b.y+d1*b.z+e1*b.x*b.y+f1*b.x*b.z+g1*b.y*b.z+h1*b.x*b.y*b.z)*sign(y1)*ob.parametrs[4];
				a1=quadrupole_fields1[i  ][j  ][k  ].y;
				b1=quadrupole_fields1[i+1][j  ][k  ].y-a1;
				c1=quadrupole_fields1[i  ][j+1][k  ].y-a1;
				d1=quadrupole_fields1[i  ][j  ][k+1].y-a1;
				e1=quadrupole_fields1[i+1][j+1][k  ].y-(a1+b1+c1);
				f1=quadrupole_fields1[i+1][j  ][k+1].y-(a1+b1+d1);
				g1=quadrupole_fields1[i  ][j+1][k+1].y-(a1+c1+d1);
				h1=quadrupole_fields1[i+1][j+1][k+1].y-(a1+b1+c1+d1+e1+f1+g1);
				H1.y=(a1+b1*b.x+c1*b.y+d1*b.z+e1*b.x*b.y+f1*b.x*b.z+g1*b.y*b.z+h1*b.x*b.y*b.z)*sign(x1)*ob.parametrs[4];
				a1=quadrupole_fields1[i  ][j  ][k  ].z;
				b1=quadrupole_fields1[i+1][j  ][k  ].z-a1;
				c1=quadrupole_fields1[i  ][j+1][k  ].z-a1;
				d1=quadrupole_fields1[i  ][j  ][k+1].z-a1;
				e1=quadrupole_fields1[i+1][j+1][k  ].z-(a1+b1+c1);
				f1=quadrupole_fields1[i+1][j  ][k+1].z-(a1+b1+d1);
				g1=quadrupole_fields1[i  ][j+1][k+1].z-(a1+c1+d1);
				h1=quadrupole_fields1[i+1][j+1][k+1].z-(a1+b1+c1+d1+e1+f1+g1);
				H.z=(a1+b1*b.x+c1*b.y+d1*b.z+e1*b.x*b.y+f1*b.x*b.z+g1*b.y*b.z+h1*b.x*b.y*b.z)*sign(x1*y1*z1)*ob.parametrs[4];
				//поворот системы координат обратно на 45 градусов
				H.x = H1.x*cos(-phi)-H1.y*sin(-phi);
				H.y = H1.x*sin(-phi)+H1.y*cos(-phi);
				return H;
			}
			if ((b.z>=(ob.parametrs[3]-(10*fl/p.step)))&&(b.z<(ob.parametrs[2]+(10*fl/p.step))))
			{
				z1 = b.z-ob.parametrs[2];
				z2 = b.z-ob.parametrs[3];
				i = fabs(x1)/(quadrupole_parametrs.size_dx);
				j = fabs(y1)/(quadrupole_parametrs.size_dy);
				k1= fabs(z1)/(quadrupole_parametrs.size_dz);
				k2= fabs(z2)/(quadrupole_parametrs.size_dz);
				b.x = (fabs(x1)/quadrupole_parametrs.size_dx)-i;
				b.y = (fabs(y1)/quadrupole_parametrs.size_dy)-j;
				b.z = (fabs(z1)/quadrupole_parametrs.size_dz)-k1;
				a1=quadrupole_fields1[i  ][j  ][k1  ].x;
				b1=quadrupole_fields1[i+1][j  ][k1  ].x-a1;
				c1=quadrupole_fields1[i  ][j+1][k1  ].x-a1;
				d1=quadrupole_fields1[i  ][j  ][k1+1].x-a1;
				e1=quadrupole_fields1[i+1][j+1][k1  ].x-(a1+b1+c1);
				f1=quadrupole_fields1[i+1][j  ][k1+1].x-(a1+b1+d1);
				g1=quadrupole_fields1[i  ][j+1][k1+1].x-(a1+c1+d1);
				h1=quadrupole_fields1[i+1][j+1][k1+1].x-(a1+b1+c1+d1+e1+f1+g1);
				H1.x=(a1+b1*b.x+c1*b.y+d1*b.z+e1*b.x*b.y+f1*b.x*b.z+g1*b.y*b.z+h1*b.x*b.y*b.z)*sign(y1)*ob.parametrs[4];
				a1=quadrupole_fields1[i  ][j  ][k1  ].y;
				b1=quadrupole_fields1[i+1][j  ][k1  ].y-a1;
				c1=quadrupole_fields1[i  ][j+1][k1  ].y-a1;
				d1=quadrupole_fields1[i  ][j  ][k1+1].y-a1;
				e1=quadrupole_fields1[i+1][j+1][k1  ].y-(a1+b1+c1);
				f1=quadrupole_fields1[i+1][j  ][k1+1].y-(a1+b1+d1);
				g1=quadrupole_fields1[i  ][j+1][k1+1].y-(a1+c1+d1);
				h1=quadrupole_fields1[i+1][j+1][k1+1].y-(a1+b1+c1+d1+e1+f1+g1);
				H1.y=(a1+b1*b.x+c1*b.y+d1*b.z+e1*b.x*b.y+f1*b.x*b.z+g1*b.y*b.z+h1*b.x*b.y*b.z)*sign(x1)*ob.parametrs[4];
				a1=quadrupole_fields1[i  ][j  ][k1  ].z;
				b1=quadrupole_fields1[i+1][j  ][k1  ].z-a1;
				c1=quadrupole_fields1[i  ][j+1][k1  ].z-a1;
				d1=quadrupole_fields1[i  ][j  ][k1+1].z-a1;
				e1=quadrupole_fields1[i+1][j+1][k1  ].z-(a1+b1+c1);
				f1=quadrupole_fields1[i+1][j  ][k1+1].z-(a1+b1+d1);
				g1=quadrupole_fields1[i  ][j+1][k1+1].z-(a1+c1+d1);
				h1=quadrupole_fields1[i+1][j+1][k1+1].z-(a1+b1+c1+d1+e1+f1+g1);
				H.z=(a1+b1*b.x+c1*b.y+d1*b.z+e1*b.x*b.y+f1*b.x*b.z+g1*b.y*b.z+h1*b.x*b.y*b.z)*sign(x1*y1*z1)*ob.parametrs[4];
				/*
				H.x = ob.parametrs[4]*quadrupole_fields1[i][j][k1].x*sign(y1)      +   ob.parametrs[5]*quadrupole_fields2[i][j][k2].x*sign(y1);
				H.y = ob.parametrs[4]*quadrupole_fields1[i][j][k1].y*sign(x1)      +   ob.parametrs[5]*quadrupole_fields2[i][j][k2].y*sign(x1);
				H.z = ob.parametrs[4]*quadrupole_fields1[i][j][k1].z*sign(x1*y1*z1)+   ob.parametrs[5]*quadrupole_fields2[i][j][k2].z*sign(x1*y1*z2);
				*/
				b.z = (fabs(z2)/quadrupole_parametrs.size_dz)-k2;
				a1=quadrupole_fields2[i  ][j  ][k2  ].x;
				b1=quadrupole_fields2[i+1][j  ][k2  ].x-a1;
				c1=quadrupole_fields2[i  ][j+1][k2  ].x-a1;
				d1=quadrupole_fields2[i  ][j  ][k2+1].x-a1;
				e1=quadrupole_fields2[i+1][j+1][k2  ].x-(a1+b1+c1);
				f1=quadrupole_fields2[i+1][j  ][k2+1].x-(a1+b1+d1);
				g1=quadrupole_fields2[i  ][j+1][k2+1].x-(a1+c1+d1);
				h1=quadrupole_fields2[i+1][j+1][k2+1].x-(a1+b1+c1+d1+e1+f1+g1);
				H1.x+=(a1+b1*b.x+c1*b.y+d1*b.z+e1*b.x*b.y+f1*b.x*b.z+g1*b.y*b.z+h1*b.x*b.y*b.z)*sign(y1)*ob.parametrs[5];
				a1=quadrupole_fields2[i  ][j  ][k2  ].y;
				b1=quadrupole_fields2[i+1][j  ][k2  ].y-a1;
				c1=quadrupole_fields2[i  ][j+1][k2  ].y-a1;
				d1=quadrupole_fields2[i  ][j  ][k2+1].y-a1;
				e1=quadrupole_fields2[i+1][j+1][k2  ].y-(a1+b1+c1);
				f1=quadrupole_fields2[i+1][j  ][k2+1].y-(a1+b1+d1);
				g1=quadrupole_fields2[i  ][j+1][k2+1].y-(a1+c1+d1);
				h1=quadrupole_fields2[i+1][j+1][k2+1].y-(a1+b1+c1+d1+e1+f1+g1);
				H1.y+=(a1+b1*b.x+c1*b.y+d1*b.z+e1*b.x*b.y+f1*b.x*b.z+g1*b.y*b.z+h1*b.x*b.y*b.z)*sign(x1)*ob.parametrs[5];
				a1=quadrupole_fields2[i  ][j  ][k2  ].z;
				b1=quadrupole_fields2[i+1][j  ][k2  ].z-a1;
				c1=quadrupole_fields2[i  ][j+1][k2  ].z-a1;
				d1=quadrupole_fields2[i  ][j  ][k2+1].z-a1;
				e1=quadrupole_fields2[i+1][j+1][k2  ].z-(a1+b1+c1);
				f1=quadrupole_fields2[i+1][j  ][k2+1].z-(a1+b1+d1);
				g1=quadrupole_fields2[i  ][j+1][k2+1].z-(a1+c1+d1);
				h1=quadrupole_fields2[i+1][j+1][k2+1].z-(a1+b1+c1+d1+e1+f1+g1);
				H.z+=(a1+b1*b.x+c1*b.y+d1*b.z+e1*b.x*b.y+f1*b.x*b.z+g1*b.y*b.z+h1*b.x*b.y*b.z)*sign(x1*y1*z2)*ob.parametrs[5];
				//поворот системы координат обратно на 45 градусов
				H.x = H1.x*cos(-phi)-H1.y*sin(-phi);
				H.y = H1.x*sin(-phi)+H1.y*cos(-phi);
				return H;
			}
			if (b.z>=(ob.parametrs[2]+10*fl/p.step))
			{
				z2 = b.z-ob.parametrs[3];
				if (fabs(z2)>quadrupole_parametrs.size_z)
				{
					H.x = 0;
					H.y = 0;
					H.z = 0;
					return H;
				}
				i = fabs(x1)/(quadrupole_parametrs.size_dx);
				j = fabs(y1)/(quadrupole_parametrs.size_dy);
				k = fabs(z2)/(quadrupole_parametrs.size_dz);
				b.x = (fabs(x1)/quadrupole_parametrs.size_dx)-i;
				b.y = (fabs(y1)/quadrupole_parametrs.size_dy)-j;
				b.z = (fabs(z2)/quadrupole_parametrs.size_dz)-k;
				/*
				H.x = 0*quadrupole_fields2[i][j][k].x*sign(y1);
				H.y = 0*quadrupole_fields2[i][j][k].y*sign(x1);
				H.z = 0*quadrupole_fields2[i][j][k].z*sign(x1*y1*z2);
				*/
				a1=quadrupole_fields2[i  ][j  ][k  ].x;
				b1=quadrupole_fields2[i+1][j  ][k  ].x-a1;
				c1=quadrupole_fields2[i  ][j+1][k  ].x-a1;
				d1=quadrupole_fields2[i  ][j  ][k+1].x-a1;
				e1=quadrupole_fields2[i+1][j+1][k  ].x-(a1+b1+c1);
				f1=quadrupole_fields2[i+1][j  ][k+1].x-(a1+b1+d1);
				g1=quadrupole_fields2[i  ][j+1][k+1].x-(a1+c1+d1);
				h1=quadrupole_fields2[i+1][j+1][k+1].x-(a1+b1+c1+d1+e1+f1+g1);
				H1.x=(a1+b1*b.x+c1*b.y+d1*b.z+e1*b.x*b.y+f1*b.x*b.z+g1*b.y*b.z+h1*b.x*b.y*b.z)*sign(y1)*ob.parametrs[5];
				a1=quadrupole_fields2[i  ][j  ][k  ].y;
				b1=quadrupole_fields2[i+1][j  ][k  ].y-a1;
				c1=quadrupole_fields2[i  ][j+1][k  ].y-a1;
				d1=quadrupole_fields2[i  ][j  ][k+1].y-a1;
				e1=quadrupole_fields2[i+1][j+1][k  ].y-(a1+b1+c1);
				f1=quadrupole_fields2[i+1][j  ][k+1].y-(a1+b1+d1);
				g1=quadrupole_fields2[i  ][j+1][k+1].y-(a1+c1+d1);
				h1=quadrupole_fields2[i+1][j+1][k+1].y-(a1+b1+c1+d1+e1+f1+g1);
				H1.y=(a1+b1*b.x+c1*b.y+d1*b.z+e1*b.x*b.y+f1*b.x*b.z+g1*b.y*b.z+h1*b.x*b.y*b.z)*sign(x1)*ob.parametrs[5];
				a1=quadrupole_fields2[i  ][j  ][k  ].z;
				b1=quadrupole_fields2[i+1][j  ][k  ].z-a1;
				c1=quadrupole_fields2[i  ][j+1][k  ].z-a1;
				d1=quadrupole_fields2[i  ][j  ][k+1].z-a1;
				e1=quadrupole_fields2[i+1][j+1][k  ].z-(a1+b1+c1);
				f1=quadrupole_fields2[i+1][j  ][k+1].z-(a1+b1+d1);
				g1=quadrupole_fields2[i  ][j+1][k+1].z-(a1+c1+d1);
				h1=quadrupole_fields2[i+1][j+1][k+1].z-(a1+b1+c1+d1+e1+f1+g1);
				H.z=(a1+b1*b.x+c1*b.y+d1*b.z+e1*b.x*b.y+f1*b.x*b.z+g1*b.y*b.z+h1*b.x*b.y*b.z)*sign(x1*y1*z2)*ob.parametrs[5];
				//поворот системы координат обратно на 45 градусов
				H.x = H1.x*cos(-phi)-H1.y*sin(-phi);
				H.y = H1.x*sin(-phi)+H1.y*cos(-phi);
				return H;
			}
		case in_quad3: case in_quad4:
			if (ob.parametrs[3]==0)  //идеальный квадруполь
			{
				phi=ob.parametrs[6]*pi/180;
				//поворот системы координат на 45 градусов
				x1 = (b.x - ob.parametrs[7])*cos(phi)-(b.y - ob.parametrs[8])*sin(phi);
				y1 = (b.x - ob.parametrs[7])*sin(phi)+(b.y - ob.parametrs[8])*cos(phi);
				H1.x=-ob.parametrs[2]*y1;// -ob.parametrs[2]*y;
				H1.y=-ob.parametrs[2]*x1;// -ob.parametrs[2]*x;
				H.z= 0;
				//поворот системы координат обратно на 45 градусов
				H.x = H1.x*cos(-phi)-H1.y*sin(-phi);
				H.y = H1.x*sin(-phi)+H1.y*cos(-phi);
				return H;
			}
			if (ob.parametrs[3]==1)  //квадруполь с гармониками
			{
				fi=atan2(b.y,b.x);
				r_1 = sqrt(b.x*b.x+b.y*b.y);
				r0= 5*fl/p.step;
				A2=An[1];
				for (i=0;i<10;i++)
				{
					An[i]*=ob.parametrs[2]*r0/A2;
					Bn[i]*=ob.parametrs[2]*r0/A2;
				}
				for (i=0;i<10;i++)
				{
					Hr_1 +=pow(r_1/r0,i)*(An[i]*sin((i+1)*fi)+Bn[i]*cos((i+1)*fi));
					Hfi+=pow(r_1/r0,i)*(An[i]*cos((i+1)*fi)-Bn[i]*sin((i+1)*fi));
				}
				H.x=Hr_1*cos(fi+pi/2)-Hfi*sin(fi+pi/2);
				H.y=Hr_1*sin(fi+pi/2)+Hfi*cos(fi+pi/2);
				H.z= 0;
				return H;
			}
			/*
			if (ob.parametrs[3]==2)  //поля из файла.
			{
				z  = z-ob.parametrs[0];
				//поворот квадруполя на 45 градусов
				x1 = x*cos(Pi()/2)-y*sin(Pi()/2);
				y1 = x*sin(Pi()/2)+y*cos(Pi()/2);
				if ((fabs(x1)>quadrupole_parametrs.size_x)||
					(fabs(y1)>quadrupole_parametrs.size_y)||
					(fabs(z )>quadrupole_parametrs.size_z))
				{
					H.x = 0;
					H.y = 0;
					H.z = 0;
					return H;
				}
				i = fabs(x1)/(quadrupole_parametrs.size_dx);
				j = fabs(y1)/(quadrupole_parametrs.size_dy);
				k = fabs(z )/(quadrupole_parametrs.size_dz);
				H.x = quadrupole_fields1[i][j][k].x*sign(y1);
				H.y = quadrupole_fields1[i][j][k].y*sign(x1);
				H.z = quadrupole_fields1[i][j][k].z*sign(x1*y1*z);
				return H;
			}
			*/
		break;
		case in_lens:
			r_1 = sqrt(b.x*b.x+b.y*b.y);
			d = ob.parametrs[1]/2;
			z1 = b.z-ob.parametrs[0];
			if(   z1<0   )
			{
				alpha=-3*(ob.parametrs[2]/d)*(  -pow(z1/d,2)  -  pow(z1/d,1)  );
				betta=-3*(ob.parametrs[2]/(d*d))*( -1-2*(z1/d)   );
				delta=-9* ob.parametrs[2]/(2*d*d*d);
				H.z=ob.parametrs[2]*( -3*pow (z1/d,2)   -2*pow( z1/d,3) + 1 )  +  betta*r_1*r_1/2;
			}
			else
			{
				alpha=-3*(ob.parametrs[2]/d)*(  pow(z1/d,2)  -  pow(z1/d,1)  );
				betta=-3*(ob.parametrs[2]/(d*d))*( -1+2*(z1/d)   );
				delta= 9* ob.parametrs[2]/(2*d*d*d);
				H.z=ob.parametrs[2]*(2*pow(z1/d,3)-3*pow(z1/d,2)+1)  +  betta*r_1*r_1/2;
			}
			H.x = b.x*(alpha  +  delta*r_1*r_1/6);
			H.y = b.y*(alpha  +  delta*r_1*r_1/6);
			return H;
		case in_corr:
			H.x=ob.parametrs[2];
			H.y=ob.parametrs[3];
			H.z=0;
			return H;
		case in_bunch:
			// 1)формируем продольную зависимость в банче
			switch (p.cb)
			{
				case 0://simple
					by_1 = ob.Y_1; by_2 = ob.Y_2;
					bx_1 = ob.X_1; bx_2 = ob.X_2;
					bz = ob.parametrs[0];
					r_1 = sqrt(  pow( b.y - by_1,2)+pow( b.z - bz - bx_1,2));
					r_2 = sqrt(  pow( b.y - by_2,2)+pow( b.z - bz - bx_2,2));
					sigma_r_1=ob.R_1;
					sigma_r_2=ob.R_2;
					switch ((int)ob.parametrs[8])  //выбираем продольное распределение
					{
						case 0://gauss
				//			if (fabs(x-p.b_v*(t-t_ob))<5*ob.parametrs[2])
							//			  nx=exp(    -pow(x-       p.b_v*(t-t_ob),2)/(2*pow(ob.parametrs[2],2))   )/(ob.parametrs[2]*sqrt(2*pi));
								nx_1 = exp(    -pow(b.x-ob.Z_1-p.b_v*(t-t_ob),2)/(2*pow(ob.L_1,2))   )/(ob.L_1*sqrt(2*pi));
								nx_2 = exp(    -pow(b.x-ob.Z_2-p.b_v*(t-t_ob),2)/(2*pow(ob.L_2,2))   )/(ob.L_2*sqrt(2*pi));
				//			else
				//				nx = 0;
							//nx=exp(    -pow(x-p.b_v*(t-t_ob),2)/(2*pow(ob.parametrs[2],2))   )/(ob.parametrs[2]*sqrt(2*pi));
						break;
						case 1://uniform
							if (fabs(b.x-ob.Z_1-p.b_v*(t-t_ob))> (ob.L_1/2)) nx_1 = 0;
							else											 nx_1 = 1.0/ob.L_1;
							if (fabs(b.x-ob.Z_2-p.b_v*(t-t_ob))<=(ob.L_2/2)) nx_2 = 1.0/ob.L_2;
							else											 nx_2 = 0;
						break;
						case 2://triangle
							if (fabs(b.x-ob.Z_1-p.b_v*(t-t_ob))>(ob.L_1/2)) nx_1=0;
							else											nx_1=(1-2*fabs(b.x-ob.Z_1-p.b_v*(t-t_ob))/ob.L_1)/(ob.L_1/2);
							if (fabs(b.x-ob.Z_2-p.b_v*(t-t_ob))>(ob.L_2/2)) nx_2=0;
							else											nx_2=(1-2*fabs(b.x-ob.Z_2-p.b_v*(t-t_ob))/ob.L_2)/(ob.L_1/2);
						break;
						case 3://user non-simmetric 2 different Gauss 
							//if (b.x-p.b_v*(t-t_ob)>=0) nx_1=exp( -pow(b.x-p.b_v*(t-t_ob),2)/(2*pow(ob.parametrs[2]/4,2)) )/(ob.parametrs[2]*sqrt(2*pi));
							//else 					   nx_1=exp( -pow(b.x-p.b_v*(t-t_ob),2)/(2*pow(ob.parametrs[2],2))   )/(ob.parametrs[2]*sqrt(2*pi));
							if (b.x-ob.Z_1-p.b_v*(t-t_ob)>=0) nx_1 = exp(    -pow(b.x-ob.Z_1-p.b_v*(t-t_ob),2)/(2*pow(ob.L_1,2))   )/(ob.L_1*sqrt(2*pi));
							else                       		  nx_1 = exp(    -pow(b.x-ob.Z_1-p.b_v*(t-t_ob),2)/(2*pow(ob.L_2,2))   )/(ob.L_1*sqrt(2*pi));
					}
				break;
				case 1://file
					if (  ((t-t_ob-(b.x/p.b_v))<bunch.t[0]) || ((t-t_ob-(b.x/p.b_v))>bunch.t[bunch.n-1])  ) return H;
					for (i=0;i<bunch.n;i++)
					{
						if (bunch.t[i]>(t-t_ob-b.x/p.b_v))
						{
							ind=i-1;
							break;
						}
					}
					if ((ind==-1)||(ind==bunch.n-2)) ind=0;
					at=(  (t-t_ob-(b.x/p.b_v)) - bunch.t[ind]  )/(  bunch.t[ind+1] - bunch.t[ind]  );
					by_1 = ob.parametrs[4] + bunch.y[ind] + (  bunch.y[ind+1] - bunch.y[ind]  )*at;
					bz = ob.parametrs[0]   + bunch.z[ind] + (  bunch.z[ind+1] - bunch.z[ind]  )*at;
					r_1 = sqrt(  pow( b.y - by_1 ,2) + pow( b.z - bz ,2)  );
					sigma_r_1 = bunch.s[ind] + (  bunch.s[ind+1] - bunch.s[ind]  )*at;
					nx_1      = bunch.i[ind] + (  bunch.i[ind+1] - bunch.i[ind]  )*at;
				break;
				case 2://advanced
					Xt=(t-t_ob)*p.step;
					by_1 = (matrix_make_bunch[2][0][0]*Xt*Xt+matrix_make_bunch[2][1][0]*Xt+matrix_make_bunch[2][2][0])*
						 (matrix_make_bunch[2][3][0]*sin(matrix_make_bunch[2][4][0]*Xt+matrix_make_bunch[2][5][0])+matrix_make_bunch[2][6][0]*exp(-(Xt-matrix_make_bunch[2][7][0])*(Xt-matrix_make_bunch[2][7][0])/(2*matrix_make_bunch[2][8][0]*matrix_make_bunch[2][8][0]))+matrix_make_bunch[2][9][0])*
						 fl/p.step;
					bz = ob.parametrs[0] + 
					     (matrix_make_bunch[2][0][1]*Xt*Xt+matrix_make_bunch[2][1][1]*Xt+matrix_make_bunch[2][2][1])*
					     (matrix_make_bunch[2][3][1]*sin(matrix_make_bunch[2][4][1]*Xt+matrix_make_bunch[2][5][1])+matrix_make_bunch[2][6][1]*exp(-(Xt-matrix_make_bunch[2][7][1])*(Xt-matrix_make_bunch[2][7][1])/(2*matrix_make_bunch[2][8][1]*matrix_make_bunch[2][8][1]))+matrix_make_bunch[2][9][1])*
					     fl/p.step;
					r_1 = sqrt(  pow( b.y - by_1,2) + pow( b.z - bz ,2));
					yz=0;
					sigma_r_1 = (matrix_make_bunch[1][0][yz]*Xt*Xt+matrix_make_bunch[1][1][yz]*Xt+matrix_make_bunch[1][2][yz])*
					          (matrix_make_bunch[1][3][yz]*sin(matrix_make_bunch[1][4][yz]*Xt+matrix_make_bunch[1][5][yz])+matrix_make_bunch[1][6][yz]*exp(-(Xt-matrix_make_bunch[1][7][yz])*(Xt-matrix_make_bunch[1][7][yz])/(2*matrix_make_bunch[1][8][yz]*matrix_make_bunch[1][8][yz]))+matrix_make_bunch[1][9][yz])*
					          fl/p.step;
					nx_1 = (Normirovka*p.step)*(matrix_make_bunch[0][0][0]*Xt*Xt+matrix_make_bunch[0][1][0]*Xt+matrix_make_bunch[0][2][0])*
					     (matrix_make_bunch[0][3][0]*sin(matrix_make_bunch[0][4][0]*Xt+matrix_make_bunch[0][5][0])+matrix_make_bunch[0][6][0]*exp(-(Xt-matrix_make_bunch[0][7][0])*(Xt-matrix_make_bunch[0][7][0])/(2*matrix_make_bunch[0][8][0]*matrix_make_bunch[0][8][0]))+matrix_make_bunch[0][9][0]);
			}
			switch ((int)ob.TRD_1) //bunch 1 считаем радиальное поле
			{
				case 0://gauss
					Hr_1 = ob.K_1*ob.parametrs[7]*(fe*fl*fl/p.step)*p.b_v*ob.parametrs[6]*2*qe*(1 - exp(    -r_1*r_1/(2*pow(sigma_r_1,2))    ))*nx_1/r_1;
				break;
				case 1://uniform
					if (r_1>sigma_r_1) Hr_1 = ob.K_1*ob.parametrs[7]*(fe*fl*fl/p.step)*p.b_v*ob.parametrs[6]*2*qe*(  1  )*nx_1/r_1;
					else			   Hr_1 = ob.K_1*ob.parametrs[7]*(fe*fl*fl/p.step)*p.b_v*ob.parametrs[6]*2*qe*(  r_1*r_1/(sigma_r_1*sigma_r_1)  )*nx_1/r_1;
				break;
				case 2://triangle
				
					if (r_1>sigma_r_1) Hr_1 = ob.K_1*ob.parametrs[7]*(fe*fl*fl/p.step)*p.b_v*ob.parametrs[6]*2*qe*(  1  )*nx_1/r_1;
					else               Hr_1 = ob.K_1*ob.parametrs[7]*(fe*fl*fl/p.step)*p.b_v*ob.parametrs[6]*2*qe*(  3*r_1*r_1/(sigma_r_1*sigma_r_1)-2*r_1*r_1*r_1/(pow(sigma_r_1,3))  )*nx_1/r_1;
				break;
				case 4: //flat horizontal
					for (i=1;i<=p.N_Flat;i++)
					{
						dx = (2*i-p.N_Flat-1)*2*sigma_r_1/2;
						r1 = sqrt(  pow( b.y - by_1,2)+pow( b.z - bz - bx_1 - dx,2));
						Hr = ob.K_1*ob.parametrs[7]*(fe*fl*fl/p.step)*p.b_v*ob.parametrs[6]*2*qe*(1 - exp(-r1*r1/(2*pow(sigma_r_1,2)))  )*nx_1/r1;
						Hr*= exp( -pow((2*i-p.N_Flat-1)*sigma_r_1,2)/(2*pow(sigma_r_1*10,2))); // gauss flat
						H1.y += Hr * (b.z-bz-bx_1-dx)/r1;
						H1.z -= Hr * (b.y-by_1)/r1;
					}
					H.y += H1.y/p.N_Flat;
					H.z += H1.z/p.N_Flat;
				break;
				case 5: //flat vertical
					for (i=1;i<=p.N_Flat;i++)
					{
						dx = (2*i-p.N_Flat-1)*2*sigma_r_1/2;
						r1 = sqrt(  pow( b.y - by_1 - dx,2)+pow( b.z - bz - bx_1,2));
						Hr = ob.K_1*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(1 - exp(-r1*r1/(2*pow(sigma_r_1,2)))  )*nx_1/r1;
						Hr*= exp( -pow((2*i-p.N_Flat-1)*sigma_r_1,2)/(2*pow(sigma_r_1*10,2))); // gauss flat
						H1.y += Hr * (b.z-bz-bx_1)/r1;
						H1.z -= Hr * (b.y-by_1-dx)/r1;
					}
					H.y += H1.y/p.N_Flat;
					H.z += H1.z/p.N_Flat;
				break;
			}
			switch ((int)ob.TRD_2) //bunch 2
			{
				case 0://gauss
					Hr_2 = ob.K_2*ob.parametrs[7]*(fe*fl*fl/p.step)*p.b_v*ob.parametrs[6]*2*qe*(1 - exp(    -r_2*r_2/(2*pow(sigma_r_2,2))    ))*nx_2/r_2;
				break;
				case 1://uniform
					if (r_2>sigma_r_2) Hr_2 = ob.K_2*ob.parametrs[7]*(fe*fl*fl/p.step)*p.b_v*ob.parametrs[6]*2*qe*(  1  )*nx_2/r_2;
					else               Hr_2 = ob.K_2*ob.parametrs[7]*(fe*fl*fl/p.step)*p.b_v*ob.parametrs[6]*2*qe*(  r_2*r_2/(sigma_r_2*sigma_r_2)  )*nx_2/r_2;
				break;
				case 2://triangle
					if (r_2>sigma_r_2) Hr_2 = ob.K_2*ob.parametrs[7]*(fe*fl*fl/p.step)*p.b_v*ob.parametrs[6]*2*qe*(  1  )*nx_2/r_2;
					else			   Hr_2 = ob.K_2*ob.parametrs[7]*(fe*fl*fl/p.step)*p.b_v*ob.parametrs[6]*2*qe*(  3*r_2*r_2/(sigma_r_2*sigma_r_2)-2*r_2*r_2*r_2/(pow(sigma_r_2,3))  )*nx_2/r_2;
				break;
				case 4: //flat horizontal
					for (i=1;i<=p.N_Flat;i++)
					{
						dx = (2*i-p.N_Flat-1)*2*sigma_r_2/2;
						r2 = sqrt(  pow( b.y - by_2,2)+pow( b.z - bz - bx_2 - dx,2));
						Hr = ob.K_2*ob.parametrs[7]*(fe*fl*fl/p.step)*p.b_v*ob.parametrs[6]*2*qe*(1 - exp(-r2*r2/(2*pow(sigma_r_2,2)))  )*nx_2/r2;
						Hr*= exp( -pow((2*i-p.N_Flat-1)*sigma_r_2,2)/(2*pow(sigma_r_2*10,2))); // gauss flat
						H2.y += Hr * (b.z-bz-bx_2-dx)/r2;
						H2.z -= Hr * (b.y-by_2)/r2;
					}
					H.y += H2.y/p.N_Flat;
					H.z += H2.z/p.N_Flat;
				break;
				case 5: //flat vertical
					for (i=1;i<=p.N_Flat;i++)
					{
						dx = (2*i-p.N_Flat-1)*2*sigma_r_2/2;
						r2 = sqrt(  pow( b.y - by_2 - dx,2)+pow( b.z - bz - bx_2,2));
						Hr = ob.K_2*ob.parametrs[7]*(fe*fl*fl/p.step)*ob.parametrs[6]*2*qe*(1 - exp(-r2*r2/(2*pow(sigma_r_2,2)))  )*nx_2/r2;
						Hr*= exp( -pow((2*i-p.N_Flat-1)*sigma_r_2,2)/(2*pow(sigma_r_2*10,2))); // gauss flat
						H2.y += Hr * (b.z-bz-bx_2)/r2;
						H2.z -= Hr * (b.y-by_2-dx)/r2;
					}
					H.y += H2.y/p.N_Flat;
					H.z += H2.z/p.N_Flat;
				break;
			}
			H.y +=  ((b.z-bx_1-bz)/  r_1)*Hr_1+((b.z-bx_2-bz)/  r_2)*Hr_2;
			H.z += -((b.y-by_1   )/  r_1)*Hr_1-((b.y-by_2   )/  r_2)*Hr_2;
			return H;
		//case in_RF:
		//if (ob.parametrs[5]==0)   //ideal
		//if (ob.parametrs[5]==1)   //real
		//if (ob.parametrs[5]==2)   //file
		//case in_ACC:
		//break;
	}
	return H;
}

struct vector E_space_charge(int k, struct particle *b, struct parametrs p)
{
	int i; double gamma;
	struct vector E={0};
	gamma=pow(1+pow(b[k].px,2)+pow(b[k].py,2)+pow(b[k].pz,2),0.5);
	for (i=0;i<p.t_number;i++)
	{
		/*if ((i!=k)&&(b[i].status!=on_tube)&&(b[i].status!=on_screen)&&(b[i].status!=reflection)&&(b[i].status!=on_diaphragm)&&(b[i].status!=in_cathode))
		{
			E.x+=gamma*qe*(b[k].x-b[i].x)/pow(pow(b[i].x-b[k].x,2)+pow(b[i].y-b[k].y,2)+pow(b[i].z-b[k].z,2),3/2);
			E.y+=gamma*qe*(b[k].y-b[i].y)/pow(pow(b[i].x-b[k].x,2)+pow(b[i].y-b[k].y,2)+pow(b[i].z-b[k].z,2),3/2);
			E.z+=      qe*(b[k].z-b[i].z)/pow(pow(b[i].x-b[k].x,2)+pow(b[i].y-b[k].y,2)+pow(b[i].z-b[k].z,2),3/2);
		}*/
		switch (b[i].status)
		{
			case in_cathode: case on_tube: case on_screen: case on_diaphragm: case reflection: case on_stop:
			continue;
		}
		if (i!=k)
		{			
			E.x+=gamma*qe*(b[k].x-b[i].x)/pow(pow(b[i].x-b[k].x,2)+pow(b[i].y-b[k].y,2)+pow(b[i].z-b[k].z,2),3/2);
			E.y+=gamma*qe*(b[k].y-b[i].y)/pow(pow(b[i].x-b[k].x,2)+pow(b[i].y-b[k].y,2)+pow(b[i].z-b[k].z,2),3/2);
			E.z+=      qe*(b[k].z-b[i].z)/pow(pow(b[i].x-b[k].x,2)+pow(b[i].y-b[k].y,2)+pow(b[i].z-b[k].z,2),3/2);
		}
	}
	E.x*=-p.t_weight*fl*fl*fe/p.step;
	E.y*=-p.t_weight*fl*fl*fe/p.step;
	E.z*=-p.t_weight*fl*fl*fe/p.step;
	return E;
}

struct vector H_space_charge(int k, struct particle *b, struct parametrs p)
{
	int i; double gamma,betta;
	struct vector H={0};
	gamma = pow(1+pow( b[k].px,2 ) + pow( b[k].py,2 ) + pow( b[k].pz,2 ),0.5);
	betta = b[k].pz/gamma;
	for (i=0;i<p.t_number;i++)
	{
		/*if ((i!=k)&&(b[i].status!=on_tube)&&(b[i].status!=on_screen)&&(b[i].status!=reflection)&&(b[i].status!=on_diaphragm)&&(b[i].status!=in_cathode))
		{
			H.x-=gamma*betta*qe*(b[k].y-b[i].y)/pow(pow(b[i].x-b[k].x,2)+pow(b[i].y-b[k].y,2)+pow(b[i].z-b[k].z,2),3/2);
			H.y+=gamma*betta*qe*(b[k].x-b[i].x)/pow(pow(b[i].x-b[k].x,2)+pow(b[i].y-b[k].y,2)+pow(b[i].z-b[k].z,2),3/2);
		}*/
		switch (b[i].status)
		{
			case in_cathode: case on_tube: case on_screen: case on_diaphragm: case reflection: case on_stop:
			continue;
		}
		if (i!=k)
		{			
			H.x-=gamma*betta*qe*(b[k].y-b[i].y)/pow(pow(b[i].x-b[k].x,2)+pow(b[i].y-b[k].y,2)+pow(b[i].z-b[k].z,2),3/2);
			H.y+=gamma*betta*qe*(b[k].x-b[i].x)/pow(pow(b[i].x-b[k].x,2)+pow(b[i].y-b[k].y,2)+pow(b[i].z-b[k].z,2),3/2);
		}
	}
	H.x*=-p.t_weight*fl*fl*fe/p.step;
	H.y*=-p.t_weight*fl*fl*fe/p.step;
	return H;
}

struct parametrs par_to_in(struct parametrs p, struct object* ob)
{
	p.t_v=sqrt((p.t_energy*1e3/me)*(2+(p.t_energy*1e3/me)))/(1+(p.t_energy*1e3/me)); // betta 
	switch ((int)ob[7].parametrs[7])
	{
		case -1:  p.b_v = sqrt(   (  ob[7].parametrs[5]*1e9/me  )*(  2+(ob[7].parametrs[5]*1e9/me)  )   ) / (1+(ob[7].parametrs[5]*1e9/me)); break;
		case 1:   p.b_v = sqrt(   (  ob[7].parametrs[5]*1e9/mp  )*(  2+(ob[7].parametrs[5]*1e9/mp)  )   ) / (1+(ob[7].parametrs[5]*1e9/mp));
	}
	p.t_radius*=fl/p.step;
	p.gun_length*=fl/p.step;
	p.t_length*=p.t_v/p.step;
	p.radius*=fl/p.step;
	p.screen*=fl/p.step;
	p.stop*=fl/p.step;
	p.diaph*=fl/p.step;
	p.diaph_r*=fl/p.step;
	p.x_min*=fl/p.step;
	p.x_max*=fl/p.step;
	p.y_min*=fl/p.step;
	p.y_max*=fl/p.step;
	p.step;
	return p;
}

struct object* obj_to_in(struct parametrs p, struct object* ob)
{
	int i;
	struct object* re_obj;
	re_obj=calloc(number_of_object,sizeof(struct object));
	for (i=0;i<number_of_object;i++)
	{
		re_obj[i].parametrs[0]=ob[i].parametrs[0]*fl/p.step;
		re_obj[i].parametrs[1]=ob[i].parametrs[1]*fl/p.step;
		re_obj[i].name        =ob[i].name;
		// re_obj[i].E            =ob[i].E; // To be realized
		// re_obj[i].H			  =ob[i].H; // To be realized
	}
	re_obj[0].parametrs[2]=ob[0].parametrs[2]*fl/p.step;					//положение первого квадруполя
	re_obj[0].parametrs[3]=ob[0].parametrs[3]*fl/p.step;					//положение второго квадруполя
	re_obj[0].parametrs[4]=ob[0].parametrs[4];								//сила первого квадруполя
	re_obj[0].parametrs[5]=ob[0].parametrs[5];								//сила второго квадруполя
	re_obj[0].parametrs[6]=ob[0].parametrs[6];								//угол системы квадруполей
	re_obj[0].parametrs[7]=ob[0].parametrs[7]*fl/p.step;					//положение x
	re_obj[0].parametrs[8]=ob[0].parametrs[8]*fl/p.step;					//положение y
	for (i=2;i<4;i++) 														//квадруполи 3 и 4 
	{
		re_obj[i].parametrs[2]=ob[i].parametrs[2]*1000*fe*p.step*p.step/fl; //градиент в квадруполе
		re_obj[i].parametrs[3]=ob[i].parametrs[3];                          //идеальность квадруполя
		re_obj[i].parametrs[4]=ob[i].parametrs[4];                          //угол квадруполя
		re_obj[i].parametrs[7]=ob[i].parametrs[7]*fl/p.step;				//положение x
		re_obj[i].parametrs[8]=ob[i].parametrs[8]*fl/p.step;				//положение y
	}
	re_obj[4].parametrs[2]=ob[4].parametrs[2]*1000*fe*p.step;			    //Solenoid field
	re_obj[5].parametrs[2]=ob[5].parametrs[2]*fl/p.step;					//scan plates distance
	re_obj[5].parametrs[3]=ob[5].parametrs[3]*100.0*fe*fl/30000.0;			//scan Ux	   (максимальное)
	re_obj[5].parametrs[4]=ob[5].parametrs[4]*100.0*fe*fl/30000.0;			//scan Uy	   (максимальное)
	re_obj[5].parametrs[5]=ob[5].parametrs[5]*1000.0/p.step;				//scan duration [ns]
	re_obj[6].parametrs[2]=ob[6].parametrs[2]*100000.0*fe*fl/30000.0;		//поле в резонаторе
	re_obj[6].parametrs[3]=ob[6].parametrs[3]/p.step;						//период ВЧ
	re_obj[6].parametrs[4]=ob[6].parametrs[4]*fl/p.step;					//длина волны
	re_obj[9].parametrs[2]=ob[9].parametrs[2]*1.0e8*fe*fl/30000.0;		    //acceleration field
	re_obj[9].parametrs[3]=ob[9].parametrs[3]/p.step;						//период ВЧ
	re_obj[9].parametrs[4]=ob[9].parametrs[4]*fl/p.step;					//длина волны
	//параметры исследуемого пучка
	re_obj[7].parametrs[2]=ob[7].parametrs[2]*fl/p.step; 					//sigma_l_1 + sigma_l_2
	//re_obj[7].parametrs[4]=ob[7].parametrs[4]*fl/p.step; 					//bunch_y
	re_obj[7].parametrs[5]=ob[7].parametrs[5];
	re_obj[7].parametrs[6]=ob[7].parametrs[6];
	re_obj[7].parametrs[7]=ob[7].parametrs[7];
	re_obj[7].parametrs[8]=ob[7].parametrs[8];
	//re_obj[7].parametrs[9]=ob[7].parametrs[9];
	re_obj[7].X_1=ob[7].X_1*fl/p.step; 										//bunch_x_1
	re_obj[7].X_2=ob[7].X_2*fl/p.step; 										//bunch_x_2
	re_obj[7].Y_1=ob[7].Y_1*fl/p.step; 										//bunch_y_1
	re_obj[7].Y_2=ob[7].Y_2*fl/p.step; 										//bunch_y_2
	re_obj[7].Z_1=ob[7].Z_1*fl/p.step; 										//bunch_z_1
	re_obj[7].Z_2=ob[7].Z_2*fl/p.step; 										//bunch_z_2
	re_obj[7].L_1=ob[7].L_1*fl/p.step; 										//bunch_z_1
	re_obj[7].L_2=ob[7].L_2*fl/p.step; 										//bunch_z_2
	re_obj[7].R_1=ob[7].R_1*fl/p.step; 										//bunch_sigma_r_1
	re_obj[7].R_2=ob[7].R_2*fl/p.step; 										//bunch_sigma_r_2
	re_obj[7].K_1=ob[7].K_1/100;
	re_obj[7].K_2=ob[7].K_2/100;
	re_obj[7].TRD_1=ob[7].TRD_1;
	re_obj[7].TRD_2=ob[7].TRD_2;
	re_obj[8].parametrs[2]=ob[8].parametrs[2]*fe*p.step; 					//поле в корректоре
	re_obj[8].parametrs[3]=ob[8].parametrs[3]*fe*p.step; 					//поле в корректоре
	return re_obj;
}
// ================================ Important function and quite time consuming ===================== //
int where_is_particle(struct particle b, struct parametrs p, struct object *ob)
{	
	int i; double r = sqrt(b.x*b.x + b.y*b.y);
	if      (    r>=p.radius ) return on_tube;
	else if (  b.z>=p.stop   ) return on_stop;
	else if (  b.z>=p.screen ) return on_screen;
	else if ( (b.z>=p.diaph)&&(b.z<p.diaph+1)&&(r>=p.diaph_r)	) return on_diaphragm;
	for ( i=0; i<number_of_object-5; i++) if (  fabs(b.z-ob[i].parametrs[0]) <= ob[i].parametrs[1]/2  ) return ob[i].name;
	return in_free_space;
}
// ====================================== Draws beamline and components ============================= //
void draw_beam(struct object *ob)						  
{
	int i_obj;
	switch(par.graph_mode)
	{
		case 1: // All particles (x-z)
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE);
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_FILL_COLOR, VAL_BLACK);
		//рисуем трубу
			CanvasDrawRect (panelS, PANEL_S_CANVAS_2, MakeRect(50-par.a_y*500*par.radius/par.screen,0,par.a_y*2*500*par.radius/par.screen+1,500), VAL_DRAW_FRAME_AND_INTERIOR);
		//рисуем резонатор
			CanvasDrawRect (panelS, PANEL_S_CANVAS_2, 
				MakeRect(50-par.a_y*30,500*(ob[6].parametrs[0]-ob[6].parametrs[1]/2)/par.screen,par.a_y*2*30+1,500*ob[6].parametrs[1]/par.screen+1), VAL_DRAW_FRAME_AND_INTERIOR);
		//рисуем линзу
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_RED);
			CanvasDrawRect (panelS, PANEL_S_CANVAS_2, 
				MakeRect(50-par.a_y*30,500*(ob[4].parametrs[0]-ob[4].parametrs[1]/2)/par.screen,par.a_y*2*30+1,500*(ob[4].parametrs[1])/par.screen+1), VAL_DRAW_FRAME);
		//рисуем квадруполь
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_WIDTH, 3);
		//первый
			if (ob[0].parametrs[4]<0)  SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_RED);
			if (ob[0].parametrs[4]>0)  SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_BLUE);
			if (ob[0].parametrs[4]==0) SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE);
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, 
				MakePoint (500*ob[0].parametrs[2]/par.screen, 3), 
				MakePoint (500*ob[0].parametrs[2]/par.screen, 50-par.a_y*500*par.radius/par.screen-1));
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, 
				MakePoint (500*ob[0].parametrs[2]/par.screen, 97), 
				MakePoint (500*ob[0].parametrs[2]/par.screen, 50+par.a_y*500*par.radius/par.screen+3));
		//второй
			if (ob[0].parametrs[5]<0)  SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_RED);
			if (ob[0].parametrs[5]>0)  SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_BLUE);
			if (ob[0].parametrs[5]==0) SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE);
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, 
				MakePoint (500*ob[0].parametrs[3]/par.screen, 3), 
				MakePoint (500*ob[0].parametrs[3]/par.screen, 50-par.a_y*500*par.radius/par.screen-1));
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, 
				MakePoint (500*ob[0].parametrs[3]/par.screen, 97), 
				MakePoint (500*ob[0].parametrs[3]/par.screen, 50+par.a_y*500*par.radius/par.screen+3));
		//третий и четвёртый
			for (i_obj=2;i_obj<4;i_obj++)
			{
				if (ob[i_obj].parametrs[2]<0)  SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_RED);
				if (ob[i_obj].parametrs[2]>0)  SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_BLUE);
				if (ob[i_obj].parametrs[2]==0) SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE);
				CanvasDrawLine (panelS, PANEL_S_CANVAS_2, 
					MakePoint (500*ob[i_obj].parametrs[0]/par.screen, 3), 
					MakePoint (500*ob[i_obj].parametrs[0]/par.screen, 50-par.a_y*500*par.radius/par.screen-1));
				CanvasDrawLine (panelS, PANEL_S_CANVAS_2, 
					MakePoint (500*ob[i_obj].parametrs[0]/par.screen, 97), 
					MakePoint (500*ob[i_obj].parametrs[0]/par.screen, 50+par.a_y*500*par.radius/par.screen+3));
			}
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_WIDTH, 1);
		//рисуем область взаимодействия
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE);
			CanvasDrawRect (panelS, PANEL_S_CANVAS_2, 
				MakeRect(-1,500*((ob[7].parametrs[0]-(ob[7].parametrs[1]/2))/par.screen),
					102,500*2*(ob[7].parametrs[1]/2)/par.screen+1), VAL_DRAW_FRAME_AND_INTERIOR);
		//стираем лишнее
		  //CanvasDrawRect (panelS, PANEL_S_CANVAS_2, MakeRect(50+2-par.a_y*500*par.radius/par.screen,0,par.a_y*2*500*par.radius/par.screen+3-4,500), VAL_DRAW_INTERIOR);
			CanvasDrawRect (panelS, PANEL_S_CANVAS_2, MakeRect(50+1-par.a_y*500*par.radius/par.screen,0,par.a_y*2*500*par.radius/par.screen+1-2,500), VAL_DRAW_INTERIOR);
		//рисуем пластины развёртки
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, MakePoint (500*(ob[5].parametrs[0]-ob[5].parametrs[1]/2)/par.screen, 50-par.a_y*250*ob[5].parametrs[2]/par.screen), 
										  			  MakePoint (500*(ob[5].parametrs[0]+ob[5].parametrs[1]/2)/par.screen, 50-par.a_y*250*ob[5].parametrs[2]/par.screen));
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, MakePoint (500*(ob[5].parametrs[0]-ob[5].parametrs[1]/2)/par.screen, 50+par.a_y*250*ob[5].parametrs[2]/par.screen), 
										  			  MakePoint (500*(ob[5].parametrs[0]+ob[5].parametrs[1]/2)/par.screen, 50+par.a_y*250*ob[5].parametrs[2]/par.screen));
		//рисуем диафрагму
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, MakePoint (500*par.diaph/par.screen, 50-par.a_y*500*par.radius /par.screen), 
										  			  MakePoint (500*par.diaph/par.screen, 50-par.a_y*500*par.diaph_r/par.screen));
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, MakePoint (500*par.diaph/par.screen, 50+par.a_y*500*par.radius /par.screen), 
										  			  MakePoint (500*par.diaph/par.screen, 50+par.a_y*500*par.diaph_r/par.screen));
		//рисуем банч
			if (ob[7].parametrs[6]!=0)
			{
				SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_MAGENTA);
				SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_FILL_COLOR, VAL_MAGENTA);
				CanvasDrawOval (panelS, PANEL_S_CANVAS_2, 
							MakeRect(50-par.a_y*500* (ob[7].L_1*10+ob[7].Z_1*10)/par.screen, 500* (ob[7].parametrs[0]-ob[7].R_1*10+ob[7].X_1*10)/par.screen, 
								  		par.a_y*500*2*ob[7].L_1*10/par.screen , 500*2*ob[7].R_1*10/par.screen), VAL_DRAW_FRAME_AND_INTERIOR);
				CanvasDrawOval (panelS, PANEL_S_CANVAS_2, 
							MakeRect(50-par.a_y*500* (ob[7].L_2*10+ob[7].Z_2*10)/par.screen, 500* (ob[7].parametrs[0]-ob[7].R_2*10+ob[7].X_2*10)/par.screen, 
								  		par.a_y*500*2*ob[7].L_2*10/par.screen , 500*2*ob[7].R_2*10/par.screen), VAL_DRAW_FRAME_AND_INTERIOR);
			}
		break;
		case 2:
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE);
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_FILL_COLOR, VAL_BLACK);
			//рисуем трубу
			CanvasDrawRect (panelS, PANEL_S_CANVAS_2, MakeRect(50-par.a_y*500*par.radius/par.screen,0,par.a_y*2*500*par.radius/par.screen+1,500), VAL_DRAW_FRAME_AND_INTERIOR);
			//рисуем резонатор
			CanvasDrawRect (panelS, PANEL_S_CANVAS_2, 
				MakeRect(50-par.a_y*30,500*(ob[6].parametrs[0]-ob[6].parametrs[1]/2)/par.screen,par.a_y*2*30+1,500*(ob[6].parametrs[1])/par.screen+1), VAL_DRAW_FRAME_AND_INTERIOR);
			//рисуем линзу
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_RED);
			CanvasDrawRect (panelS, PANEL_S_CANVAS_2, 
				MakeRect(50-par.a_y*30,500*(ob[4].parametrs[0]-ob[4].parametrs[1]/2)/par.screen,par.a_y*2*30+1,500*(ob[4].parametrs[1])/par.screen+1), VAL_DRAW_FRAME);
			//рисуем квадруполь
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_WIDTH, 3);
			//первый
			if (ob[0].parametrs[4]<0)  SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_RED);
			if (ob[0].parametrs[4]>0)  SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_BLUE);
			if (ob[0].parametrs[4]==0) SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE);
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, 
				MakePoint (500*ob[0].parametrs[2]/par.screen, 3), 
				MakePoint (500*ob[0].parametrs[2]/par.screen, 50-par.a_y*500*par.radius/par.screen-1));
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, 
				MakePoint (500*ob[0].parametrs[2]/par.screen, 97), 
				MakePoint (500*ob[0].parametrs[2]/par.screen, 50+par.a_y*500*par.radius/par.screen+3));
			//второй
			if (ob[0].parametrs[5]<0)  SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_RED);
			if (ob[0].parametrs[5]>0)  SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_BLUE);
			if (ob[0].parametrs[5]==0) SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE);
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, 
				MakePoint (500*ob[0].parametrs[3]/par.screen, 3), 
				MakePoint (500*ob[0].parametrs[3]/par.screen, 50-par.a_y*500*par.radius/par.screen-1));
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, 
				MakePoint (500*ob[0].parametrs[3]/par.screen, 97), 
				MakePoint (500*ob[0].parametrs[3]/par.screen, 50+par.a_y*500*par.radius/par.screen+3));
			//третий и четвёртый
			for (i_obj=2;i_obj<4;i_obj++)
			{
				if (ob[i_obj].parametrs[2]<0)  SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_RED);
				if (ob[i_obj].parametrs[2]>0)  SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_BLUE);
				if (ob[i_obj].parametrs[2]==0) SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE);
				CanvasDrawLine (panelS, PANEL_S_CANVAS_2, 
					MakePoint (500*ob[i_obj].parametrs[0]/par.screen, 3), 
					MakePoint (500*ob[i_obj].parametrs[0]/par.screen, 50-par.a_y*500*par.radius/par.screen-1));
				CanvasDrawLine (panelS, PANEL_S_CANVAS_2, 
					MakePoint (500*ob[i_obj].parametrs[0]/par.screen, 97), 
					MakePoint (500*ob[i_obj].parametrs[0]/par.screen, 50+par.a_y*500*par.radius/par.screen+3));
			}
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_WIDTH, 1);
			//рисуем область взаимодействия
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE);
			CanvasDrawOval (panelS, PANEL_S_CANVAS_2, 
				MakeRect(50-par.a_y*500*( (ob[7].parametrs[1]/2)/par.screen), 500*( (ob[7].parametrs[0]-(ob[7].parametrs[1]/2))/par.screen),
							par.a_y*500*2*(ob[7].parametrs[1]/2)/par.screen+1,500*2*(ob[7].parametrs[1]/2)/par.screen+1), VAL_DRAW_FRAME_AND_INTERIOR);
			//стираем лишнее
			//CanvasDrawRect (panelS, PANEL_S_CANVAS_2, MakeRect(50+2-par.a_y*500*par.radius/par.screen,0,par.a_y*2*500*par.radius/par.screen+3-4,500), VAL_DRAW_INTERIOR);
			CanvasDrawRect (panelS, PANEL_S_CANVAS_2, MakeRect(50+1-par.a_y*500*par.radius/par.screen,0,par.a_y*2*500*par.radius/par.screen+1-2,500), VAL_DRAW_INTERIOR);
			//рисуем пластины развёртки
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_FILL_COLOR, VAL_WHITE);
			CanvasDrawRect (panelS, PANEL_S_CANVAS_2, 
				MakeRect (50-par.a_y*250*ob[5].parametrs[2]/par.screen+0, 500*(ob[5].parametrs[0]-ob[5].parametrs[1]/2)/par.screen,
							 par.a_y*500*ob[5].parametrs[2]/par.screen+1, 500* ob[5].parametrs[1]/par.screen+1), 
							   VAL_DRAW_FRAME_AND_INTERIOR); 
			//рисуем диафрагму
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, MakePoint (500*par.diaph/par.screen, 50-par.a_y*500*par.radius/par.screen), 
										  			  MakePoint (500*par.diaph/par.screen, 50-par.a_y*500*par.diaph_r/par.screen));
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, MakePoint (500*par.diaph/par.screen, 50+par.a_y*500*par.radius/par.screen), 
										  			  MakePoint (500*par.diaph/par.screen, 50+par.a_y*500*par.diaph_r/par.screen));
			//рисуем банч
			if ((ob[7].parametrs[6]!=0) && (1))
			{
				SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_MAGENTA);
				SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_FILL_COLOR, VAL_MAGENTA);
				CanvasDrawOval (panelS, PANEL_S_CANVAS_2, 
					MakeRect(50-par.a_y*500*( ob[7].R_1*1+ob[7].Y_1)/par.screen, 500*(ob[7].parametrs[0]-ob[7].R_1*1)/par.screen,
								par.a_y*500*2*ob[7].R_1*1/par.screen,			     500*2*ob[7].R_1*10/par.screen), VAL_DRAW_FRAME_AND_INTERIOR);
				CanvasDrawOval (panelS, PANEL_S_CANVAS_2, 
					MakeRect(50-par.a_y*500*( ob[7].R_2*1+ob[7].Y_2)/par.screen, 500*(ob[7].parametrs[0]-ob[7].R_2*1)/par.screen,
								par.a_y*500*2*ob[7].R_2*1/par.screen,			     500*2*ob[7].R_2*10/par.screen), VAL_DRAW_FRAME_AND_INTERIOR);
			}
		break;
		case 3: case 4:
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE);
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_FILL_COLOR, VAL_BLACK);
			//рисуем трубу
			CanvasDrawRect (panelS, PANEL_S_CANVAS_2, MakeRect(50-par.a_y*500*par.radius/par.screen,0,par.a_y*2*500*par.radius/par.screen+1,500), VAL_DRAW_FRAME_AND_INTERIOR);
			//рисуем резонатор
			CanvasDrawRect (panelS, PANEL_S_CANVAS_2, 
				MakeRect(50-par.a_y*30,500*(ob[6].parametrs[0]-ob[6].parametrs[1]/2)/par.screen,par.a_y*2*30+1,500*ob[6].parametrs[1]/par.screen+1), VAL_DRAW_FRAME_AND_INTERIOR);
			//рисуем линзу
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_RED);
			CanvasDrawRect (panelS, PANEL_S_CANVAS_2, 
				MakeRect(50-par.a_y*30,500*(ob[4].parametrs[0]-ob[4].parametrs[1]/2)/par.screen,par.a_y*2*30+1,500*(ob[4].parametrs[1])/par.screen+1), VAL_DRAW_FRAME);
			//рисуем квадруполь
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_WIDTH, 3);
			//первый // so dificult because 1 and 2 are double quadrupole
			if (ob[0].parametrs[4]<0)  SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_RED);
			if (ob[0].parametrs[4]>0)  SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_BLUE);
			if (ob[0].parametrs[4]==0) SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE);
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, 
				MakePoint (500*ob[0].parametrs[2]/par.screen, 3), 
				MakePoint (500*ob[0].parametrs[2]/par.screen, 50-par.a_y*500*par.radius/par.screen-1));
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, 
				MakePoint (500*ob[0].parametrs[2]/par.screen, 97), 
				MakePoint (500*ob[0].parametrs[2]/par.screen, 50+par.a_y*500*par.radius/par.screen+3));
			//второй
			if (ob[0].parametrs[5]<0)  SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_RED);
			if (ob[0].parametrs[5]>0)  SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_BLUE);
			if (ob[0].parametrs[5]==0) SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE);
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, 
				MakePoint (500*ob[0].parametrs[3]/par.screen, 3), 
				MakePoint (500*ob[0].parametrs[3]/par.screen, 50-par.a_y*500*par.radius/par.screen-1));
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, 
				MakePoint (500*ob[0].parametrs[3]/par.screen, 97), 
				MakePoint (500*ob[0].parametrs[3]/par.screen, 50+par.a_y*500*par.radius/par.screen+3));
			//третий и четвёртый
			for (i_obj=2;i_obj<4;i_obj++)
			{
				if (ob[i_obj].parametrs[2]<0)  SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_RED);
				if (ob[i_obj].parametrs[2]>0)  SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_BLUE);
				if (ob[i_obj].parametrs[2]==0) SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE);
				CanvasDrawLine (panelS, PANEL_S_CANVAS_2, 
					MakePoint (500*ob[i_obj].parametrs[0]/par.screen, 3), 
					MakePoint (500*ob[i_obj].parametrs[0]/par.screen, 50-par.a_y*500*par.radius/par.screen-1));
				CanvasDrawLine (panelS, PANEL_S_CANVAS_2, 
					MakePoint (500*ob[i_obj].parametrs[0]/par.screen, 97), 
					MakePoint (500*ob[i_obj].parametrs[0]/par.screen, 50+par.a_y*500*par.radius/par.screen+3));
			}
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_WIDTH, 1);
			//стираем лишнее
			CanvasDrawRect (panelS, PANEL_S_CANVAS_2, MakeRect(50+1-par.a_y*500*par.radius/par.screen,0,par.a_y*2*500*par.radius/par.screen+1-2,500), VAL_DRAW_INTERIOR);
			//рисуем пластины развёртки
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE);
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, MakePoint (500*(ob[5].parametrs[0]-ob[5].parametrs[1]/2)/par.screen, 50-par.a_y*250*ob[5].parametrs[2]/par.screen), 
										  			  MakePoint (500*(ob[5].parametrs[0]+ob[5].parametrs[1]/2)/par.screen, 50-par.a_y*250*ob[5].parametrs[2]/par.screen));
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, MakePoint (500*(ob[5].parametrs[0]-ob[5].parametrs[1]/2)/par.screen, 50+par.a_y*250*ob[5].parametrs[2]/par.screen), 
										  			  MakePoint (500*(ob[5].parametrs[0]+ob[5].parametrs[1]/2)/par.screen, 50+par.a_y*250*ob[5].parametrs[2]/par.screen));
			//рисуем диафрагму
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, MakePoint (500*par.diaph/par.screen, 50-par.a_y*500*par.radius/par.screen), 
										  			  MakePoint (500*par.diaph/par.screen, 50-par.a_y*500*par.diaph_r/par.screen));
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, MakePoint (500*par.diaph/par.screen, 50+par.a_y*500*par.radius/par.screen), 
										  			  MakePoint (500*par.diaph/par.screen, 50+par.a_y*500*par.diaph_r/par.screen));
	}
}
// =============== particles tracking with particles drawing ===================
int run_12(struct particle* b,struct parametrs p, struct object* ob)
{
	FILE *F = 0;
	struct vector E,H,H1,f,v;
	double gamma, b1, b2, b3;
	int i, k, t=0, tik=p.tikalka, t_ob;
	SetSleepPolicy(VAL_SLEEP_NONE);
	num = p.t_number;
	if ( p.cb==1 ) bunch = loadbunch(p);
	if ( p.sv )    F = fopen ("graph.dat", "wb"); //save beam particles each step
	do
	{
		for ( i=0; i<p.t_number; i++)
		{	
			switch( b[i].status )
			{
				case in_cathode:
					gamma = pow(1+pow( b[i].px,2 ) + pow( b[i].py,2 ) + pow( b[i].pz,2 ),0.5);            
					v.z = b[i].pz/gamma; //3 считаем скорости
					b[i].z += 2*v.z; 	 //4 считаем новые координаты
					if (b[i].z>0)                                                 
					{                                                             
						switch ( p.gun_mode )
						{
							case 0:  b[i].status=in_free_space; break;
							case 1:  b[i].status=in_gun;
									 b[i].pz=0;
						}
					}
				continue;
				case in_gun:
					E = E_gun(b[i],p);
					b[i].px -= 2*E.x;                                                            
					b[i].py -= 2*E.y;				                                              
					b[i].pz -= 2*E.z;
					if (b[i].pz<=0) { b[i].status=reflection; num--; continue; } //reflection
					gamma = pow( 1 + pow( b[i].px,2 ) + pow( b[i].py,2 ) + pow( b[i].pz,2 ),0.5);
				    v.x = b[i].px/gamma; //3 считаем скорости
				    v.y = b[i].py/gamma;
				    v.z = b[i].pz/gamma;
			    	b[i].x += 2*v.x;	 //4 считаем новые координаты
				    b[i].y += 2*v.y;
				    b[i].z += 2*v.z;
					if (  b[i].z>=p.gun_length  )  b[i].status = in_free_space;
				continue;
				case on_tube: case on_screen: case on_diaphragm: case reflection: continue;
				default:
					switch (b[i].status) // switch works faster than //E = E_obj(t, b[i], p, ob[b[i].status]);
					{
						case in_scan: E = E_obj(t, b[i], p, ob[5]); break;
						case in_RF:   E = E_obj(t, b[i], p, ob[6]); break;
						case in_bunch:E = E_obj(t, b[i], p, ob[7]); break;
						case in_acc:  E = E_obj(t, b[i], p, ob[9]); break;
						default:      E.x=0; E.y=0; E.z=0;
					}
					if (p.space_charge==1)
					{
						E.x += E_space_charge(i,b,p).x;
						E.y += E_space_charge(i,b,p).y;
						E.z += E_space_charge(i,b,p).z;
					}
					b[i].px -= 2*E.x; //2 приращение импульса для электронов
					b[i].py -= 2*E.y;
					b[i].pz -= 2*E.z;
					if ( b[i].pz <=0 ) { b[i].status = reflection; num--; continue;	} //reflection
					gamma = pow( 1 + pow( b[i].px,2 ) + pow( b[i].py, 2) + pow( b[i].pz, 2),0.5);
				    v.x = b[i].px/gamma; //3 считаем скорости
				    v.y = b[i].py/gamma;                                                     
				    v.z = b[i].pz/gamma;                                                     
				    b[i].x += v.x; //4 считаем новые координаты
				    b[i].y += v.y;
				    b[i].z += v.z;
					switch ( b[i].status ) //switch works faster than //H = H_obj(t, b[i].x, b[i].y, b[i].z, p, ob[b[i].status]);
					{
						case in_quad:  H = H_obj(t, b[i], p, ob[0]); break;
						case in_quad2: H = H_obj(t, b[i], p, ob[1]); break;
						case in_quad3: H = H_obj(t, b[i], p, ob[2]); break;
						case in_quad4: H = H_obj(t, b[i], p, ob[3]); break;
						case in_lens:  H = H_obj(t, b[i], p, ob[4]); break;
						case in_RF:	   H = H_obj(t, b[i], p, ob[6]); break;
						case in_bunch: H = H_obj(t, b[i], p, ob[7]); break;
						case in_corr:  H = H_obj(t, b[i], p, ob[8]); break;
						default: 	   H.x=0; H.y=0; H.z=0;
					}
					if (p.space_charge==1)
					{
						H.x+=H_space_charge(i,b,p).x;
						H.y+=H_space_charge(i,b,p).y;
						H.z+=H_space_charge(i,b,p).z;
					}
					H1.x=H.x/gamma;                                                         
				    H1.y=H.y/gamma;                                                         
				    H1.z=H.z/gamma;                                                         
				    b2=1+H1.x*H1.x+H1.y*H1.y+H1.z*H1.z;                                         
				    b1=2-b2;                                                              
				    b3=2*(v.x*H1.x+v.y*H1.y+v.z*H1.z);
				    f.x=-2*(v.y*H1.z-v.z*H1.y); //для электронов
				    f.y=-2*(v.z*H1.x-v.x*H1.z);
				    f.z=-2*(v.x*H1.y-v.y*H1.x);
				    v.x=(v.x*b1+f.x+H1.x*b3)/b2;
				    v.y=(v.y*b1+f.y+H1.y*b3)/b2;
				    v.z=(v.z*b1+f.z+H1.z*b3)/b2;
				    if (v.z<=0) { b[i].status=reflection; num--; continue; }
				    b[i].x+=v.x; //8 считаем новые координаты
				    b[i].y+=v.y;
				    b[i].z+=v.z;
					b[i].status=where_is_particle(b[i],p,ob);
					if	((b[i].status==on_tube)||(b[i].status==on_screen)||(b[i].status==on_diaphragm)) { num--; continue; }
					if	(b[i].status==on_stop) { num=0; break; }
				    b[i].px=v.x*gamma; //10 считаем импульсы
				    b[i].py=v.y*gamma;
				    b[i].pz=v.z*gamma;
			}
		}
		t+=2;
		tik--;
		if (tik<=0) 
		{
			tik=p.tikalka; //Обновление экрана
			CanvasStartBatchDraw (panelS, PANEL_S_CANVAS_2);
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE);
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_WIDTH, 1);
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_FILL_COLOR, VAL_BLACK);
			// clear path in the tube if set
			if (!p.clear_graph) CanvasDrawRect (panelS, PANEL_S_CANVAS_2, MakeRect(50+1-p.a_y*500*p.radius/p.screen,0,p.a_y*2*500*p.radius/p.screen+1-2,500), VAL_DRAW_INTERIOR);
			//рисуем диафрагму
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, MakePoint (500*p.diaph/p.screen, 50-p.a_y*500*p.radius/ p.screen), 
										  			  MakePoint (500*p.diaph/p.screen, 50-p.a_y*500*p.diaph_r/p.screen));
			CanvasDrawLine (panelS, PANEL_S_CANVAS_2, MakePoint (500*p.diaph/p.screen, 50+p.a_y*500*p.radius/ p.screen), 
										  			  MakePoint (500*p.diaph/p.screen, 50+p.a_y*500*p.diaph_r/p.screen));
			t_ob = (ob[7].parametrs[0]+p.t_length/2)/p.t_v;
			switch(p.graph_mode)
			{
				case 1:
					SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE); //рисуем пластины развёртки
					CanvasDrawLine (panelS, PANEL_S_CANVAS_2, MakePoint (500*(ob[5].parametrs[0]-ob[5].parametrs[1]/2)/p.screen, 50-p.a_y*250*ob[5].parametrs[2]/p.screen), 
												  			  MakePoint (500*(ob[5].parametrs[0]+ob[5].parametrs[1]/2)/p.screen, 50-p.a_y*250*ob[5].parametrs[2]/p.screen));
					CanvasDrawLine (panelS, PANEL_S_CANVAS_2, MakePoint (500*(ob[5].parametrs[0]-ob[5].parametrs[1]/2)/p.screen, 50+p.a_y*250*ob[5].parametrs[2]/p.screen), 
										  			  		  MakePoint (500*(ob[5].parametrs[0]+ob[5].parametrs[1]/2)/p.screen, 50+p.a_y*250*ob[5].parametrs[2]/p.screen));
					if (ob[7].parametrs[6]!=0) //рисуем банч
					{
						SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_MAGENTA);
						SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_FILL_COLOR, VAL_MAGENTA);
						//CanvasDrawOval (panelS, PANEL_S_CANVAS_2, 
						//	MakeRect(50-p.a_y*500*( ob[7].parametrs[2]/p.screen), 500*( ob[7].parametrs[0]-ob[7].parametrs[2]+ob[7].X_1)/p.screen,
						//				p.a_y*500*2*ob[7].parametrs[2]/p.screen,  500*2*ob[7].parametrs[2]/p.screen), VAL_DRAW_FRAME_AND_INTERIOR);
						CanvasDrawOval (panelS, PANEL_S_CANVAS_2, 
							MakeRect(50-p.a_y*500* ((ob[7].L_1*10+ob[7].Z_1*10+t-t_ob)/p.screen), 500* (ob[7].parametrs[0]-ob[7].R_1*10+ob[7].X_1*10)/p.screen, 
								  		p.a_y*500*2* ob[7].L_1*10/p.screen , 		 			  500*2*ob[7].R_1*10/p.screen), VAL_DRAW_FRAME_AND_INTERIOR);
						CanvasDrawOval (panelS, PANEL_S_CANVAS_2, 
							MakeRect(50-p.a_y*500* ((ob[7].L_2*10+ob[7].Z_2*10+t-t_ob)/p.screen), 500* (ob[7].parametrs[0]-ob[7].R_2*10+ob[7].X_2*10)/p.screen, 
								  		p.a_y*500*2* ob[7].L_2*10/p.screen , 		 			  500*2*ob[7].R_2*10/p.screen), VAL_DRAW_FRAME_AND_INTERIOR);
					}
				break;
				case 2:
					SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE); //рисуем пластины развёртки
					SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_FILL_COLOR, VAL_WHITE);
					CanvasDrawRect (panelS, PANEL_S_CANVAS_2, 
							MakeRect (50-p.a_y*250*ob[5].parametrs[2]/p.screen+0, 500*(ob[5].parametrs[0]-ob[5].parametrs[1]/2)/p.screen,
							 			 p.a_y*500*ob[5].parametrs[2]/p.screen+1, 500* ob[5].parametrs[1]/p.screen+1), VAL_DRAW_FRAME_AND_INTERIOR);
					SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_MAGENTA);
					SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_FILL_COLOR, VAL_MAGENTA);
					if ((ob[7].parametrs[6]!=0)&&(abs(t-t_ob)<ob[7].L_1*1)) //рисуем банч
					{
						CanvasDrawOval (panelS, PANEL_S_CANVAS_2, 
							MakeRect(50-p.a_y*500*( ob[7].R_1*1+ob[7].Y_1*1)/p.screen, 500*( ob[7].parametrs[0]-ob[7].R_1*1)/p.screen,
										p.a_y*500*2*ob[7].R_1*1/p.screen,			   500*2*ob[7].R_1*10/p.screen),VAL_DRAW_FRAME_AND_INTERIOR);
					}
				break;
				case 3:
					SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE); //рисуем пластины развёртки
					CanvasDrawLine (panelS, PANEL_S_CANVAS_2, MakePoint (500*(ob[5].parametrs[0]-ob[5].parametrs[1]/2)/p.screen, 50-p.a_y*250*ob[5].parametrs[2]/p.screen), 
												  			  MakePoint (500*(ob[5].parametrs[0]+ob[5].parametrs[1]/2)/p.screen, 50-p.a_y*250*ob[5].parametrs[2]/p.screen));
					CanvasDrawLine (panelS, PANEL_S_CANVAS_2, MakePoint (500*(ob[5].parametrs[0]-ob[5].parametrs[1]/2)/p.screen, 50+p.a_y*250*ob[5].parametrs[2]/p.screen), 
										  			  		  MakePoint (500*(ob[5].parametrs[0]+ob[5].parametrs[1]/2)/p.screen, 50+p.a_y*250*ob[5].parametrs[2]/p.screen));
				break;
				case 4:
					SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_WHITE); //рисуем пластины развёртки
					CanvasDrawLine (panelS, PANEL_S_CANVAS_2, MakePoint (500*(ob[5].parametrs[0]-ob[5].parametrs[1]/2)/p.screen, 50-p.a_y*250*ob[5].parametrs[2]/p.screen), 
												  			  MakePoint (500*(ob[5].parametrs[0]+ob[5].parametrs[1]/2)/p.screen, 50-p.a_y*250*ob[5].parametrs[2]/p.screen));
					CanvasDrawLine (panelS, PANEL_S_CANVAS_2, MakePoint (500*(ob[5].parametrs[0]-ob[5].parametrs[1]/2)/p.screen, 50+p.a_y*250*ob[5].parametrs[2]/p.screen), 
										  			  		  MakePoint (500*(ob[5].parametrs[0]+ob[5].parametrs[1]/2)/p.screen, 50+p.a_y*250*ob[5].parametrs[2]/p.screen));
			}
			SetCtrlAttribute (panelS, PANEL_S_CANVAS_2, ATTR_PEN_COLOR, VAL_GREEN);
			for(k=0;k<p.t_number_draw;k++)
			{
				if((b[k].status!=in_cathode)&&(b[k].status!=on_screen))
				{
					switch(p.graph_mode) //draw beam particles
					{
						case 1: CanvasDrawPoint (panelS, PANEL_S_CANVAS_2, MakePoint(500*b[k].z/p.screen,50-p.a_y*500*b[k].x/p.screen)); break;// x-z plane   
						case 2: CanvasDrawPoint (panelS, PANEL_S_CANVAS_2, MakePoint(500*b[k].z/p.screen,50-p.a_y*500*b[k].y/p.screen)); break;// y-z plane
						case 3: CanvasDrawPoint (panelS, PANEL_S_CANVAS_2, MakePoint(500*b[k].z/p.screen,50-p.a_y*500*( b[k].x+b[k].y)*(sqrt(2)/2)/p.screen)); break;// xy-z plane
						case 4: CanvasDrawPoint (panelS, PANEL_S_CANVAS_2, MakePoint(500*b[k].z/p.screen,50-p.a_y*500*(-b[k].x+b[k].y)*(sqrt(2)/2)/p.screen)); // xy-z plane
					}
				}
			}
			if(p.plot_3D) draw_3D(b,p);
			CanvasEndBatchDraw (panelS, PANEL_S_CANVAS_2);
			if(p.plot_2D)
			{
				CanvasStartBatchDraw (panelS, PANEL_S_CANVAS);
				make_image(beam,par,0);
				make_grid(par);
				CanvasEndBatchDraw (panelS, PANEL_S_CANVAS);
			}
			ProcessDrawEvents ();
			if(p.sv) fwrite (b, sizeof(struct particle), p.t_number, F); //save beam particles for each step in file
			ProcessSystemEvents ();										 //работа кнопки стоп
		}
	}
	while (num>0);
	if(p.sv) fclose (F);
	return t;
}
// =============== particles tracking without particles drawing ==================
int run_0 ( struct particle* b, struct parametrs p, struct object* ob )
{
	double gamma, b1, b2, b3;
	int i, t=0, particle_stop = 0;
	struct vector H, H1, E = {0}, E1, v, f;
	SetSleepPolicy ( VAL_SLEEP_NONE ); //use max processor time
	num = p.t_number;
	if ( p.cb==1 ) bunch = loadbunch(p);
	for ( i=0; i<p.t_number; i++ ) // do run for each particle separately, one by one
	{
		t = 0; particle_stop = 0;
		do
		{
			switch ( b[i].status )
			{
				case in_cathode:
					gamma = pow(1 + pow(b[i].px,2) + pow(b[i].py,2) + pow(b[i].pz,2), 0.5);            
					v.z = b[i].pz/gamma; //3 считаем скорости
					t = (- b[i].z/v.z) + 1;
					b[i].z += v.z*t; 	 //4 считаем новые координаты
					switch ( p.gun_mode )
					{
						case 0: b[i].status = in_free_space; break;
						case 1: b[i].status = in_gun;
								b[i].pz = 0;
					}
				continue;
				case in_gun:
					E = E_gun( b[i], p );
					b[i].px -= 2*E.x;                                                            
					b[i].py -= 2*E.y;				                                              
					b[i].pz -= 2*E.z;
					if ( b[i].pz <= 0 ) { b[i].status = reflection; continue; }
					gamma = pow(1 + pow(b[i].px,2) + pow(b[i].py,2) + pow(b[i].pz,2), 0.5);
					v.x = b[i].px/gamma; //3 считаем скорости
					v.y = b[i].py/gamma;
					v.z = b[i].pz/gamma;
			    	b[i].x += 2*v.x; 	 //4 считаем новые координаты
					b[i].y += 2*v.y;
					b[i].z += 2*v.z;
					if (  b[i].z >= p.gun_length  )  b[i].status = in_free_space;
				continue;
				case on_tube: case on_screen: case reflection: case on_diaphragm: case on_stop:
					particle_stop = 1;
					num --; // particle reached final position or stopped or ...
					if ( num<=0 ) return t; // all particle lost or ...
				break;
				default: //E = ob.E(t, b[i]);  // initial idea. How to realize it?
					switch ( b[i].status ) //switch works faster than //E = E_obj(t, b[i], p, ob[b[i].status]);
					{	
						case in_scan:  E = E_obj(t, b[i], p, ob[5]); break;
						case in_RF:    E = E_obj(t, b[i], p, ob[6]); break;
						case in_bunch: E = E_obj(t, b[i], p, ob[7]); break;
						case in_acc:   E = E_obj(t, b[i], p, ob[9]); break;
						default:	   E.x = 0; E.y = 0; E.z = 0;
					}
					if ( p.space_charge==1 )
					{
						E1 = E_space_charge( i, b, p);
						E.x += E1.x;
						E.y += E1.y;
						E.z += E1.z;
					}
					b[i].px -= 2*E.x; //2 приращение импульса для электронов
					b[i].py -= 2*E.y; 
					b[i].pz -= 2*E.z;
					if (b[i].pz <= 0) { b[i].status = reflection; continue;	}
					gamma = pow(1 + pow(b[i].px,2) + pow(b[i].py,2) + pow(b[i].pz,2), 0.5);
			        v.x = b[i].px/gamma; //3 считаем скорости
			        v.y = b[i].py/gamma;
			        v.z = b[i].pz/gamma;
					b[i].x += v.x; 		 //4 считаем новые координаты
			        b[i].y += v.y;
			        b[i].z += v.z;
					//b[i].status = where_is_particle(b[i],p,ob); //new +24 second for 10000 particles
					switch ( b[i].status ) //switch works faster than //H = H_obj(t, b[i], p, ob[b[i].status]);
					{
						case in_quad:  H = H_obj( t, b[i], p, ob[0]); break;
						case in_quad2: H = H_obj( t, b[i], p, ob[1]); break;
						case in_quad3: H = H_obj( t, b[i], p, ob[2]); break;
						case in_quad4: H = H_obj( t, b[i], p, ob[3]); break;
						case in_lens:  H = H_obj( t, b[i], p, ob[4]); break;
						case in_RF:    H = H_obj( t, b[i], p, ob[6]); break;
						case in_bunch: H = H_obj( t, b[i], p, ob[7]); break;
						case in_corr:  H = H_obj( t, b[i], p, ob[8]); break;
						default:	   H.x = 0; H.y = 0; H.z = 0;
					}
					if ( p.space_charge==1 )
					{
						H1 = H_space_charge( i, b, p);
						H.x += H1.x;
						H.y += H1.y;
						H.z += H1.z;
					}
					H.x = H.x/gamma;                                                         
			        H.y = H.y/gamma;                                                         
			        H.z = H.z/gamma;                                                         
			        b2 = 1 + H.x*H.x + H.y*H.y + H.z*H.z;                                         
			        b1 = 2 - b2;                                                              
			        b3 = 2*(v.x*H.x + v.y*H.y + v.z*H.z);                                           
			        f.x = -2*(v.y*H.z - v.z*H.y); //для электронов
			        f.y = -2*(v.z*H.x - v.x*H.z);
					f.z = -2*(v.x*H.y - v.y*H.x);
			        v.x = (v.x*b1 + f.x + H.x*b3)/b2;
			        v.y = (v.y*b1 + f.y + H.y*b3)/b2;
			        v.z = (v.z*b1 + f.z + H.z*b3)/b2;
			        if ( v.z<=0 ) { b[i].status = reflection; continue; }
			        b[i].x += v.x; //8 считаем новые координаты
			        b[i].y += v.y;
			        b[i].z += v.z;
					b[i].status = where_is_particle( b[i], p, ob); //проверка где находится частица
			        b[i].px = v.x*gamma; //10 считаем импульсы
			        b[i].py = v.y*gamma;
			        b[i].pz = v.z*gamma;
			} // from "switch (b[i].status)"
			t += 2; 
		} // from "do"
		while(!particle_stop);//while (  !((b[i].status==on_tube)||(b[i].status==on_screen)||(b[i].status==reflection)||(b[i].status==on_diaphragm)||(b[i].status==on_stop))  );
		ProcessSystemEvents (); //работа кнопки стоп
	} // from "for (i=0;i<p.t_number;i++)"
	return t; // number of steps*2 for last particle
}

int CVICALLBACK Analyse (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{   // ========== image analysis for SNS: first derivative from curve ========
	int k,l,i; double x_min,x_max, *X,*Y, *x_mass, *y_mass;
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal(panelA,PANEL_A_N,&N);
			GetCtrlVal(panelA,PANEL_A_XMIN,&x_min);
			GetCtrlVal(panelA,PANEL_A_XMAX,&x_max);
			K=malloc(sizeof(double)*N);
			KX=malloc(sizeof(double)*N);
			X=malloc(sizeof(double)*par.t_number);
			Y=malloc(sizeof(double)*par.t_number);
			x_mass=malloc(sizeof(double)*(N+1));
			y_mass=malloc(sizeof(double)*(N+1));
			for (k=0;k<N;k++)
			{
				l=0;
				for (i=0;i<par.t_number;i++)
				{
					if ((beam[i].status==on_screen)&&(beam[i].x*(par.step/fl)>x_min+k*(x_max-x_min)/N)&&(beam[i].x*(par.step/fl)<x_min+(k+1)*(x_max-x_min)/N))
					{
						X[l]=beam[i].x*(par.step/fl);
						Y[l]=beam[i].y*(par.step/fl);
						l++;
					}
				}
				x_mass[k]=0;
				y_mass[k]=0;
				for (i=0;i<l;i++) { x_mass[k]+= X[i]/l; y_mass[k]+= Y[i]/l; }
				SetCtrlAttribute (panelS, PANEL_S_CANVAS, ATTR_PEN_WIDTH, 3);
				SetCtrlAttribute (panelS, PANEL_S_CANVAS, ATTR_PEN_COLOR, 0xff);
				
				CanvasDrawPoint (panelS, PANEL_S_CANVAS, MakePoint((x_mass[k]-par.x_min)*width/(par.x_max-par.x_min),
					height-(y_mass[k]-par.y_min)*height/(par.y_max-par.y_min)));				
			}
			for (k=0;k<N-1;k++)
			{
				if (x_mass[k+1]!=x_mass[k]) K[k]=1-(y_mass[k+1]-y_mass[k])/(x_mass[k+1]-x_mass[k]);
				KX[k]=x_min+(k+1)*(x_max-x_min)/N;
			}
			DeleteGraphPlot (panelA, PANEL_A_GRAPH, plotHandle, VAL_IMMEDIATE_DRAW);
			plotHandle = PlotXY (panelA, PANEL_A_GRAPH, KX, K, N-1, VAL_DOUBLE, VAL_DOUBLE, VAL_FAT_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_WHITE);
			free(KX); free(K); free(X);  free(Y); free(y_mass); free(x_mass);
	}
	return 0;
}

/*   //++++++++++++++++ Written by Malyutin +++++++++++++++++++
int CVICALLBACK Analyse (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{ 
	int k,l,i; double x_min,x_max, *X,*Y, A,B,C,D;
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal(panelA,PANEL_A_N,&N);
			GetCtrlVal(panelA,PANEL_A_XMIN,&x_min);
			GetCtrlVal(panelA,PANEL_A_XMAX,&x_max);
			K=malloc(sizeof(double)*N);
			KX=malloc(sizeof(double)*N);
			X=malloc(sizeof(double)*par.t_number);
			Y=malloc(sizeof(double)*par.t_number);
			for (k=0;k<N;k++)
			{
				l=0;
				for (i=0;i<par.t_number;i++)
				{
					if ((beam[i].status==on_screen)&&(beam[i].x*(par.step/fl)>x_min+k*(x_max-x_min)/N)&&(beam[i].x*(par.step/fl)<x_min+(k+1)*(x_max-x_min)/N))
					{
						X[l]=beam[i].x*(par.step/fl);
						Y[l]=beam[i].y*(par.step/fl);
						l++;
					}
				}
				A=B=C=D=0;
				for (i=0;i<l;i++)
				{
					A+=X[i]*Y[i];
					B+=Y[i];
					C+=X[i];
					D+=X[i]*X[i];
				}
				K[k]=1-(l*A-B*C)/(l*D-C*C);
				KX[k]=x_min+(k+0.5)*(x_max-x_min)/N;
			}
			DeleteGraphPlot (panelA, PANEL_A_GRAPH, plotHandle, VAL_IMMEDIATE_DRAW);
			plotHandle = PlotXY (panelA, PANEL_A_GRAPH, KX, K, N, VAL_DOUBLE, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_GREEN);
			free(KX);free(K);free(X);free(Y);
		break;
	}
	return 0;
}*/

int CVICALLBACK Save_a (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	FILE* f; int i;
	switch (event)
	{
		case EVENT_COMMIT:
			f = fopen ("save\\diff.txt","w");		
			for (i=0;i<N-1;i++) fprintf(f,"%lf %lf\n", KX[i], K[i]);
			fclose (f);
	}
	return 0;
}

int CVICALLBACK Polyfit (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	int i,j,ii=0,num_p_fit=0,order; double *x_for_fit,*y_for_fit,*y_fit,*dy_fit,*c_out,msn,x_min,x_max;
	switch (event)
	{
		case EVENT_COMMIT:
			for(i=0;i<par.t_number;i++) if(beam[i].status==on_screen) num_p_fit++;
			if (num_p_fit==0) return 0; 
			GetCtrlVal(panelA,PANEL_A_ORDER,&order);
			GetCtrlVal(panelA,PANEL_A_XMIN,&x_min);
			GetCtrlVal(panelA,PANEL_A_XMAX,&x_max);
			if(order>num_p_fit) order=num_p_fit;
			x_for_fit=malloc(sizeof(double)*(num_p_fit+1));
			y_for_fit=malloc(sizeof(double)*(num_p_fit+1));
			y_fit=malloc(sizeof(double)*(num_p_fit+1));
			dy_fit=malloc(sizeof(double)*(num_p_fit+1));
			c_out=malloc(sizeof(double)*(order+1));
			for(i=0;i<par.t_number;i++)
				if(beam[i].status==on_screen)
				{
					x_for_fit[ii]=beam[i].x*par.step/fl;					
					y_for_fit[ii]=beam[i].y*par.step/fl;
					ii++;
				}
			PolyFit (x_for_fit, y_for_fit, ii, order, y_fit, c_out, &msn);
			for(i=0;i<ii;i++) //считаем производную
			{
				dy_fit[i]=0;
				if((x_for_fit[i]>x_min)&&(x_for_fit[i]<x_max))
				{
					for(j=0;j<order+1;j++) dy_fit[i]+=c_out[j]*j*pow(x_for_fit[i],j-1);
					dy_fit[i]=1-dy_fit[i];
				}
				else dy_fit[i]=0;
			}
			DeleteGraphPlot (panelA, PANEL_A_GRAPH, -1, VAL_IMMEDIATE_DRAW);
			PlotXY (panelA, PANEL_A_GRAPH, x_for_fit, dy_fit, ii, VAL_DOUBLE, VAL_DOUBLE, VAL_SCATTER, VAL_SOLID_SQUARE,VAL_SOLID, 1, VAL_RED);			
			SetCtrlAttribute (panelS, PANEL_S_CANVAS, ATTR_PEN_WIDTH, 3);
			SetCtrlAttribute (panelS, PANEL_S_CANVAS, ATTR_PEN_COLOR, VAL_RED);
			for(i=0;i<ii;i++) CanvasDrawPoint (panelS, PANEL_S_CANVAS, MakePoint((x_for_fit[i]-par.x_min)*width/(par.x_max-par.x_min), height-(y_fit[i]-par.y_min)*height/(par.y_max-par.y_min)));
			free(x_for_fit); 
			free(y_for_fit); 
			free(y_fit);
			free(dy_fit);
			free(c_out);
	}
	return 0;
}
