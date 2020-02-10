#include "3DGraphCtrl.h"
#include "toolbox.h"
#include <utility.h>
#include <ansi_c.h>
#include <cvirte.h>		
#include <userint.h>
#include "simulator_beam.h"
#include "simulator_vars.h"
#include "simulator_draw.h"
#include "simulator_calc.h"
#include "make_bunch_my.h"
#include "lDSystem.h"

int Get_quadrupole_fields(void) //считываем поля из файлов (поля в файле в кГс)
{
	FILE *F1,*F2;
	int i,j,k; char str[100]; double x,y,z,Bx,By,Bz;
	if ( quadrupole_parametrs.sw==1) return 0;
	F1=fopen("quadrupole1.dat","r"); F2=fopen("quadrupole2.dat","r");
	for (k=0;k<=QNZ;k++)
	{
		for (j=0;j<=QNY;j++)
		{
			for (i=0;i<=QNX;i++)
			{
				fgets(str,sizeof(str),F1);
				sscanf(str,"%lf %lf %lf %lf %lf %lf", &x, &y, &z, &Bx, &By, &Bz);
				quadrupole_fields1[i][j][k].x=Bx*1000*fe*par.step;
				quadrupole_fields1[i][j][k].y=By*1000*fe*par.step;
				quadrupole_fields1[i][j][k].z=Bz*1000*fe*par.step;
				fgets(str,sizeof(str),F2);
				sscanf(str,"%lf %lf %lf %lf %lf %lf", &x, &y, &z, &Bx, &By, &Bz);
				quadrupole_fields2[i][j][k].x=Bx*1000*fe*par.step;
				quadrupole_fields2[i][j][k].y=By*1000*fe*par.step;
				quadrupole_fields2[i][j][k].z=Bz*1000*fe*par.step;
			}
		}
	}
	quadrupole_parametrs.sw=1;
	quadrupole_parametrs.size_x  = 12.5*fl/par.step;  //размер полей 5.0 см
	quadrupole_parametrs.size_y  = 12.5*fl/par.step;  //размер полей 5.0 см
	quadrupole_parametrs.size_z  = 10.0*fl/par.step;  //размер полей 10 см
	quadrupole_parametrs.size_dx = quadrupole_parametrs.size_x/QNX;
	quadrupole_parametrs.size_dy = quadrupole_parametrs.size_y/QNY;
	quadrupole_parametrs.size_dz = quadrupole_parametrs.size_z/QNZ;
	return 1;
}

void Get_parametrs(void)
{
	GetCtrlVal (panelSet, PANEL_SET_RING_T,      	&par.t_distr);
	GetCtrlVal (panelSet, PANEL_SET_SIZE,        	&par.t_radius);
	GetCtrlVal (panelSet, PANEL_SET_LENGTH,      	&par.t_length);
	GetCtrlVal (panelSet, PANEL_SET_TEMP,        	&par.t_temp);
	GetCtrlVal (panelSet, PANEL_SET_NUM,         	&par.t_number);
	GetCtrlVal (panelSet, PANEL_SET_ENERGY_T,    	&par.t_energy);
	GetCtrlVal (panelSet, PANEL_SET_SPREAD,      	&par.t_spread);
	GetCtrlVal (panelSet, PANEL_SET_WEIGHT_T,    	&par.t_weight);
	SetCtrlVal (panelSet, PANEL_SET_CURRENT,     	1.6e-4*par.t_weight*par.t_number/par.t_length);
	GetCtrlVal (panelSet, PANEL_SET_SPACE_CHARGE,	&par.space_charge);
	GetCtrlVal (panelSet, PANEL_SET_LENGTH_G,    	&par.gun_length);
	GetCtrlVal (panelSet, PANEL_SET_RADIUS,      	&par.radius);
	GetCtrlVal (panelSet, PANEL_SET_SCREEN,      	&par.screen);
	GetCtrlVal (panelSet, PANEL_SET_STOP,        	&par.stop);
	GetCtrlVal (panelSet, PANEL_SET_DIAPH,       	&par.diaph);
	GetCtrlVal (panelSet, PANEL_SET_DIAPH_R,     	&par.diaph_r);
	GetCtrlVal (panelSet, PANEL_SET_STEP,        	&par.step);
	GetCtrlVal (panelSet, PANEL_SET_RING_B_S,    	&par.cb);
	GetCtrlVal (panelSet, PANEL_SET_N_FLAT,	    	&par.N_Flat);
	//GetCtrlVal (panelSet, PANEL_SET_RING_SCAN,    &par.scan_mode);
	GetCtrlVal (panelS, PANEL_S_X_MIN,           	&par.x_min);
	GetCtrlVal (panelS, PANEL_S_X_MAX,           	&par.x_max);
	GetCtrlVal (panelS, PANEL_S_Y_MIN,           	&par.y_min);
	GetCtrlVal (panelS, PANEL_S_Y_MAX,           	&par.y_max);
	GetCtrlVal (panelS, PANEL_S_GRAPH_MODE,      	&par.graph_mode);
	GetCtrlVal (panelS, PANEL_S_TIKALKA,         	&par.tikalka);
	GetCtrlVal (panelS, PANEL_S_NP,					&par.t_number_draw);
	GetCtrlVal (panelS, PANEL_S_A_Y,             	&par.a_y);
	GetCtrlVal (panelS, PANEL_S_SV,              	&par.sv);
	GetCtrlVal (panelS, PANEL_S_CLEAR_GRAPH,     	&par.clear_graph);
	GetCtrlVal (panelS, PANEL_S_PLOT_2D,	     	&par.plot_2D);
	GetCtrlVal (panelS, PANEL_S_PLOT_3D,	     	&par.plot_3D);
	GetCtrlVal (panelSet, PANEL_SET_QUAD_Z_1,    	&(obj[0].parametrs[2]));		 //положение первого квадруполя
	//GetCtrlVal (panelSet, PANEL_SET_QUAD_L_1,    &(obj[0].parametrs[3]));
	GetCtrlVal (panelSet, PANEL_SET_QUAD_G_1,    	&(obj[0].parametrs[4]));		 //сила первого квадруполя
	//GetCtrlVal (panelSet, PANEL_SET_RING_QUAD_1, &(obj[0].parametrs[3]));
	GetCtrlVal (panelSet, PANEL_SET_QUAD_Z_2,    	&(obj[0].parametrs[3]));		 //положение второго квадруполя
	//GetCtrlVal (panelSet, PANEL_SET_QUAD_L_2,    &(obj[0].parametrs[5]));
	GetCtrlVal (panelSet, PANEL_SET_QUAD_G_2,    	&(obj[0].parametrs[5]));		 //сила второго квадруполя
	//GetCtrlVal (panelSet, PANEL_SET_RING_QUAD_2, &(obj[1].parametrs[3]));
	GetCtrlVal (panelSet, PANEL_SET_QUAD_PHI_1,     &(obj[0].parametrs[6]));		 //угол системы квадруполей
	obj[0].parametrs[0]=(obj[0].parametrs[3]+obj[0].parametrs[2])/2;			 	 //центр между двумя квадруполями
	obj[0].parametrs[1]= obj[0].parametrs[3]-obj[0].parametrs[2]+20;			 	 //расстояние между квадруполями + 20 см
	GetCtrlVal (panelSet, PANEL_SET_QUAD_X_1,    	&(obj[0].parametrs[7]));
	GetCtrlVal (panelSet, PANEL_SET_QUAD_Y_1,    	&(obj[0].parametrs[8]));
	GetCtrlVal (panelSet, PANEL_SET_QUAD_Z_3,    	&(obj[2].parametrs[0]));
	GetCtrlVal (panelSet, PANEL_SET_QUAD_L_3,    	&(obj[2].parametrs[1]));
	GetCtrlVal (panelSet, PANEL_SET_QUAD_G_3,    	&(obj[2].parametrs[2]));
	GetCtrlVal (panelSet, PANEL_SET_RING_QUAD_3, 	&(obj[2].parametrs[3]));
	GetCtrlVal (panelSet, PANEL_SET_QUAD_PHI_3,  	&(obj[2].parametrs[6]));
	GetCtrlVal (panelSet, PANEL_SET_QUAD_X_1,    	&(obj[2].parametrs[7]));
	GetCtrlVal (panelSet, PANEL_SET_QUAD_Y_1,    	&(obj[2].parametrs[8]));
	GetCtrlVal (panelSet, PANEL_SET_QUAD_Z_4,    	&(obj[3].parametrs[0]));
	GetCtrlVal (panelSet, PANEL_SET_QUAD_L_4,    	&(obj[3].parametrs[1]));
	GetCtrlVal (panelSet, PANEL_SET_QUAD_G_4,    	&(obj[3].parametrs[2]));
	GetCtrlVal (panelSet, PANEL_SET_RING_QUAD_4, 	&(obj[3].parametrs[3]));
	GetCtrlVal (panelSet, PANEL_SET_QUAD_PHI_4,  	&(obj[3].parametrs[6]));
	GetCtrlVal (panelSet, PANEL_SET_QUAD_X_1,    	&(obj[3].parametrs[7]));
	GetCtrlVal (panelSet, PANEL_SET_QUAD_Y_1,    	&(obj[3].parametrs[8]));
	GetCtrlVal (panelSet, PANEL_SET_LENS,        	&(obj[4].parametrs[0]));
	GetCtrlVal (panelSet, PANEL_SET_LENS_L,      	&(obj[4].parametrs[1]));
	GetCtrlVal (panelSet, PANEL_SET_B_MAX,       	&(obj[4].parametrs[2]));
	GetCtrlVal (panelSet, PANEL_SET_SCAN,        	&(obj[5].parametrs[0]));
	GetCtrlVal (panelSet, PANEL_SET_SCAN_L,      	&(obj[5].parametrs[1]));
	GetCtrlVal (panelSet, PANEL_SET_SCAN_W,      	&(obj[5].parametrs[2]));
	GetCtrlVal (panelSet, PANEL_SET_UX_SCAN,     	&(obj[5].parametrs[3]));
	GetCtrlVal (panelSet, PANEL_SET_UY_SCAN,     	&(obj[5].parametrs[4]));
	GetCtrlVal (panelSet, PANEL_SET_TIME_S,      	&(obj[5].parametrs[5]));
	GetCtrlVal (panelSet, PANEL_SET_RING_SCAN,   	&(obj[5].parametrs[6]));
	GetCtrlVal (panelSet, PANEL_SET_RF,          	&(obj[6].parametrs[0]));
	GetCtrlVal (panelSet, PANEL_SET_RF_L,        	&(obj[6].parametrs[1]));									     
	GetCtrlVal (panelSet, PANEL_SET_U_RF,        	&(obj[6].parametrs[2]));
	GetCtrlVal (panelSet, PANEL_SET_PERIOD,      	&(obj[6].parametrs[3]));
	GetCtrlVal (panelSet, PANEL_SET_LAMBDA,      	&(obj[6].parametrs[4]));
	GetCtrlVal (panelSet, PANEL_SET_RING_RF,     	&(obj[6].parametrs[5]));
	GetCtrlVal (panelSet, PANEL_SET_BUNCH,       	&(obj[7].parametrs[0]));
	GetCtrlVal (panelSet, PANEL_SET_BUNCH_RADIUS,	&(obj[7].parametrs[1]));
	GetCtrlVal (panelSet, PANEL_SET_SIGMA_L_1,   	&(obj[7].L_1));
	GetCtrlVal (panelSet, PANEL_SET_SIGMA_L_2,   	&(obj[7].L_2));				obj[7].parametrs[2] = obj[7].L_1 + obj[7].L_2;
	GetCtrlVal (panelSet, PANEL_SET_BUNCH_ENERGY,	&(obj[7].parametrs[5]));
	GetCtrlVal (panelSet, PANEL_SET_BUNCH_N,     	&(obj[7].parametrs[6]));
	GetCtrlVal (panelSet, PANEL_SET_BUNCH_TYPE,  	&(obj[7].parametrs[7]));
	GetCtrlVal (panelSet, PANEL_SET_RING_B_LN_D, 	&(obj[7].parametrs[8]));
	GetCtrlVal (panelSet, PANEL_SET_RING_B_TRD_1,	&(obj[7].TRD_1));
	GetCtrlVal (panelSet, PANEL_SET_RING_B_TRD_2,	&(obj[7].TRD_2));
	GetCtrlVal (panelSet, PANEL_SET_SIGMA_R_1,   	&(obj[7].R_1));
	GetCtrlVal (panelSet, PANEL_SET_SIGMA_R_2,   	&(obj[7].R_2));
	GetCtrlVal (panelSet, PANEL_SET_BUNCH_X_1,   	&(obj[7].X_1));
	GetCtrlVal (panelSet, PANEL_SET_BUNCH_X_2,   	&(obj[7].X_2));
	GetCtrlVal (panelSet, PANEL_SET_BUNCH_Y_1,   	&(obj[7].Y_1));
	GetCtrlVal (panelSet, PANEL_SET_BUNCH_Y_2,   	&(obj[7].Y_2));
	GetCtrlVal (panelSet, PANEL_SET_BUNCH_Z_1,   	&(obj[7].Z_1));
	GetCtrlVal (panelSet, PANEL_SET_BUNCH_Z_2,   	&(obj[7].Z_2));
	GetCtrlVal (panelSet, PANEL_SET_NUMERICSLIDE_1, &(obj[7].K_1));
	GetCtrlVal (panelSet, PANEL_SET_NUMERICSLIDE_2, &(obj[7].K_2));
	GetCtrlVal (panelSet, PANEL_SET_CR,          	&(obj[8].parametrs[0]));
	GetCtrlVal (panelSet, PANEL_SET_CR_L,        	&(obj[8].parametrs[1]));
	GetCtrlVal (panelSet, PANEL_SET_CR_BX,       	&(obj[8].parametrs[2]));
	GetCtrlVal (panelSet, PANEL_SET_CR_BY,       	&(obj[8].parametrs[3]));
	GetCtrlVal (panelSet, PANEL_SET_ACC,         	&(obj[9].parametrs[0]));
	GetCtrlVal (panelSet, PANEL_SET_ACC_L,       	&(obj[9].parametrs[1]));
	GetCtrlVal (panelSet, PANEL_SET_ACC_G,       	&(obj[9].parametrs[2]));
	GetCtrlVal (panelSet, PANEL_SET_ACC_PERIOD,  	&(obj[9].parametrs[3]));
	Get_quadrupole_fields();
}

void save_config(struct parametrs* p, struct object* ob, char *fn)
{
	int i_obj;
	FILE *f, *f2;
	f = fopen (fn,"w");		
	fprintf(f,"\n     Параметры тестового пучка\n\n");
	fprintf(f,"%-10.6lf радиус тестового пучка [см]\n",               p->t_radius);
	fprintf(f,"%-10.1lf длина тестового пучка [пс]\n",                p->t_length);
	fprintf(f,"%-10.3lf температура тестового пучка [эВ]\n",          p->t_temp);
	fprintf(f,"%-10d число частиц в тестовом пучке\n",                p->t_number);
	fprintf(f,"%-10.0lf энергия тестового пучка [кэВ]\n",             p->t_energy);
	fprintf(f,"%-10d вид распределения частиц тестового пучка\n",     p->t_distr);
	fprintf(f,"%-10.1lf длина пушки\n",p->gun_length);
	fprintf(f,"\n     Параметры исследуемого сгустка\n\n");
	fprintf(f,"%-10.0lf тип частиц сгустка\n",                        ob[7].parametrs[7]);
	fprintf(f,"%-10.2e число частиц в сгустках\n",                    ob[7].parametrs[6]);
	fprintf(f,"%-10.1lf процент заряда в первом сгустке [%%]\n",      ob[7].K_1);
	fprintf(f,"%-10.2lf энергия сгустка [ГэВ]\n",                     ob[7].parametrs[5]);
	fprintf(f,"%-10.4lf вертикальное смещение сгустка 1 [см]\n",      ob[7].Y_1);
	fprintf(f,"%-10.4lf горизонтальное смещение сгустка 1 [см]\n",    ob[7].X_1);
	fprintf(f,"%-10.4lf продольное смещение сгустка 1 [см]\n",    	  ob[7].Z_1);
	fprintf(f,"%-10.4lf поперечная сигма сгустка 1 [см]\n",           ob[7].R_1);
	fprintf(f,"%-10.4lf продольная сигма сгустка 1 [см]\n",           ob[7].L_1);
	fprintf(f,"%-10.4lf вертикальное смещение сгустка 2 [см]\n",      ob[7].Y_2);
	fprintf(f,"%-10.4lf горизонтальное смещение сгустка 2 [см]\n",    ob[7].X_2);
	fprintf(f,"%-10.4lf продольное смещение сгустка 2 [см]\n",    	  ob[7].Z_2);
	fprintf(f,"%-10.4lf поперечная сигма сгустка 2 [см]\n",           ob[7].R_2);
	fprintf(f,"%-10.4lf продольная сигма сгустка 2 [см]\n",           ob[7].L_2);
	fprintf(f,"\n     Геометрия системы\n\n");
	fprintf(f,"%-10.1lf радиус трубы [см]\n",                         p->radius);
	fprintf(f,"%-10.1lf расстояние до линзы [см]\n",                  ob[4].parametrs[0]);
	fprintf(f,"%-10.1lf толщина линзы [см]\n",                        ob[4].parametrs[1]);
	fprintf(f,"%-10.1lf расстояние до резонатора [см]\n",             ob[6].parametrs[0]);
	fprintf(f,"%-10.1lf толщина резонатора [см]\n",                   ob[6].parametrs[1]);
	fprintf(f,"%-10.1lf расстояние до пластин развертки [см]\n",      ob[5].parametrs[0]);
	fprintf(f,"%-10.1lf длина пластин развёртки [см]\n",              ob[5].parametrs[1]);
	fprintf(f,"%-10.1lf расстояние между пластин развертки [см]\n",   ob[5].parametrs[2]);
	fprintf(f,"%-10.1lf расстояние до %d квадруполя [см]\n",     	  ob[0].parametrs[2],1);
	fprintf(f,"%-10.1lf толщина %d квадруполя [см]\n",				  0.0,1);
	fprintf(f,"%-10.1lf расстояние до %d квадруполя [см]\n",       	  ob[0].parametrs[3],2);
	fprintf(f,"%-10.1lf толщина %d квадруполя [см]\n",                0.0,2);
	for (i_obj=2;i_obj<4;i_obj++)
	{
		fprintf(f,"%-10.1lf расстояние до %d квадруполя [см]\n",      ob[i_obj].parametrs[0],i_obj+1);
		fprintf(f,"%-10.1lf толщина %d квадруполя [см]\n",            ob[i_obj].parametrs[1],i_obj+1);
	}
	fprintf(f,"%-10.1lf расстояние до сгустка [см]\n",                ob[7].parametrs[0]);
	fprintf(f,"%-10.1lf ширина электромагнитных полей сгутка [см]\n", ob[7].parametrs[1]);
	fprintf(f,"%-10.1lf расстояние до экрана [см]\n",                 p->screen);
	fprintf(f,"%-10.1lf расстояние до диафрагмы [см]\n",              p->diaph);
	fprintf(f,"%-10.3lf радиус диафрагмы [см]\n",                     p->diaph_r);
	fprintf(f,"\n     Параметры\n\n");	
	fprintf(f,"%-10.3lf шаг по времени [пс]\n",                       p->step);
	fprintf(f,"%-10.1lf максимальное напряжение в резонаторе [В]\n",  ob[6].parametrs[2]);
	fprintf(f,"%-10.1lf период ВЧ в резонаторе [нс]\n",               ob[6].parametrs[3]);
	fprintf(f,"%-10.3lf максимальное магнитное поле в линзе [кГс]\n", ob[4].parametrs[2]);
	fprintf(f,"%-10.4lf Градиент в %d квадруполе [кГс/см]\n",         ob[0].parametrs[4],1);
	fprintf(f,"%-10.4lf Градиент в %d квадруполе [кГс/см]\n",         ob[0].parametrs[5],2);
	for (i_obj=2;i_obj<4;i_obj++)
	fprintf(f,"%-10.4lf Градиент в %d квадруполе [кГс/см]\n",    	  ob[i_obj].parametrs[2],i_obj);
	fprintf(f,"%-10.1lf максимальное напряжение на пластинах развертки X [В]\n",ob[5].parametrs[3]);
	fprintf(f,"%-10.1lf максимальное напряжение на пластинах развертки Y [В]\n",ob[5].parametrs[4]);
	fprintf(f,"%-10.1lf время развертки [нс]\n",ob[5].parametrs[5]);
	fprintf(f,"\n     Окно\n\n");	
	fprintf(f,"%-10.2lf x min [см]\n",p->x_min);
	fprintf(f,"%-10.2lf x max [см]\n",p->x_max);
	fprintf(f,"%-10.2lf y min [см]\n",p->y_min);
	fprintf(f,"%-10.2lf y max [см]\n",p->y_max);
	fprintf(f,"%-10.1lf расстояние до корректора [см]\n",             ob[8].parametrs[0]);
	fprintf(f,"%-10.1lf длина корректора [см]\n",                     ob[8].parametrs[1]);
	fprintf(f,"%-10.2lf поле в корректоре Bx [Гс]\n", 				  ob[8].parametrs[2]);
	fprintf(f,"%-10.2lf поле в корректоре By [Гс]\n", 				  ob[8].parametrs[3]);
	fclose (f);
	f2 = fopen ("config.dat","w");
	fwrite (p,sizeof(*p),1,f2);
	fwrite (ob,number_of_object*sizeof(*ob),1,f2);
	fclose (f2);
}

void load_config(struct parametrs* p, char *file_name) //FILE* stream)
{
	int i_obj; char str[70];
	FILE *f2, *stream;
	stream = fopen (file_name,"r");
	if(stream)
	{
		f2 = fopen ("config.dat","r");
		if (f2!=NULL)
		{
			fread (p,sizeof(*p),1,f2);
			fread (obj,number_of_object*sizeof(*obj),1,f2);
			fclose (f2);
		}
		fgets(str,sizeof(str),stream);
		fgets(str,sizeof(str),stream);
		fgets(str,sizeof(str),stream);
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(p->t_radius));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(p->t_length));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(p->t_temp));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%i", 								&(p->t_number));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(p->t_energy));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%d", 								&(p->t_distr));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(p->gun_length));
		fgets(str,sizeof(str),stream);
		fgets(str,sizeof(str),stream);
		fgets(str,sizeof(str),stream);
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[7].parametrs[7]));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[7].parametrs[6]));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[7].K_1)); obj[7].K_2 = 100-obj[7].K_1;
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[7].parametrs[5]));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[7].Y_1));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[7].X_1));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[7].Z_1));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[7].R_1));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[7].L_1));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[7].Y_2));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[7].X_2));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[7].Z_2));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[7].R_2));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[7].L_2));
		fgets(str,sizeof(str),stream);
		fgets(str,sizeof(str),stream);
		fgets(str,sizeof(str),stream);
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(p->radius));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[4].parametrs[0]));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[4].parametrs[1]));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[6].parametrs[0]));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[6].parametrs[1]));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[5].parametrs[0]));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[5].parametrs[1]));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[5].parametrs[2]));
		for (i_obj=0;i_obj<4;i_obj++)
		{
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[i_obj].parametrs[0]));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[i_obj].parametrs[1]));
		}
		obj[0].parametrs[2]=obj[0].parametrs[0];
		obj[0].parametrs[3]=obj[1].parametrs[0];
		obj[0].parametrs[0]=(obj[0].parametrs[3]+obj[0].parametrs[2])/2;
		obj[0].parametrs[1]= obj[0].parametrs[3]-obj[0].parametrs[2]+20;
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[7].parametrs[0]));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[7].parametrs[1]));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(p->screen));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(p->diaph));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(p->diaph_r));
		fgets(str,sizeof(str),stream);
		fgets(str,sizeof(str),stream);
		fgets(str,sizeof(str),stream);
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(p->step));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[6].parametrs[2]));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[6].parametrs[3]));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[4].parametrs[2]));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[0].parametrs[4]));		
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[0].parametrs[5]));
		for (i_obj=2;i_obj<4;i_obj++)
		{
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[i_obj].parametrs[2]));
		}
		for (i_obj=1;i_obj<4;i_obj++)
		{
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[5].parametrs[i_obj+2]));
		}		
		fgets(str,sizeof(str),stream);
		fgets(str,sizeof(str),stream);
		fgets(str,sizeof(str),stream);
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(p->x_min));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(p->x_max));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(p->y_min));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(p->y_max));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[8].parametrs[0]));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[8].parametrs[1]));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[8].parametrs[2]));
		fgets(str,sizeof(str),stream);
		sscanf(str,"%lf",								&(obj[8].parametrs[3]));
	}
	fclose (stream);
}

void set_config(struct parametrs p, struct object *obj)
{
	SetCtrlVal (panelSet, PANEL_SET_RING_T,   	 par.t_distr);
	SetCtrlVal (panelSet, PANEL_SET_SIZE,        par.t_radius);
	SetCtrlVal (panelSet, PANEL_SET_LENGTH,      par.t_length);
	SetCtrlVal (panelSet, PANEL_SET_TEMP,        par.t_temp);
	SetCtrlVal (panelSet, PANEL_SET_NUM,         par.t_number);
	SetCtrlVal (panelSet, PANEL_SET_ENERGY_T,    par.t_energy);
	SetCtrlVal (panelSet, PANEL_SET_LENGTH_G,    par.gun_length);
	SetCtrlVal (panelSet, PANEL_SET_RADIUS,      par.radius);
	SetCtrlVal (panelSet, PANEL_SET_SCREEN,      par.screen);
	SetCtrlVal (panelSet, PANEL_SET_DIAPH,       par.diaph);
	SetCtrlVal (panelSet, PANEL_SET_DIAPH_R,     par.diaph_r);
	SetCtrlVal (panelSet, PANEL_SET_STEP,        par.step);
	SetCtrlVal (panelS,   PANEL_S_X_MIN,         par.x_min);
	SetCtrlVal (panelS,   PANEL_S_X_MAX,         par.x_max);
	SetCtrlVal (panelS,   PANEL_S_Y_MIN,      	 par.y_min);
	SetCtrlVal (panelS,   PANEL_S_Y_MAX,      	 par.y_max);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_Z_1,    obj[0].parametrs[2]);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_L_1,    0.0);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_G_1,    obj[0].parametrs[4]);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_Z_2,    obj[0].parametrs[3]);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_L_2,    0.0);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_G_2,    obj[0].parametrs[5]);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_PHI_3,  obj[0].parametrs[6]);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_X_1,    obj[0].parametrs[7]);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_Y_1,    obj[0].parametrs[8]);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_Z_3,    obj[2].parametrs[0]);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_L_3,    obj[2].parametrs[1]);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_G_3,    obj[2].parametrs[2]);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_PHI_3,  obj[0].parametrs[6]);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_X_1,    obj[2].parametrs[7]);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_Y_1,    obj[2].parametrs[8]);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_Z_4,    obj[3].parametrs[0]);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_L_4,    obj[3].parametrs[1]);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_G_4,    obj[3].parametrs[2]);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_PHI_3,  obj[0].parametrs[6]);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_X_1,    obj[3].parametrs[7]);
	SetCtrlVal (panelSet, PANEL_SET_QUAD_Y_1,    obj[3].parametrs[8]);
	SetCtrlVal (panelSet, PANEL_SET_LENS,        obj[4].parametrs[0]);
	SetCtrlVal (panelSet, PANEL_SET_LENS_L,      obj[4].parametrs[1]);
	SetCtrlVal (panelSet, PANEL_SET_B_MAX,       obj[4].parametrs[2]);
	SetCtrlVal (panelSet, PANEL_SET_SCAN,        obj[5].parametrs[0]);
	SetCtrlVal (panelSet, PANEL_SET_SCAN_L,      obj[5].parametrs[1]);
	SetCtrlVal (panelSet, PANEL_SET_SCAN_W,      obj[5].parametrs[2]);
	SetCtrlVal (panelSet, PANEL_SET_UX_SCAN,     obj[5].parametrs[3]);
	SetCtrlVal (panelSet, PANEL_SET_UY_SCAN,     obj[5].parametrs[4]);
	SetCtrlVal (panelSet, PANEL_SET_TIME_S,      obj[5].parametrs[5]);
	SetCtrlVal (panelSet, PANEL_SET_RF,          obj[6].parametrs[0]);
	SetCtrlVal (panelSet, PANEL_SET_RF_L,        obj[6].parametrs[1]);									     
	SetCtrlVal (panelSet, PANEL_SET_U_RF,        obj[6].parametrs[2]);
	SetCtrlVal (panelSet, PANEL_SET_PERIOD,      obj[6].parametrs[3]);
	SetCtrlVal (panelSet, PANEL_SET_LAMBDA,      obj[6].parametrs[4]);
	SetCtrlVal (panelSet, PANEL_SET_BUNCH,       obj[7].parametrs[0]);
	SetCtrlVal (panelSet, PANEL_SET_BUNCH_RADIUS,obj[7].parametrs[1]);
	SetCtrlVal (panelSet, PANEL_SET_SIGMA_L_1,   obj[7].L_1);
	SetCtrlVal (panelSet, PANEL_SET_SIGMA_L_2,   obj[7].L_2);
	SetCtrlVal (panelSet, PANEL_SET_BUNCH_ENERGY,obj[7].parametrs[5]);
	SetCtrlVal (panelSet, PANEL_SET_BUNCH_N,     obj[7].parametrs[6]);
	SetCtrlVal (panelSet, PANEL_SET_BUNCH_N_2,   obj[7].parametrs[6]*1.6e-7);
	SetCtrlVal (panelSet, PANEL_SET_BUNCH_TYPE,  obj[7].parametrs[7]);
	SetCtrlVal (panelSet, PANEL_SET_RING_B_LN_D, obj[7].parametrs[8]);
	SetCtrlVal (panelSet, PANEL_SET_RING_B_TRD_1,obj[7].TRD_1);
	SetCtrlVal (panelSet, PANEL_SET_RING_B_TRD_1,obj[7].TRD_2);
	SetCtrlVal (panelSet, PANEL_SET_NUMERICSLIDE_1,obj[7].K_1);
	SetCtrlVal (panelSet, PANEL_SET_NUMERICSLIDE_2,obj[7].K_2);
	SetCtrlVal (panelSet, PANEL_SET_SIGMA_R_1,   obj[7].R_1);
	SetCtrlVal (panelSet, PANEL_SET_SIGMA_R_2,   obj[7].R_2);
	SetCtrlVal (panelSet, PANEL_SET_BUNCH_X_1,   obj[7].X_1);
	SetCtrlVal (panelSet, PANEL_SET_BUNCH_X_2,   obj[7].X_2);
	SetCtrlVal (panelSet, PANEL_SET_BUNCH_Y_1,   obj[7].Y_1);
	SetCtrlVal (panelSet, PANEL_SET_BUNCH_Y_2,   obj[7].Y_2);
	SetCtrlVal (panelSet, PANEL_SET_BUNCH_Z_1,   obj[7].Z_1);
	SetCtrlVal (panelSet, PANEL_SET_BUNCH_Z_2,   obj[7].Z_2);
	SetCtrlVal (panelSet, PANEL_SET_CR,          obj[8].parametrs[0]);
	SetCtrlVal (panelSet, PANEL_SET_CR_L,        obj[8].parametrs[1]);
	SetCtrlVal (panelSet, PANEL_SET_CR_BX,       obj[8].parametrs[2]);
	SetCtrlVal (panelSet, PANEL_SET_CR_BY,       obj[8].parametrs[3]);
	/*SetCtrlVal (panelSet, PANEL_SET_ACC,         obj[9].parametrs[0]); % position
	SetCtrlVal (panelSet, PANEL_SET_ACC_L,       obj[9].parametrs[1]); % length
	SetCtrlVal (panelSet, PANEL_SET_ACC_G,       obj[9].parametrs[2]); % gradient
	SetCtrlVal (panelSet, PANEL_SET_ACC_PERIOD,  obj[9].parametrs[3]); % period RF*/
}

void save_info(void)
{
	int bitmapID, year,month,day,hours,minutes,second; char file_name[50];
	GetCtrlBitmap (panelS, PANEL_S_CANVAS, 0, &bitmapID);
	GetSystemDate (&month, &day, &year);
	GetSystemTime (&hours, &minutes, &second);
	sprintf (file_name, "save\\%d_%02d_%02d_%02d_%02d_%02d.png",year,month,day,hours,minutes,second);
	SaveBitmapToPNGFile (bitmapID, file_name);
	sprintf (file_name, "save\\%d_%02d_%02d_%02d_%02d_%02d.txt",year,month,day,hours,minutes,second);
	save_config(&par,obj,file_name);
	if (par.cb==2)
	{
		sprintf (file_name, "save\\%d_%02d_%02d_%02d_%02d_%02d.mkb",year,month,day,hours,minutes,second);
		SavePS (panelM_B, file_name);
	}
}

void save_data(void)
{
	FILE *F; int i;
	F=fopen("save\\1.dat","w");
	for (i=0;i<par.t_number;i++) if (beam[i].status==on_screen) fprintf(F,"%lf %lf\n", beam[i].x*par.step/fl,beam[i].y*par.step/fl);
	fclose (F);
}

int main (int argc, char *argv[])
{
	int i_obj,j; char file_name[MAX_PATHNAME_LEN];
	SetSleepPolicy(VAL_SLEEP_MORE);   	// DisableBreakOnLibraryErrors();
	if (InitCVIRTE (0, argv, 0) == 0)  return -1;	/* out of memory */
	if ((panelS      = LoadPanel (0, "simulator_beam.uir", PANEL_S   )) < 0) return -1;
	if ((panelSet    = LoadPanel (0, "simulator_beam.uir", PANEL_SET )) < 0) return -1;
	if ((panelIM_SET = LoadPanel (0, "simulator_beam.uir", P_IM_SET  )) < 0) return -1;
	if ((panelM_B    = LoadPanel (0, "simulator_beam.uir", PANEL_M_B )) < 0) return -1;
	if ((panelCHNG   = LoadPanel (0, "simulator_beam.uir", PANEL_CHNG)) < 0) return -1;
	if ((panelA      = LoadPanel (0, "simulator_beam.uir", PANEL_A   )) < 0) return -1;
	if ((panelG      = LoadPanel (0, "simulator_beam.uir", PANEL_G   )) < 0) return -1;
	if ((panelPHSP   = LoadPanel (0, "simulator_beam.uir", PANEL_PHSP)) < 0) return -1;
	if ((panel_3D    = LoadPanel (0, "simulator_beam.uir", PANEL_3D  )) < 0) return -1;
	change_switch=0;
	num=0;
	obj=(struct object*)malloc(number_of_object*sizeof(struct object));
	obj[ 0].name=in_quad;  	 // magnetic
	obj[ 1].name=in_quad2; 	 // magnetic
	obj[ 2].name=in_quad3; 	 // magnetic
	obj[ 3].name=in_quad4; 	 // magnetic
	obj[ 4].name=in_lens;  	 // magnetic
	obj[ 5].name=in_scan;  	 // electric
	obj[ 6].name=in_RF;	   	 // electric magnetic
	obj[ 7].name=in_bunch; 	 // electric magnetic
	obj[ 8].name=in_corr;  	 // magnetic
	obj[ 9].name=in_acc;  	 // electric magnetic
	obj[10].name=in_cathode;
	obj[11].name=in_gun;	 // electric magnetic
	obj[12].name=in_free_space;
	obj[13].name=on_tube;
	obj[14].name=on_screen;
	obj[15].name=on_diaphragm;
	obj[16].name=reflection;
	RecallPanelState (panelIM_SET, "simulator_draw.cfg", 0);
	sprintf (file_name, "config.txt");
	load_config(&par,file_name); 	//f=fopen ("config.txt", "r");
	set_config(par,obj);
	for (i_obj=0;i_obj<number_of_object;i_obj++)
	{
		for (j=0;j<9;j++) obj[i_obj].parametrs[j]=0;
		// obj[i_obj].E(int, struct particle)=E_obj(,,par,obj[i_obj]);   // To be realized
		// obj[i_obj].H(,)=H_obj(,,par,obj[i_obj]);						 // To be realized
	}//Get_quadrupole_fields();
	quadrupole_parametrs.sw=0;
	DisplayPanel (panelS);
	RunUserInterface ();
	DiscardPanel (panelS		);
	DiscardPanel (panelSet		);
	DiscardPanel (panelIM_SET	);
	DiscardPanel (panelM_B		);
	DiscardPanel (panelCHNG		);
	DiscardPanel (panelA		);
	DiscardPanel (panelG		);
	DiscardPanel (panelPHSP		);
	DiscardPanel (panel_3D		);
	free(obj);
	return 0;
}

int CVICALLBACK Start (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	int i,t=0, st[9]={0}; char str[50],str2[50];
	struct system_time
	{
		int h,m,s;
	} time_start,time_finish;
	struct parametrs par_in;
	struct object*   obj_in;
	switch (event)
	{
		case EVENT_COMMIT:
			
			GetCtrlAttribute (panelS, PANEL_S_START, ATTR_LABEL_TEXT, str);
			if (!strcmp(str,"Stop"))
			{
				SetCtrlAttribute (panelS, PANEL_S_START, ATTR_LABEL_TEXT, "Start");
				SetCtrlAttribute (panelS, PANEL_S_START, ATTR_CMD_BUTTON_COLOR, VAL_GREEN);
				SetSleepPolicy(VAL_SLEEP_MORE);
				num=0; // global variable = number of flying particles
				return 0;
			}
			//CanvasClear (panelS, PANEL_S_CANVAS_2, VAL_ENTIRE_OBJECT);
			SetCtrlAttribute (panelS, PANEL_S_START, ATTR_LABEL_TEXT, "Stop");
			SetCtrlAttribute (panelS, PANEL_S_START, ATTR_CMD_BUTTON_COLOR, VAL_RED);
			Get_parametrs();
			ProcessDrawEvents ();
			save_config(&par,obj,"config.txt");
			switch((int)par.scan_mode)
			{
				case 1: DSReadAnsysData();
			}
			if (change_switch!=0) //меняем значение параметра программно, change one parameter in loop
			{
				change_number++;
				do
				{
					GetSystemTime (&time_start.h, &time_start.m, &time_start.s);
					change_number--;
					SetCtrlVal (change_panel, change_ctrl, change_ini);
					if (change_panel==panelM_B) GetCtrlVal (change_panel, change_ctrl, &matrix_make_bunch[(change_ctrl-2)/10][(change_ctrl-2)%10][yz]);
					else Get_parametrs();
					obj_in = obj_to_in(par,obj);
					par_in = par_to_in(par,obj);
					beam = (struct particle*) malloc (par.t_number*sizeof(struct particle));
					loadbeam(beam,par_in);
					t=run_0(beam,par_in,obj_in);
					make_image(beam,par_in,1); //рисуем картинку     //
					save_info(); 			   //сохраняем результаты//
					GetSystemTime (&time_finish.h, &time_finish.m, &time_finish.s);
					t=(time_finish.h*3600+time_finish.m*60+time_finish.s)-(time_start.h*3600+time_start.m*60+time_start.s);
					t*=(change_number);
					sprintf(str,"Remaining time -> %d h : %d m : %d s \nCurrent value  -> %lf",
						t/3600, (t-3600*(t/3600))/60, t-3600*(t/3600)-60*((t-3600*(t/3600))/60), change_ini);
					ResetTextBox (panelS, PANEL_S_TEXTBOX, "");
					InsertTextBoxLine (panelS, PANEL_S_TEXTBOX, 0, str);
					ProcessDrawEvents ();
					change_ini+=change_step;
				}
				while (change_number!=0);
				SetCtrlAttribute (change_panel, change_ctrl, ATTR_TEXT_BGCOLOR, VAL_WHITE);
				change_switch=0;
			}
			else
			{
				GetSystemTime (&time_start.h, &time_start.m, &time_start.s);
				//ResetTextBox (panelS, PANEL_S_TEXTBOX, "");
				sprintf(str2,"Started at %02d:%02d:%02d ...",time_start.h, time_start.m, time_start.s);
				InsertTextBoxLine (panelS, PANEL_S_TEXTBOX, 0, str2);
				InsertTextBoxLine (panelS, PANEL_S_TEXTBOX, 1, "");
				//считаем параметры во внутренних единицах
				obj_in=obj_to_in(par,obj);
				par_in=par_to_in(par,obj);
				//задаём начальное распределение частиц тестового пучка
				beam = (struct particle*) malloc (par.t_number*sizeof(struct particle));
				loadbeam(beam,par_in);
				//считаем движение тестового пучка
				switch (par.graph_mode)
				{
					case 0: t=run_0(beam,par_in,obj_in); break;
					case 1: case 2: case 3: case 4: 	 t=run_12(beam,par_in,obj_in);
				}
				GetCtrlAttribute (panelS, PANEL_S_START, ATTR_LABEL_TEXT, str); //check button status
				if(!strcmp(str,"Stop")) //рисуем картинку
				{
					make_image(beam,par,1);
					make_grid(par);
				}
				for (i=0;i<par.t_number;i++) //смотрим где частицы
				{
					switch (beam[i].status)
					{
						case on_tube:      st[5]++; break;
						case on_screen:    st[6]++; break;
						case reflection:   st[7]++; break;
						case on_diaphragm: st[8]++;
					}
				} //ResetTextBox (panelS, PANEL_S_TEXTBOX, "");
				sprintf(str,"On tube -> %d",st[5]);         InsertTextBoxLine (panelS, PANEL_S_TEXTBOX, 1, str);
				sprintf(str,"On diaphragm -> %d",st[8]);	InsertTextBoxLine (panelS, PANEL_S_TEXTBOX, 2, str);
				sprintf(str,"On screen -> %d",st[6]);		InsertTextBoxLine (panelS, PANEL_S_TEXTBOX, 3, str);
				sprintf(str,"Reflection -> %d",st[7]);		InsertTextBoxLine (panelS, PANEL_S_TEXTBOX, 4, str);
				sprintf(str,"Number of step -> %d",t);		InsertTextBoxLine (panelS, PANEL_S_TEXTBOX, 5, str);
				GetSystemTime (&time_finish.h, &time_finish.m, &time_finish.s);
				t=(time_finish.h*3600+time_finish.m*60+time_finish.s)-(time_start.h*3600+time_start.m*60+time_start.s);
				sprintf(str,"All time -> %d h : %d m : %d s",t/3600, (t-3600*(t/3600))/60, t-3600*(t/3600)-60*((t-3600*(t/3600))/60));
				InsertTextBoxLine (panelS, PANEL_S_TEXTBOX, 6, str);
				sprintf(str,"Finished at %02d:%02d:%02d",time_finish.h, time_finish.m, time_finish.s);
				InsertTextBoxLine (panelS, PANEL_S_TEXTBOX, 7, str);
			}
			SetCtrlAttribute (panelS, PANEL_S_START, ATTR_LABEL_TEXT, "Start");
			SetCtrlAttribute (panelS, PANEL_S_START, ATTR_CMD_BUTTON_COLOR, VAL_GREEN);
			SetSleepPolicy(VAL_SLEEP_MORE);
			free(obj_in);
	}
	return 0;
}

int CVICALLBACK Redraw (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal (panelS, PANEL_S_X_MIN, &par.x_min);
			GetCtrlVal (panelS, PANEL_S_X_MAX, &par.x_max);
			GetCtrlVal (panelS, PANEL_S_Y_MIN, &par.y_min);
			GetCtrlVal (panelS, PANEL_S_Y_MAX, &par.y_max);
			make_image(beam,par,1);
			make_grid(par);
	}
	return 0;
}

int CVICALLBACK QuitCallback (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			SavePanelState (panelIM_SET, "simulator_draw.cfg", 0);
			QuitUserInterface (0);
	}
	return 0;
}

int CVICALLBACK Canvas_callback (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
  	double dx,dy;
	dx = par.x_max-par.x_min;
	dy = par.y_max-par.y_min;
	switch (event)
	{
		case EVENT_RIGHT_CLICK:	DisplayPanel (panelIM_SET); 
		break;
		case EVENT_LEFT_CLICK:
			par.x_min+=((double)eventData2/width-0.5)*dx;
			par.x_max+=((double)eventData2/width-0.5)*dx;
			par.y_min-=((double)eventData1/height-0.5)*dy;
			par.y_max-=((double)eventData1/height-0.5)*dy;
			SetCtrlVal (panelS, PANEL_S_X_MIN, par.x_min);
			SetCtrlVal (panelS, PANEL_S_X_MAX, par.x_max);
			SetCtrlVal (panelS, PANEL_S_Y_MIN, par.y_min);
			SetCtrlVal (panelS, PANEL_S_Y_MAX, par.y_max);
			make_image(beam,par,1);
			make_grid(par);
		break;
		case EVENT_LEFT_DOUBLE_CLICK:
			par.x_min=-dx/2;
			par.x_max= dx/2;
			par.y_min=-dy/2;
			par.y_max= dy/2;
			SetCtrlVal (panelS, PANEL_S_X_MIN, par.x_min);
			SetCtrlVal (panelS, PANEL_S_X_MAX, par.x_max);
			SetCtrlVal (panelS, PANEL_S_Y_MIN, par.y_min);
			SetCtrlVal (panelS, PANEL_S_Y_MAX, par.y_max);
			make_image(beam,par,1);
			make_grid(par);
		break;
		case EVENT_MOUSE_WHEEL_SCROLL:
			GetCtrlVal (panelS, PANEL_S_X_MIN, &par.x_min);
			GetCtrlVal (panelS, PANEL_S_X_MAX, &par.x_max);
			GetCtrlVal (panelS, PANEL_S_Y_MIN, &par.y_min);
			GetCtrlVal (panelS, PANEL_S_Y_MAX, &par.y_max);
			if (eventData1==MOUSE_WHEEL_SCROLL_UP)
			{
				par.x_min-=dx/2;
				par.x_max+=dx/2;
				par.y_min-=dy/2;
				par.y_max+=dy/2;
			}
			if (eventData1==MOUSE_WHEEL_SCROLL_DOWN)
			{
				par.x_min+=dx/4;
				par.x_max-=dx/4;
				par.y_min+=dy/4;
				par.y_max-=dy/4;
			}
			SetCtrlVal (panelS, PANEL_S_X_MIN, par.x_min);
			SetCtrlVal (panelS, PANEL_S_X_MAX, par.x_max);
			SetCtrlVal (panelS, PANEL_S_Y_MIN, par.y_min);
			SetCtrlVal (panelS, PANEL_S_Y_MAX, par.y_max);
			make_image(beam,par,1);
			make_grid(par);	
	}
	return 0;
}

int CVICALLBACK Settings (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event) { case EVENT_COMMIT: DisplayPanel (panelSet); }
	return 0;
}
int CVICALLBACK OkCallback (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event) { case EVENT_COMMIT: HidePanel (panelSet); }
	return 0;
}

int CVICALLBACK Graph_mode (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			CanvasStartBatchDraw (panelS, PANEL_S_CANVAS_2);
			CanvasClear (panelS, PANEL_S_CANVAS_2, VAL_ENTIRE_OBJECT);
			GetCtrlVal  (panelS, PANEL_S_GRAPH_MODE, &par.graph_mode);
			GetCtrlVal  (panelS, PANEL_S_A_Y, &par.a_y);
			draw_beam(obj);
			CanvasEndBatchDraw (panelS, PANEL_S_CANVAS_2);
			ProcessDrawEvents ();
	}		
	return 0;
}

int CVICALLBACK Refresh (int panel, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_GOT_FOCUS:
			if (num==0)
			{
				CanvasStartBatchDraw (panelS, PANEL_S_CANVAS_2);
				CanvasClear (panelS, PANEL_S_CANVAS_2, VAL_ENTIRE_OBJECT);
				Get_parametrs();
				draw_beam(obj);
				CanvasEndBatchDraw (panelS, PANEL_S_CANVAS_2);
				ProcessDrawEvents ();
			}
		break;
		case EVENT_LOST_FOCUS: break;
		case EVENT_CLOSE:	   break;
	}
	return 0;
}

int CVICALLBACK Save (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			save_info();
			save_data();
	}
	return 0;
}

int  CVICALLBACK Panel_set(int panel, int event, void *callbackData, int eventData1, int eventData2)
{
	int ctrl,active_panel, data_type;
	double ini;
	switch (event)
	{
		case EVENT_LEFT_DOUBLE_CLICK:
			active_panel = GetActivePanel ();
			ctrl = GetActiveCtrl (active_panel);
			GetCtrlAttribute (active_panel, ctrl, ATTR_DATA_TYPE, &data_type);
			if (data_type==VAL_INTEGER) return 0;
			if ((change_switch==1)&&(ctrl==change_ctrl))
			{
				SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_INI,change_ini);
				SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_FIN,change_fin);
				SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_STEP,change_step);
				SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_NUMBER,change_number);
			}
			else
			{
				GetCtrlVal(active_panel,ctrl,&ini);
				SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_INI,ini);
				SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_FIN,ini);
				SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_STEP,0.0);
				SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_NUMBER,0);
			}
			DisplayPanel (panelCHNG);
	}
	return 0;
}

int CVICALLBACK OkCallbackCHNG (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	int ctrl,active_panel, number;
	double ini,fin,step;
	switch (event)
	{
		case EVENT_COMMIT:
			HidePanel(panelCHNG);
			GetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_INI,&ini);
			GetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_FIN,&fin);
			GetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_STEP,&step);
			GetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_NUMBER,&number);
			active_panel = GetActivePanel ();
			ctrl = GetActiveCtrl (active_panel);
			SetCtrlVal(active_panel,ctrl,ini);
			if (change_switch!=0)
			{
				if ((ini!=fin)&&(step!=0))
				{
					if (change_panel==active_panel)
					{
						if (change_ctrl!=ctrl)
						{
							SetCtrlAttribute (active_panel, change_ctrl, ATTR_TEXT_BGCOLOR, VAL_WHITE);
							SetCtrlAttribute (active_panel, ctrl, ATTR_TEXT_BGCOLOR, VAL_RED);
						}
					}
					else
					{
						SetCtrlAttribute (change_panel, change_ctrl, ATTR_TEXT_BGCOLOR, VAL_WHITE);
						SetCtrlAttribute (active_panel, ctrl, ATTR_TEXT_BGCOLOR, VAL_RED);
						change_panel=active_panel;
					}
					change_ctrl=ctrl;
					change_ini=ini;
					change_fin=fin;
					change_number=number;
				}
				else if ((change_ctrl==ctrl)&&(change_panel==active_panel))
				{
					SetCtrlAttribute (active_panel, change_ctrl, ATTR_TEXT_BGCOLOR, VAL_WHITE);
					change_switch=0;
				}
			}
			else
			{
				if ((ini!=fin)&&(step!=0))
				{
					SetCtrlAttribute (active_panel, ctrl, ATTR_TEXT_BGCOLOR, VAL_RED);
					change_ctrl=ctrl;
					change_ini=ini;
					change_fin=fin;
					change_step=step;
					change_number=number;
					if (active_panel==panelSet) change_switch=1;
					if (active_panel==panelM_B) change_switch=2;
					change_panel=active_panel;
				}
			}
	}
	return 0;
}

int CVICALLBACK Set_range (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	double ini,fin;
	int number;
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_INI,&ini);
			GetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_FIN,&fin);
			GetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_NUMBER,&number);
			if (number!=0) SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_STEP,((fin-ini)/number));
			else           SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_STEP,0.0);
	}
	return 0;
}

int CVICALLBACK Set_step (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	double ini,fin,step;
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_INI,&ini);
			GetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_FIN,&fin);
			GetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_STEP,&step);
			if (step!=0)
			{
				 SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_NUMBER,(int)ceil((fin-ini)/step));
				 SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_STEP,((fin-ini)/ceil((fin-ini)/step)));
			}
			else SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_NUMBER,0);
	}
	return 0;
}

int CVICALLBACK Set_number (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	double ini,fin;
	int number;
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_INI,&ini);
			GetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_FIN,&fin);
			GetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_NUMBER,&number);
			if (number!=0) SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_STEP,((fin-ini)/number));
			else           SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_STEP,0.0);
	}
	return 0;
}

int CVICALLBACK Analysis (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{  // =========== calculates X and Y profiles and RMS values ============== //
	int fy[500] = {0}, fx[500] = {0}, i, j, k, bg_white, max, nn;
	double y,dy,x,dx, mean_x, mean_y, std_x, std_y;
	char txt[50];
	switch (event)
	{
		case EVENT_COMMIT:
			SetCtrlAttribute (panelS, PANEL_S_CANVAS, ATTR_PEN_COLOR, VAL_RED);
			SetCtrlAttribute (panelS, PANEL_S_CANVAS, ATTR_PEN_WIDTH, 2);
			CanvasSetPenPosition (panelS, PANEL_S_CANVAS, MakePoint(0,500));
			max = 0;
			dy = (fl/par.step)*(par.y_max - par.y_min)/(height/1);
			for (i=0; i<height/1; i++)
			{
				y = par.y_min*(fl/par.step) + i*dy;
				for(k=0;k<par.t_number;k++)
				{
					if(beam[k].status==on_screen)
					{
						if((beam[k].y>=y)&&(beam[k].y<y+dy)&&(beam[k].x>=par.x_min*fl/par.step)&&(beam[k].x<par.x_max*fl/par.step))
						{
							fy[i]+=1;
							if (max<fy[i]) max = fy[i];
						}
					}
				}
			}
			for (i=0; i<height/1; i++) CanvasDrawLineTo (panelS, PANEL_S_CANVAS, MakePoint(100.0*fy[i]/max,height-i));
			SetCtrlAttribute (panelS, PANEL_S_CANVAS, ATTR_PEN_COLOR, VAL_BLUE);
			CanvasSetPenPosition (panelS, PANEL_S_CANVAS, MakePoint(0,500));
			max = 0;
			dx = (fl/par.step)*(par.x_max - par.x_min)/(width/1);
			for (j=0; j<width/1; j++)
			{
				x = par.x_min*(fl/par.step) + j*dx;
				for(k=0;k<par.t_number;k++)
				{
					if(beam[k].status==on_screen)
					{
						if((beam[k].x>=x)&&(beam[k].x<x+dx)&&(beam[k].y>=par.y_min*fl/par.step)&&(beam[k].y<par.y_max*fl/par.step))
						{
							fx[j]+=1;
							if (max<fx[j]) max = fx[j];
						}
					}
				}
			}
			for (j=0; j<height/1; j++) CanvasDrawLineTo (panelS, PANEL_S_CANVAS, MakePoint(j,500-100.0*fx[j]/max));
			mean_x = 0; mean_y = 0; nn = 0;
			for(k=0;k<par.t_number;k++)
			{
				if(beam[k].status==on_screen)
				{
					mean_x+=beam[k].x;
					mean_y+=beam[k].y;
					nn += 1;
				}
			}
			mean_x/=nn;mean_y/=nn;
			std_x = 0; std_y = 0;
			for(k=0;k<par.t_number;k++)
			{
				if(beam[k].status==on_screen)
				{
					std_x+=pow(beam[k].x-mean_x,2);
					std_y+=pow(beam[k].y-mean_y,2);
				}
			}
			std_x = sqrt(std_x/nn); std_y = sqrt(std_y/nn);
			sprintf(txt, " xRMS = %.3f mm \n yRMS = %.3f mm ", 10*std_x*par.step/fl, 10*std_y*par.step/fl);
			GetCtrlVal(panelIM_SET,P_IM_SET_GRID_2,&bg_white);
			if (bg_white) SetCtrlAttribute (panelS, PANEL_S_CANVAS, ATTR_PEN_COLOR, VAL_BLACK);
			else 		  SetCtrlAttribute (panelS, PANEL_S_CANVAS, ATTR_PEN_COLOR, VAL_WHITE);
			SetCtrlAttribute (panelS, PANEL_S_CANVAS, ATTR_PEN_FILL_COLOR, VAL_TRANSPARENT);
			CanvasDrawTextAtPoint (panelS, PANEL_S_CANVAS, txt, VAL_MENU_META_FONT,MakePoint(10,10), VAL_UPPER_LEFT);
	}
	return 0;
}

int CVICALLBACK Test (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	double X[100]={0},Y[100]={0},k,d; int i;
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal(panelA,PANEL_A_K,&k);
			GetCtrlVal(panelA,PANEL_A_D,&d);
			for (i=0;i<100;i++)
			{
				X[i]=-obj[7].parametrs[3]+i*2*obj[7].parametrs[3]/99+d;
				Y[i]=k*sqrt(pow(obj[7].parametrs[3],2)-pow(X[i]-d,2))/obj[7].parametrs[3];
			}
			PlotXY (panelA, PANEL_A_GRAPH, X, Y, 100, VAL_DOUBLE, VAL_DOUBLE, VAL_FAT_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_CYAN);
	}
	return 0;
}

int CVICALLBACK QuitCallback_A (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			HidePanel(panelA);
			DeleteGraphPlot (panelA, PANEL_A_GRAPH, -1, VAL_IMMEDIATE_DRAW);
	}
	return 0;
}

int CVICALLBACK Ring_beam (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	FILE *F;
	char str[100]={0}; int i,j; double a;
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal(panelSet, PANEL_SET_RING_T_S,&par.gun_mode);
			switch (par.gun_mode)
			{
				case 0: SetCtrlAttribute (panelSet, PANEL_SET_LENGTH_G, ATTR_DIMMED, 1); break; // no gun
				case 1: // gun field from file
					SetCtrlAttribute (panelSet, PANEL_SET_LENGTH_G, ATTR_DIMMED, 0);
					F=fopen("gun_field.dat","r");
					for (j=0;j<2;j++)
					{
						fgets(str,100,F);
						for (i=0;i<100;i++)
						{
							fgets(str,100,F);
							sscanf(str,"%lf %lf %lf %lf",&a,&a,&Er_gun[i][j],&Ez_gun[i][j]); 
						}
					}
			}							  
	}
	return 0;
}

int CVICALLBACK Plot_graph (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{ // =============== plot static electric and magnetic fields ===================== //
	#define number 1000 						// number of point along Z
	double X_arr[number],Y_arr[number], t_ob;   // time to the object
	int a, i, st; 								// particle_status
	struct particle b = {0};
	struct parametrs par_in;
	struct object*   obj_in;
	switch (event)
	{
		case EVENT_COMMIT:
			Get_parametrs();
			obj_in=obj_to_in(par,obj);
			par_in=par_to_in(par,obj);
			GetCtrlVal(panelG,PANEL_G_RING_FIELDS,&a); // which fields to plot
			GetCtrlVal(panelG,PANEL_G_FIELD_X,&b.x);   // x coordinate
			GetCtrlVal(panelG,PANEL_G_FIELD_Y,&b.y);   // y coordinate
			b.x*=fl/par.step;
			b.y*=fl/par.step;
			for (i=0;i<number;i++)
			{
				X_arr[i]=i*par.screen/number;
				b.z=i*par_in.screen/number;
				st=where_is_particle(b,par_in,obj_in);
				t_ob = (obj_in[st].parametrs[0]+par_in.t_length/2)/par_in.t_v;
				switch (a)
				{ // which fields to plot
					case 0: Y_arr[i]=H_obj(t_ob, b, par_in, obj_in[st]).x	   /(1000.0*fe*par_in.step); break;
					case 1: Y_arr[i]=H_obj(t_ob, b, par_in, obj_in[st]).y	   /(1000.0*fe*par_in.step); break;
					case 2: Y_arr[i]=H_obj(t_ob, b, par_in, obj_in[st]).z	   /(1000.0*fe*par_in.step); break;
					case 3: Y_arr[i]=E_obj(t_ob, b, par_in, obj_in[st]).x*30000/(1000.0*fe*par_in.step); break;
					case 4: Y_arr[i]=E_obj(t_ob, b, par_in, obj_in[st]).y*30000/(1000.0*fe*par_in.step); break;
					case 5: Y_arr[i]=E_obj(t_ob, b, par_in, obj_in[st]).z*30000/(1000.0*fe*par_in.step); 
				}
			}
			DeleteGraphPlot (panelG, PANEL_G_GRAPH, -1, VAL_IMMEDIATE_DRAW);
			PlotXY (panelG, PANEL_G_GRAPH, X_arr, Y_arr, number, VAL_DOUBLE, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_BLUE);
			DisplayPanel(panelG);
	}
	return 0;
}

int CVICALLBACK QuitCallback_graph (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event) { case EVENT_COMMIT: HidePanel(panelG); }
	return 0;
}

int CVICALLBACK Slide (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	double data;
	switch (event)
	{
		case EVENT_VAL_CHANGED:
			GetCtrlVal(panel, PANEL_SET_NUMERICSLIDE_1,&data); 		obj[7].K_1 = data;
			SetCtrlVal(panel, PANEL_SET_NUMERICSLIDE_2,100-data); 	obj[7].K_2 = 100 - data; 
	}
	return 0;
}

int CVICALLBACK Charge_1 (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	double np1;
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal(panel, PANEL_SET_BUNCH_N  ,&np1);
			SetCtrlVal(panel, PANEL_SET_BUNCH_N_2,np1*1.6e-7);
	}
	return 0;
}

int CVICALLBACK Charge_2 (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	double charge;
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal(panel, PANEL_SET_BUNCH_N_2,&charge);
			SetCtrlVal(panel, PANEL_SET_BUNCH_N  ,charge/1.6e-7);
	}
	return 0;
}

int CVICALLBACK Save_config (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	int year,month,day,hours,minutes,second;
	char file_name[MAX_PATHNAME_LEN];
	switch (event)
	{
		case EVENT_COMMIT:
			GetSystemDate (&month, &day, &year);
			GetSystemTime (&hours, &minutes, &second);
			sprintf (file_name, "save\\%d_%02d_%02d_%02d_%02d_%02d.txt",year,month,day,hours,minutes,second);
			FileSelectPopup ("", file_name, "", "Save config", VAL_SAVE_BUTTON, 0, 0, 1, 1, file_name);
			save_config(&par,obj,file_name);
	}
	return 0;
}

int CVICALLBACK Load_config2 (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	char file_name[MAX_PATHNAME_LEN];
	switch (event)
	{
		case EVENT_COMMIT:
			FileSelectPopup ("", "config.txt", "", "Load config", VAL_LOAD_BUTTON, 0, 0, 1, 1, file_name);
			load_config(&par,file_name);
			set_config(par,obj);
	}
	return 0;
}

int CVICALLBACK Phase_space (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{// ============== plot phase space ===================== //
	int i, nn=0; double mean=0;
	switch (event)
	{
		case EVENT_COMMIT:
			DeleteGraphPlot (panelPHSP, PANEL_PHSP_GRAPH_X,-1, VAL_IMMEDIATE_DRAW);
			DeleteGraphPlot (panelPHSP, PANEL_PHSP_GRAPH_Y,-1, VAL_IMMEDIATE_DRAW);
			DeleteGraphPlot (panelPHSP, PANEL_PHSP_GRAPH_Z,-1, VAL_IMMEDIATE_DRAW);
			for (i=0;i<par.t_number;i++)
			{
				PlotPoint (panelPHSP, PANEL_PHSP_GRAPH_X, 10*beam[i].x/(fl/par.step), 1000*beam[i].px/beam[i].pz, VAL_SIMPLE_DOT, VAL_GREEN);
				PlotPoint (panelPHSP, PANEL_PHSP_GRAPH_Y, 10*beam[i].y/(fl/par.step), 1000*beam[i].py/beam[i].pz, VAL_SIMPLE_DOT, VAL_GREEN);
			}
			for (i=0;i<par.t_number;i++)
			{
				mean += beam[i].z;//  /(fl/par.step);
				nn   += 1;
			}
			mean = (mean/nn)*par.step/fl;
			for (i=0;i<par.t_number;i++) PlotPoint (panelPHSP, PANEL_PHSP_GRAPH_Z, 10*(beam[i].z/(fl/par.step)-mean), beam[i].pz*me/1e3, VAL_SIMPLE_DOT, VAL_GREEN);
			DisplayPanel(panelPHSP);
	}
	return 0;
}

int CVICALLBACK QuitCallback_PHSP (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event) { case EVENT_COMMIT: HidePanel(panelPHSP); }
	return 0;
}

int CVICALLBACK Bunch_transverse_ring (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			GetCtrlVal (panelSet, PANEL_SET_RING_B_TRD_1,	&(obj[7].TRD_1));
			switch ( obj[7].TRD_1 )
			{
				case 0: case 1: case 2: case 3: SetCtrlAttribute(panelSet, PANEL_SET_N_FLAT, ATTR_VISIBLE, 0);
				break;
				case 4: case 5:					SetCtrlAttribute(panelSet, PANEL_SET_N_FLAT, ATTR_VISIBLE, 1);
			}
			break;
	}
	return 0;
}
