#define fl 33.356409519815204957557671447492
#define fe 1.758820175e-5
#define fp 9.578834102e-9
#define vc 29979245800.0
#define qe 4.803204400786e-10
#define pi 3.141592653589793
#define me 510998.918
#define mp 938272029.0

#define height 500  //vertical screen size
#define width  500  //horizontal screen size
#define number_of_object 17

#define in_quad 		 0
#define in_quad2 		 1
#define in_quad3 		 2
#define in_quad4 		 3
#define in_lens 		 4
#define in_scan 		 5
#define in_RF 			 6
#define in_bunch 		 7
#define in_corr	   		 8
#define in_acc	 		 9
#define in_cathode 		10
#define in_gun 			11
#define in_free_space	12
#define on_tube 		13
#define on_screen 		14
#define on_diaphragm	15
#define reflection 		16
#define on_stop         17
/*typedef enum {in_cathode, in_gun, in_lens, in_RF, in_scan, in_quad, in_quad2, in_quad3, in_quad4, 
				in_bunch, in_corr, in_free_space, on_tube, on_screen, 
				on_diaphragm, reflection, in_acc} particle_status; */
double Er_gun[100][2],Ez_gun[100][2];

int num;  // flying particles
int panelS;
int panelIM_SET;
int panelSet;
int panelM_B;
int panelCHNG;
int panelA;
int panelG; 
int panelPHSP;
int panel_3D;

int change_ctrl, change_number, change_switch, change_panel;
double change_ini, change_fin, change_step;

struct particle  
{
	double x,y,z;     //struct vector c;
	double px,py,pz;  //struct vector p;
	int status;
} *beam;

struct vector
{
	double x, y, z;
};

struct gauss_bunch //летит по "x"
{
	int n;
	double *t, *i, *s, *y, *z;
};

struct parametrs
{
	//testing beam
	int    gun_mode;			//способ задания пушки
	double gun_length;			//длина пушки (cm)
	double t_radius;		   	//поперечный размер тестового пучка (cm)
	double t_length;		   	//продольный размер тестового пучка (ps)
	double t_temp;			   	//температура тестового пучка (eV)
	int    t_number;			//число частиц тестового пучка
	double t_energy;		   	//энергия тестового пучка (keV)
	double t_spread;			//energy spread of the testing beam [keV]
	unsigned int t_distr;		//форма распределения тестового пучка
	double t_v;					//скорость тестового пучка
	double space_charge;		//пространственный заряд
	double t_weight;			//вес частицы (сколько электронов)
	//bunch
	int cb;						//Способ задания банча
	unsigned short int N_Flat;	//Number of Gaussian distributions
	double b_v;					//скорость сгустка
	//geometry				   
	double radius;			    //радиус трубы (cm)
	double screen;			    //положение экрана (cm)
	double stop;				//stop all particles
	double diaph;				//положение диафрагмы
	double diaph_r;				//радиус диафрагмы
	//option				    
	double step;			    //шаг по времени (ns)
	double scan_mode;	    	//тип развёртки
	//image
	double x_min, x_max;	    //screen range
	double y_min, y_max;
	int tikalka;  				//Число шагов без обновления окна
	int t_number_draw;			//number of particles for drawing
	int graph_mode;				//Вид отображения
	double a_y;					//Коэффициент растяжения по вертикали
	int sv;						//сохранение движения пучка
	int clear_graph;			//стирание следа частицы
	int plot_2D;				//plot particles in 3D
	int plot_3D;				//plot particles in 3D
} par;

struct object
{
	int name;
	double parametrs[10];
	//struct vector (*E)(int,struct particle,struct parametrs,struct object);
	//struct vector (*H)(int,struct particle,struct parametrs,struct object);
	// struct vector (*E)(int,struct particle); // To be realized
	// struct vector (*H)(int,struct particle); // To be realized
	int    TRD_1;
	int    TRD_2;
	double L_1, L_2;
	double X_1, X_2;
	double Y_1, Y_2;
	double Z_1, Z_2;
	double R_1, R_2;
	double K_1, K_2;
} *obj;

#define QNX 25
#define QNY 25
#define QNZ 20

struct vector quadrupole_fields1 [QNX+1][QNY+1][QNZ+1];
struct vector quadrupole_fields2 [QNX+1][QNY+1][QNZ+1];
struct vector rf_magnetic_fields [11][101];
struct vector rf_electric_fields [11][101];

struct q_p
{
	double size_x, size_y, size_z;
	double size_dx,size_dy,size_dz;
	int sw; // прочитали или нет квадруполь из файла
} quadrupole_parametrs;
