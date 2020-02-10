#include <ansi_c.h>
#include <utility.h>
#include <userint.h>
#include "simulator_beam.h"
#include "simulator_vars.h"
#include "make_bunch_my.h"

#define heightmb 350
#define widthmb 500

static int plot_sigma,plot_delta,auto_refresh;  
static unsigned char imbits[heightmb*widthmb];


void Set_zero_color(void)
{
int NumCtrl,Ctrl,i,q,tmp_int;
double tmp;
	GetPanelAttribute (panelM_B, ATTR_NUM_CTRLS, &NumCtrl);
	GetPanelAttribute (panelM_B, ATTR_PANEL_FIRST_CTRL, &Ctrl);
	for (i=0;i<33;i++)
	{
		GetCtrlAttribute (panelM_B, Ctrl, ATTR_DATA_TYPE,&q);
		if(q==VAL_DOUBLE)
		{
			GetCtrlVal (panelM_B, Ctrl,&tmp);
			if (fabs(tmp)<1e-10) 
				SetCtrlAttribute (panelM_B,Ctrl, ATTR_TEXT_BGCOLOR,VAL_WHITE);
			else
				SetCtrlAttribute (panelM_B,Ctrl, ATTR_TEXT_BGCOLOR,VAL_GREEN);
		}
		if(q==VAL_UNSIGNED_INTEGER)
		{
			GetCtrlVal (panelM_B, Ctrl,&tmp_int);
			if (tmp_int==0) 
				SetCtrlAttribute (panelM_B,Ctrl, ATTR_TEXT_BGCOLOR,VAL_WHITE);
			else
				SetCtrlAttribute (panelM_B,Ctrl, ATTR_TEXT_BGCOLOR,VAL_GREEN);
		}
		GetCtrlAttribute (panelM_B, Ctrl, ATTR_NEXT_CTRL, &Ctrl);
	}

}

void Set_zero_color_one_control(void)
{
int AP,AC,q,tmp_int;
double tmp;
	
	AP = GetActivePanel ();
	AC = GetActiveCtrl (AP);
	
	GetCtrlAttribute (AP,AC, ATTR_DATA_TYPE,&q);
	if(q==VAL_DOUBLE)
	{
		GetCtrlVal (AP, AC,&tmp);
		if (fabs(tmp)<1e-10)  SetCtrlAttribute (AP,AC, ATTR_TEXT_BGCOLOR,VAL_WHITE);
		else SetCtrlAttribute (AP,AC, ATTR_TEXT_BGCOLOR,VAL_GREEN);
	}
	if(q==VAL_UNSIGNED_INTEGER)
	{
		GetCtrlVal (AP, AC,&tmp_int);
		if (tmp_int==0) SetCtrlAttribute (AP,AC, ATTR_TEXT_BGCOLOR,VAL_WHITE);
		else SetCtrlAttribute (AP,AC, ATTR_TEXT_BGCOLOR,VAL_GREEN);
	}
}


int CVICALLBACK RefreshM_B (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event) 
	{
		case EVENT_COMMIT:
			//auto_display_format(panel,control,0.001,1000);
			Refresh_rho();
			Refresh_sigma();
			Refresh_delta();
			Draw_bunch();
		break;
	}
	return 0;
}

int CVICALLBACK Refresh_Rho (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event) 
	{
		case EVENT_COMMIT:
			Set_zero_color_one_control(); 
			Refresh_rho();
			GetCtrlVal(panelM_B,PANEL_M_B_AUTO_REFRESH,&auto_refresh);
			if (auto_refresh)Draw_bunch();
	}
	return 0;
}
int CVICALLBACK Refresh_Sigma (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event) {
		case EVENT_COMMIT:
			Set_zero_color_one_control(); 
			Refresh_sigma();
			GetCtrlVal(panelM_B,PANEL_M_B_AUTO_REFRESH,&auto_refresh);
			if (auto_refresh)Draw_bunch();
	}
	return 0;
}
int CVICALLBACK Refresh_Delta (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)       
{                                                                 
	switch (event) 
	{
		case EVENT_COMMIT:  
			Set_zero_color_one_control();
			Refresh_delta();
			GetCtrlVal(panelM_B,PANEL_M_B_AUTO_REFRESH,&auto_refresh);
			if (auto_refresh)Draw_bunch();
	}                                                             
	return 0;                                                     
}                                                                 

void Refresh_rho()
{
	double X_i;
	int i;
	
	GetCtrlVal(panelM_B,PANEL_M_B_NUM,&NUM);
	GetCtrlVal(panelM_B,PANEL_M_B_XMIN,&XMIN);
	GetCtrlVal(panelM_B,PANEL_M_B_XMAX,&XMAX);
	GetCtrlVal(panelM_B,PANEL_M_B_A_R,&matrix_make_bunch[0][0][0]);
	GetCtrlVal(panelM_B,PANEL_M_B_B_R,&matrix_make_bunch[0][1][0]);
	GetCtrlVal(panelM_B,PANEL_M_B_C_R,&matrix_make_bunch[0][2][0]);	
	GetCtrlVal(panelM_B,PANEL_M_B_D_R,&matrix_make_bunch[0][3][0]);
	GetCtrlVal(panelM_B,PANEL_M_B_O_R,&matrix_make_bunch[0][4][0]);
	GetCtrlVal(panelM_B,PANEL_M_B_Fi_R,&matrix_make_bunch[0][5][0]);	
	GetCtrlVal(panelM_B,PANEL_M_B_E_R,&matrix_make_bunch[0][6][0]);
	GetCtrlVal(panelM_B,PANEL_M_B_X0_R,&matrix_make_bunch[0][7][0]);
	GetCtrlVal(panelM_B,PANEL_M_B_S_R,&matrix_make_bunch[0][8][0]);
	GetCtrlVal(panelM_B,PANEL_M_B_F_R,&matrix_make_bunch[0][9][0]);	
	if (NUM<2) return;
	for (i=0;i<NUM;i++)
	{
		X[i]=XMIN+i*(XMAX-XMIN)/(NUM-1);
		X_i=X[i];//-1.0+i*2.0/(NUM-1);
		RHO[i]=(matrix_make_bunch[0][0][0]*X_i*X_i+matrix_make_bunch[0][1][0]*X_i+matrix_make_bunch[0][2][0])*(matrix_make_bunch[0][3][0]*sin(matrix_make_bunch[0][4][0]*X_i+matrix_make_bunch[0][5][0])+matrix_make_bunch[0][6][0]*exp(-(X_i-matrix_make_bunch[0][7][0])*(X_i-matrix_make_bunch[0][7][0])/(2*matrix_make_bunch[0][8][0]*matrix_make_bunch[0][8][0]))+matrix_make_bunch[0][9][0]);
		if (RHO[i]<0) RHO[i]=0;
	}
	DeleteGraphPlot (panelM_B, PANEL_M_B_GRAPH, -1, VAL_IMMEDIATE_DRAW);
	
	PlotXY (panelM_B, PANEL_M_B_GRAPH, X, RHO, NUM, VAL_DOUBLE, VAL_DOUBLE,
			VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_RED);
			
	Normirovka=0;
	X_i=XMIN;
	while (X_i<XMAX)
	{
		Normirovka+=par.step*(matrix_make_bunch[0][0][0]*X_i*X_i+matrix_make_bunch[0][1][0]*X_i+matrix_make_bunch[0][2][0])*
		                     (matrix_make_bunch[0][3][0]*sin(matrix_make_bunch[0][4][0]*X_i+matrix_make_bunch[0][5][0])+matrix_make_bunch[0][6][0]*exp(-(X_i-matrix_make_bunch[0][7][0])*(X_i-matrix_make_bunch[0][7][0])/(2*matrix_make_bunch[0][8][0]*matrix_make_bunch[0][8][0]))+matrix_make_bunch[0][9][0]);
		X_i+=par.step;
	}
	Normirovka=1/Normirovka;
}

void Refresh_sigma()
{
	double X_i;
	int i;
	
	GetCtrlVal(panelM_B,PANEL_M_B_Y_Z,&yz);
	GetCtrlVal(panelM_B,PANEL_M_B_NUM,&NUM);
	GetCtrlVal(panelM_B,PANEL_M_B_XMIN,&XMIN);
	GetCtrlVal(panelM_B,PANEL_M_B_XMAX,&XMAX);
	GetCtrlVal(panelM_B,PANEL_M_B_A_S,&matrix_make_bunch[1][0][yz]);
	GetCtrlVal(panelM_B,PANEL_M_B_B_S,&matrix_make_bunch[1][1][yz]);
	GetCtrlVal(panelM_B,PANEL_M_B_C_S,&matrix_make_bunch[1][2][yz]);	
	GetCtrlVal(panelM_B,PANEL_M_B_D_S,&matrix_make_bunch[1][3][yz]);
	GetCtrlVal(panelM_B,PANEL_M_B_O_S,&matrix_make_bunch[1][4][yz]);
	GetCtrlVal(panelM_B,PANEL_M_B_Fi_S,&matrix_make_bunch[1][5][yz]);	
	GetCtrlVal(panelM_B,PANEL_M_B_E_S,&matrix_make_bunch[1][6][yz]);
	GetCtrlVal(panelM_B,PANEL_M_B_X0_S,&matrix_make_bunch[1][7][yz]);
	GetCtrlVal(panelM_B,PANEL_M_B_S_S,&matrix_make_bunch[1][8][yz]);	
	GetCtrlVal(panelM_B,PANEL_M_B_F_S,&matrix_make_bunch[1][9][yz]);		
	if (NUM<2) return;
	for (i=0;i<NUM;i++)
	{
		X[i]=XMIN+i*(XMAX-XMIN)/(NUM-1);
		X_i=X[i];//-1.0+i*2.0/(NUM-1);		
		SIGMA[yz][i]=(matrix_make_bunch[1][0][yz]*X_i*X_i+matrix_make_bunch[1][1][yz]*X_i+matrix_make_bunch[1][2][yz])*(matrix_make_bunch[1][3][yz]*sin(matrix_make_bunch[1][4][yz]*X_i+matrix_make_bunch[1][5][yz])+matrix_make_bunch[1][6][yz]*exp(-(X_i-matrix_make_bunch[1][7][yz])*(X_i-matrix_make_bunch[1][7][yz])/(2*matrix_make_bunch[1][8][yz]*matrix_make_bunch[1][8][yz]))+matrix_make_bunch[1][9][yz]);
		if (SIGMA[yz][i]<0) SIGMA[yz][i]=0;
	}
	if (plot_sigma!=0) DeleteGraphPlot (panelM_B, PANEL_M_B_GRAPH_2, plot_sigma, VAL_IMMEDIATE_DRAW);
	
	plot_sigma=PlotXY (panelM_B, PANEL_M_B_GRAPH_2, X, SIGMA[yz], NUM, VAL_DOUBLE, VAL_DOUBLE,
			VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_BLUE);
}
void Refresh_delta()
{
	double X_i;
	int i;
	
	GetCtrlVal(panelM_B,PANEL_M_B_Y_Z,&yz);
	GetCtrlVal(panelM_B,PANEL_M_B_NUM,&NUM);
	GetCtrlVal(panelM_B,PANEL_M_B_XMIN,&XMIN);
	GetCtrlVal(panelM_B,PANEL_M_B_XMAX,&XMAX);
	GetCtrlVal(panelM_B,PANEL_M_B_A_D,&matrix_make_bunch[2][0][yz]);
	GetCtrlVal(panelM_B,PANEL_M_B_B_D,&matrix_make_bunch[2][1][yz]);
	GetCtrlVal(panelM_B,PANEL_M_B_C_D,&matrix_make_bunch[2][2][yz]);	
	GetCtrlVal(panelM_B,PANEL_M_B_D_D,&matrix_make_bunch[2][3][yz]);
	GetCtrlVal(panelM_B,PANEL_M_B_O_D,&matrix_make_bunch[2][4][yz]);
	GetCtrlVal(panelM_B,PANEL_M_B_Fi_D,&matrix_make_bunch[2][5][yz]);	
	GetCtrlVal(panelM_B,PANEL_M_B_E_D,&matrix_make_bunch[2][6][yz]);
	GetCtrlVal(panelM_B,PANEL_M_B_X0_D,&matrix_make_bunch[2][7][yz]);
	GetCtrlVal(panelM_B,PANEL_M_B_S_D,&matrix_make_bunch[2][8][yz]);	
	GetCtrlVal(panelM_B,PANEL_M_B_F_D,&matrix_make_bunch[2][9][yz]);		
	if (NUM<2) return;
	for (i=0;i<NUM;i++)
	{
		X[i]=XMIN+i*(XMAX-XMIN)/(NUM-1);
		X_i=X[i];//-1.0+i*2.0/(NUM-1);		
		DELTA[yz][i]=(matrix_make_bunch[2][0][yz]*X_i*X_i+matrix_make_bunch[2][1][yz]*X_i+matrix_make_bunch[2][2][yz])*(matrix_make_bunch[2][3][yz]*sin(matrix_make_bunch[2][4][yz]*X_i+matrix_make_bunch[2][5][yz])+matrix_make_bunch[2][6][yz]*exp(-(X_i-matrix_make_bunch[2][7][yz])*(X_i-matrix_make_bunch[2][7][yz])/(2*matrix_make_bunch[2][8][yz]*matrix_make_bunch[2][8][yz]))+matrix_make_bunch[2][9][yz]);
	}
	if (plot_delta!=0) DeleteGraphPlot (panelM_B, PANEL_M_B_GRAPH_2, plot_delta, VAL_IMMEDIATE_DRAW);
	
	plot_delta=PlotXY (panelM_B, PANEL_M_B_GRAPH_2, X, DELTA[yz], NUM, VAL_DOUBLE, VAL_DOUBLE,
			VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, 0x00e000);
}

void Draw_bunch(void)
{

int ctmp_imbits;
int i,j;

static float tmp_imbits[heightmb*widthmb];
int bitmp,index_x,s_yz;
double delta_plus_2_sigma_max,delta_plus_2_sigma_min,y,tmp_imbitsmax,num=3;
	
	delta_plus_2_sigma_max= fabs(DELTA[yz][0])+fabs(num*SIGMA[yz][0]);
	delta_plus_2_sigma_min=-fabs(DELTA[yz][0])-fabs(num*SIGMA[yz][0]);	
			
	for(i=1;i<widthmb;i++) 
	{
		if(delta_plus_2_sigma_max<( fabs(DELTA[yz][i])+fabs(num*SIGMA[yz][i]))) delta_plus_2_sigma_max= fabs(DELTA[yz][i])+fabs(num*SIGMA[yz][i]);
		if(delta_plus_2_sigma_min>(-fabs(DELTA[yz][i])-fabs(num*SIGMA[yz][i]))) delta_plus_2_sigma_min=-fabs(DELTA[yz][i])-fabs(num*SIGMA[yz][i]);
	}
	for(i=0;i<widthmb;i++)
	{
		index_x=(int)(i*NUM/widthmb);
		for(j=0;j<heightmb;j++)
		{
			y=delta_plus_2_sigma_max+j*(delta_plus_2_sigma_min-delta_plus_2_sigma_max)/(1.0*heightmb);
			tmp_imbits[i+widthmb*j]=RHO[index_x]*exp(-(y-DELTA[yz][index_x])*(y-DELTA[yz][index_x])/(2*SIGMA[yz][index_x]*SIGMA[yz][index_x]));
		}
	}
	tmp_imbitsmax=0;
	for(i=0;i<widthmb*heightmb;i++) if(tmp_imbitsmax<tmp_imbits[i]) tmp_imbitsmax=tmp_imbits[i];
	if(tmp_imbitsmax!=0) ctmp_imbits=255/tmp_imbitsmax; 
	else ctmp_imbits=1;
	
	for(i=0;i<widthmb*heightmb;i++)	imbits[i]=(unsigned char)(ctmp_imbits*tmp_imbits[i]);
	
	//DiscardBitmap (bitmp);
	NewBitmap (widthmb, 8, widthmb, heightmb, CT, imbits, NULL, &bitmp);
	CanvasDrawBitmap (panelM_B, PANEL_M_B_CANVAS, bitmp, VAL_ENTIRE_OBJECT,VAL_ENTIRE_OBJECT);
	GetCtrlVal(panelM_B,PANEL_M_B_Y_Z,&s_yz);
	//GetCtrlAttribute (panelM_B, PANEL_M_B_CANVAS, ATTR_HEIGHT, &c_height);
	if (s_yz==0)
	{
		
		SetCtrlAttribute (panelM_B, PANEL_M_B_CANVAS, ATTR_PEN_FILL_COLOR, 0);
		SetCtrlAttribute (panelM_B, PANEL_M_B_CANVAS, ATTR_PEN_COLOR,0xffff00);
		//CanvasDrawLine (panelM_B, PANEL_M_B_CANVAS, MakePoint(10,0),MakePoint(10,c_height));
		CanvasDrawTextAtPoint (panelM_B, PANEL_M_B_CANVAS, "X-Y",VAL_APP_META_FONT, MakePoint (20, 20),VAL_CENTER);
	}
	else
	{
		SetCtrlAttribute (panelM_B, PANEL_M_B_CANVAS, ATTR_PEN_FILL_COLOR, 0); 
		SetCtrlAttribute (panelM_B, PANEL_M_B_CANVAS, ATTR_PEN_COLOR,0xff00ff);
		CanvasDrawTextAtPoint (panelM_B, PANEL_M_B_CANVAS, "X-Z",VAL_APP_META_FONT, MakePoint (20, 20),VAL_CENTER);
	}
	DiscardBitmap (bitmp);	
}

int CVICALLBACK Auto_refresh (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	int status;
	
	switch (event) 
	{
		case EVENT_COMMIT:
			GetCtrlVal(panelM_B,PANEL_M_B_AUTO_REFRESH,&status);
			if (status==1)
			{
				Set_zero_color(); 
				Refresh_rho();
				Refresh_sigma();
				Refresh_delta();
				Draw_bunch();
			}
	}
	return 0;
}

void auto_display_format(int panel,int control,double min,double max)
{
	int style,data_type,val_int;
	double val_double;
	
	GetCtrlAttribute (panel, control, ATTR_CTRL_STYLE, &style);
	if ((style==CTRL_NUMERIC)||(style==CTRL_NUMERIC_LS))
	{
		GetCtrlAttribute (panel, control, ATTR_DATA_TYPE, &data_type);
		switch (data_type)
		{
			case VAL_INTEGER:
			case VAL_SHORT_INTEGER:
			case VAL_UNSIGNED_INTEGER:
			case VAL_UNSIGNED_SHORT_INTEGER:
				GetCtrlVal(panel,control,&val_int);
				if(val_int==0) 									 SetCtrlAttribute (panel, control, ATTR_FORMAT, VAL_FLOATING_PT_FORMAT);
				else if ((abs(val_int)<min)||(abs(val_int)>max)) SetCtrlAttribute (panel, control, ATTR_FORMAT, VAL_SCIENTIFIC_FORMAT);
				else 											 SetCtrlAttribute (panel, control, ATTR_FORMAT, VAL_FLOATING_PT_FORMAT);
			break;
			case  (VAL_DOUBLE):
				GetCtrlVal(panel,control,&val_double);
				if(val_double==0) 										 SetCtrlAttribute (panel, control, ATTR_FORMAT, VAL_FLOATING_PT_FORMAT);
				else if ((fabs(val_double)<min)||(fabs(val_double)>max)) SetCtrlAttribute (panel, control, ATTR_FORMAT, VAL_SCIENTIFIC_FORMAT);
				else 													 SetCtrlAttribute (panel, control, ATTR_FORMAT, VAL_FLOATING_PT_FORMAT);
		}
	}
}
