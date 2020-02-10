#include <ansi_c.h>
#include <userint.h>
#include "simulator_beam.h"
#include "simulator_vars.h"
#include "make_bunch_my.h"

void Set_y_z(void)
{
	int status_y_z;
	GetCtrlVal(panelM_B,PANEL_M_B_Y_Z,&status_y_z);
	if (status_y_z==0)
	{
		SetCtrlAttribute (panelM_B, PANEL_M_B_TEXTMSG_SIGMA, ATTR_TEXT_COLOR,0xffff00);
		SetCtrlVal(panelM_B,PANEL_M_B_TEXTMSG_SIGMA,"Y =");
		SetCtrlAttribute (panelM_B, PANEL_M_B_TEXTMSG_DELTA, ATTR_TEXT_COLOR,0xffff00);		
		SetCtrlVal(panelM_B,PANEL_M_B_TEXTMSG_DELTA,"Y =");
	}
	else
	{
		SetCtrlAttribute (panelM_B, PANEL_M_B_TEXTMSG_SIGMA, ATTR_TEXT_COLOR,0xff00ff);
		SetCtrlVal(panelM_B,PANEL_M_B_TEXTMSG_SIGMA,"Z =");
		SetCtrlAttribute (panelM_B, PANEL_M_B_TEXTMSG_DELTA, ATTR_TEXT_COLOR,0xff00ff);
		SetCtrlVal(panelM_B,PANEL_M_B_TEXTMSG_DELTA,"Z =");
	}
}

int CVICALLBACK y_z (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event) 
	{
		case EVENT_COMMIT:
			GetCtrlVal(panelM_B,PANEL_M_B_Y_Z,&yz);
			Set_y_z();
			SetCtrlVal(panelM_B,PANEL_M_B_A_S,matrix_make_bunch[1][0][yz]);
			SetCtrlVal(panelM_B,PANEL_M_B_B_S,matrix_make_bunch[1][1][yz]);   	
			SetCtrlVal(panelM_B,PANEL_M_B_C_S,matrix_make_bunch[1][2][yz]);
			SetCtrlVal(panelM_B,PANEL_M_B_D_S,matrix_make_bunch[1][3][yz]);
			SetCtrlVal(panelM_B,PANEL_M_B_O_S,matrix_make_bunch[1][4][yz]);
			SetCtrlVal(panelM_B,PANEL_M_B_Fi_S,matrix_make_bunch[1][5][yz]);
			SetCtrlVal(panelM_B,PANEL_M_B_E_S,matrix_make_bunch[1][6][yz]);   	
			SetCtrlVal(panelM_B,PANEL_M_B_X0_S,matrix_make_bunch[1][7][yz]);
			SetCtrlVal(panelM_B,PANEL_M_B_S_S,matrix_make_bunch[1][8][yz]);
			SetCtrlVal(panelM_B,PANEL_M_B_F_S,matrix_make_bunch[1][9][yz]);
			SetCtrlVal(panelM_B,PANEL_M_B_A_D,matrix_make_bunch[2][0][yz]);
			SetCtrlVal(panelM_B,PANEL_M_B_B_D,matrix_make_bunch[2][1][yz]);   	
			SetCtrlVal(panelM_B,PANEL_M_B_C_D,matrix_make_bunch[2][2][yz]);
			SetCtrlVal(panelM_B,PANEL_M_B_D_D,matrix_make_bunch[2][3][yz]);
			SetCtrlVal(panelM_B,PANEL_M_B_O_D,matrix_make_bunch[2][4][yz]);
			SetCtrlVal(panelM_B,PANEL_M_B_Fi_D,matrix_make_bunch[2][5][yz]);
			SetCtrlVal(panelM_B,PANEL_M_B_E_D,matrix_make_bunch[2][6][yz]);
			SetCtrlVal(panelM_B,PANEL_M_B_X0_D,matrix_make_bunch[2][7][yz]);
			SetCtrlVal(panelM_B,PANEL_M_B_S_D,matrix_make_bunch[2][8][yz]);
			SetCtrlVal(panelM_B,PANEL_M_B_F_D,matrix_make_bunch[2][9][yz]);
			Set_zero_color();			
			Refresh_rho();
			Refresh_sigma();
			Refresh_delta();
			Draw_bunch();
		break;
	}
	return 0;
}

void RecallPS (int P, char* fn)
{
	FILE* F_cfg;
	int auto_refresh;

	F_cfg = fopen (fn, "r");
	fscanf (F_cfg, "%d\n", &NUM);
	fscanf(F_cfg,"%lf\n",&XMIN);
	fscanf(F_cfg,"%lf\n",&XMAX);
	fscanf (F_cfg, "%lf %lf %lf %lf %lf\n", &matrix_make_bunch[0][0][0], &matrix_make_bunch[1][0][0], &matrix_make_bunch[2][0][0], &matrix_make_bunch[1][0][1], &matrix_make_bunch[2][0][1]);   
	fscanf (F_cfg, "%lf %lf %lf %lf %lf\n", &matrix_make_bunch[0][1][0], &matrix_make_bunch[1][1][0], &matrix_make_bunch[2][1][0], &matrix_make_bunch[1][1][1], &matrix_make_bunch[2][1][1]);   
	fscanf (F_cfg, "%lf %lf %lf %lf %lf\n", &matrix_make_bunch[0][2][0], &matrix_make_bunch[1][2][0], &matrix_make_bunch[2][2][0], &matrix_make_bunch[1][2][1], &matrix_make_bunch[2][2][1]);	
	fscanf (F_cfg, "%lf %lf %lf %lf %lf\n", &matrix_make_bunch[0][3][0], &matrix_make_bunch[1][3][0], &matrix_make_bunch[2][3][0], &matrix_make_bunch[1][3][1], &matrix_make_bunch[2][3][1]);   
	fscanf (F_cfg, "%lf %lf %lf %lf %lf\n", &matrix_make_bunch[0][4][0], &matrix_make_bunch[1][4][0], &matrix_make_bunch[2][4][0], &matrix_make_bunch[1][4][1], &matrix_make_bunch[2][4][1]);   
	fscanf (F_cfg, "%lf %lf %lf %lf %lf\n", &matrix_make_bunch[0][5][0], &matrix_make_bunch[1][5][0], &matrix_make_bunch[2][5][0], &matrix_make_bunch[1][5][1], &matrix_make_bunch[2][5][1]);	
	fscanf (F_cfg, "%lf %lf %lf %lf %lf\n", &matrix_make_bunch[0][6][0], &matrix_make_bunch[1][6][0], &matrix_make_bunch[2][6][0], &matrix_make_bunch[1][6][1], &matrix_make_bunch[2][6][1]);   
	fscanf (F_cfg, "%lf %lf %lf %lf %lf\n", &matrix_make_bunch[0][7][0], &matrix_make_bunch[1][7][0], &matrix_make_bunch[2][7][0], &matrix_make_bunch[1][7][1], &matrix_make_bunch[2][7][1]); 
	fscanf (F_cfg, "%lf %lf %lf %lf %lf\n", &matrix_make_bunch[0][8][0], &matrix_make_bunch[1][8][0], &matrix_make_bunch[2][8][0], &matrix_make_bunch[1][8][1], &matrix_make_bunch[2][8][1]);   
	fscanf (F_cfg, "%lf %lf %lf %lf %lf\n", &matrix_make_bunch[0][9][0], &matrix_make_bunch[1][9][0], &matrix_make_bunch[2][9][0], &matrix_make_bunch[1][9][1], &matrix_make_bunch[2][9][1]);
	fscanf (F_cfg, "%d\n", &yz);
	fscanf (F_cfg, "%d\n", &auto_refresh);
	SetCtrlVal(panelM_B,PANEL_M_B_Y_Z,yz);   	
	SetCtrlVal(panelM_B,PANEL_M_B_NUM,NUM);   	
	SetCtrlVal(panelM_B,PANEL_M_B_XMIN,XMIN); 	
	SetCtrlVal(panelM_B,PANEL_M_B_XMAX,XMAX);
	SetCtrlVal(panelM_B,PANEL_M_B_A_R, matrix_make_bunch[0][0][0]);   	
	SetCtrlVal(panelM_B,PANEL_M_B_B_R, matrix_make_bunch[0][1][0]);   	
	SetCtrlVal(panelM_B,PANEL_M_B_C_R, matrix_make_bunch[0][2][0]);		
	SetCtrlVal(panelM_B,PANEL_M_B_D_R, matrix_make_bunch[0][3][0]);   	
	SetCtrlVal(panelM_B,PANEL_M_B_O_R, matrix_make_bunch[0][4][0]);   	
	SetCtrlVal(panelM_B,PANEL_M_B_Fi_R,matrix_make_bunch[0][5][0]);		
	SetCtrlVal(panelM_B,PANEL_M_B_E_R, matrix_make_bunch[0][6][0]);   	
	SetCtrlVal(panelM_B,PANEL_M_B_X0_R,matrix_make_bunch[0][7][0]); 	
	SetCtrlVal(panelM_B,PANEL_M_B_S_R, matrix_make_bunch[0][8][0]);   	
	SetCtrlVal(panelM_B,PANEL_M_B_F_R, matrix_make_bunch[0][9][0]);		
	SetCtrlVal(panelM_B,PANEL_M_B_A_S, matrix_make_bunch[1][0][yz]);   	
	SetCtrlVal(panelM_B,PANEL_M_B_B_S, matrix_make_bunch[1][1][yz]);   	
	SetCtrlVal(panelM_B,PANEL_M_B_C_S, matrix_make_bunch[1][2][yz]);		
	SetCtrlVal(panelM_B,PANEL_M_B_D_S, matrix_make_bunch[1][3][yz]);   	
	SetCtrlVal(panelM_B,PANEL_M_B_O_S, matrix_make_bunch[1][4][yz]);   	
	SetCtrlVal(panelM_B,PANEL_M_B_Fi_S,matrix_make_bunch[1][5][yz]);		
	SetCtrlVal(panelM_B,PANEL_M_B_E_S, matrix_make_bunch[1][6][yz]);   	
	SetCtrlVal(panelM_B,PANEL_M_B_X0_S,matrix_make_bunch[1][7][yz]); 	
	SetCtrlVal(panelM_B,PANEL_M_B_S_S, matrix_make_bunch[1][8][yz]);		
	SetCtrlVal(panelM_B,PANEL_M_B_F_S, matrix_make_bunch[1][9][yz]);	
	SetCtrlVal(panelM_B,PANEL_M_B_A_D, matrix_make_bunch[2][0][yz]);   	
	SetCtrlVal(panelM_B,PANEL_M_B_B_D, matrix_make_bunch[2][1][yz]);   	
	SetCtrlVal(panelM_B,PANEL_M_B_C_D, matrix_make_bunch[2][2][yz]);	
	SetCtrlVal(panelM_B,PANEL_M_B_D_D, matrix_make_bunch[2][3][yz]);   	
	SetCtrlVal(panelM_B,PANEL_M_B_O_D, matrix_make_bunch[2][4][yz]);   	
	SetCtrlVal(panelM_B,PANEL_M_B_Fi_D,matrix_make_bunch[2][5][yz]);	
	SetCtrlVal(panelM_B,PANEL_M_B_E_D, matrix_make_bunch[2][6][yz]);   	
	SetCtrlVal(panelM_B,PANEL_M_B_X0_D,matrix_make_bunch[2][7][yz]); 	
	SetCtrlVal(panelM_B,PANEL_M_B_S_D, matrix_make_bunch[2][8][yz]);	
	SetCtrlVal(panelM_B,PANEL_M_B_F_D, matrix_make_bunch[2][9][yz]);
	SetCtrlVal(panelM_B,PANEL_M_B_AUTO_REFRESH,auto_refresh);
}	


void SavePS (int P, char* fn)
{
	int auto_refresh;
	FILE* F_cfg;
	F_cfg = fopen (fn, "w");
	fprintf(F_cfg,"%d\n",NUM);
	fprintf(F_cfg,"%f\n",XMIN);
	fprintf(F_cfg,"%f\n",XMAX);
	fprintf(F_cfg, "%f %f %f %f %f\n", matrix_make_bunch[0][0][0], matrix_make_bunch[1][0][0], matrix_make_bunch[2][0][0], matrix_make_bunch[1][0][1], matrix_make_bunch[2][0][1]);   
	fprintf(F_cfg, "%f %f %f %f %f\n", matrix_make_bunch[0][1][0], matrix_make_bunch[1][1][0], matrix_make_bunch[2][1][0], matrix_make_bunch[1][1][1], matrix_make_bunch[2][1][1]);   
	fprintf(F_cfg, "%f %f %f %f %f\n", matrix_make_bunch[0][2][0], matrix_make_bunch[1][2][0], matrix_make_bunch[2][2][0], matrix_make_bunch[1][2][1], matrix_make_bunch[2][2][1]);
	fprintf(F_cfg, "%f %f %f %f %f\n", matrix_make_bunch[0][3][0], matrix_make_bunch[1][3][0], matrix_make_bunch[2][3][0], matrix_make_bunch[1][3][1], matrix_make_bunch[2][3][1]);
	fprintf(F_cfg, "%f %f %f %f %f\n", matrix_make_bunch[0][4][0], matrix_make_bunch[1][4][0], matrix_make_bunch[2][4][0], matrix_make_bunch[1][4][1], matrix_make_bunch[2][4][1]);
	fprintf(F_cfg, "%f %f %f %f %f\n", matrix_make_bunch[0][5][0], matrix_make_bunch[1][5][0], matrix_make_bunch[2][5][0], matrix_make_bunch[1][5][1], matrix_make_bunch[2][5][1]);
	fprintf(F_cfg, "%f %f %f %f %f\n", matrix_make_bunch[0][6][0], matrix_make_bunch[1][6][0], matrix_make_bunch[2][6][0], matrix_make_bunch[1][6][1], matrix_make_bunch[2][6][1]);
	fprintf(F_cfg, "%f %f %f %f %f\n", matrix_make_bunch[0][7][0], matrix_make_bunch[1][7][0], matrix_make_bunch[2][7][0], matrix_make_bunch[1][7][1], matrix_make_bunch[2][7][1]);
	fprintf(F_cfg, "%f %f %f %f %f\n", matrix_make_bunch[0][8][0], matrix_make_bunch[1][8][0], matrix_make_bunch[2][8][0], matrix_make_bunch[1][8][1], matrix_make_bunch[2][8][1]);
	fprintf(F_cfg, "%f %f %f %f %f\n", matrix_make_bunch[0][9][0], matrix_make_bunch[1][9][0], matrix_make_bunch[2][9][0], matrix_make_bunch[1][9][1], matrix_make_bunch[2][9][1]);
	fprintf(F_cfg,"%d\n",yz);
	GetCtrlVal(panelM_B,PANEL_M_B_AUTO_REFRESH,&auto_refresh);
	fprintf(F_cfg,"%d\n",auto_refresh);
	fclose(F_cfg);
}

int CVICALLBACK Ring_bunch (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	int i;
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal (panelSet, PANEL_SET_RING_B_S, &par.cb);
			switch (par.cb)
			{
				case 0: //банч из простых уставок
					SetCtrlAttribute (panelSet, PANEL_SET_TUNING, ATTR_DIMMED, 1);
					SetCtrlAttribute (panelSet, PANEL_SET_SIGMA_R_1, ATTR_DIMMED, 0);
					SetCtrlAttribute (panelSet, PANEL_SET_SIGMA_L_1, ATTR_DIMMED, 0);
					SetCtrlAttribute (panelSet, PANEL_SET_BUNCH_Y_1, ATTR_DIMMED, 0);
				break;
				case 1: //банч из файла
					SetCtrlAttribute (panelSet, PANEL_SET_TUNING, ATTR_DIMMED, 1);
					SetCtrlAttribute (panelSet, PANEL_SET_SIGMA_R_1, ATTR_DIMMED, 0);
					SetCtrlAttribute (panelSet, PANEL_SET_SIGMA_L_1, ATTR_DIMMED, 0);
					SetCtrlAttribute (panelSet, PANEL_SET_BUNCH_Y_1, ATTR_DIMMED, 0);
				break;
				case 2: //банч из сложных уставок
					SetCtrlAttribute (panelSet, PANEL_SET_TUNING, ATTR_DIMMED, 0);
					SetCtrlAttribute (panelSet, PANEL_SET_SIGMA_R_1, ATTR_DIMMED, 1);
					SetCtrlAttribute (panelSet, PANEL_SET_SIGMA_L_1, ATTR_DIMMED, 1);
					SetCtrlAttribute (panelSet, PANEL_SET_BUNCH_Y_1, ATTR_DIMMED, 1);
					DisplayPanel (panelM_B);
					RecallPS (panelM_B, "make_bunch.cfg");
					Set_zero_color();
					Set_y_z();
					for(i=0;i<256;i++)
					{
						CT[i]=i+i*256+i*256*256;
					}
					SetCtrlAttribute (panelM_B, PANEL_M_B_GRAPH_2, ATTR_YLOOSE_FIT_AUTOSCALING, 1);
					Refresh_rho();	
					Refresh_sigma();
					Refresh_delta();
					Draw_bunch();
			}
	}
	return 0;
}

int CVICALLBACK Tuning (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			DisplayPanel (panelM_B);
			Refresh_rho();
			Refresh_sigma();
			Refresh_delta();
			Draw_bunch();
			break;
	}
	return 0;
}

int CVICALLBACK Ok_M_B (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			HidePanel (panelM_B);
			SavePS(panelM_B, "make_bunch.cfg");
			break;
	}
	return 0;
}

int CVICALLBACK Panel_make_bunch (int panel, int event, void *callbackData, int eventData1, int eventData2)
{
	int ctrl,data_type;
	double ini;
	switch (event)
	{
		case EVENT_LEFT_DOUBLE_CLICK:
			ctrl = GetActiveCtrl (panelM_B);
			GetCtrlAttribute (panelM_B, ctrl, ATTR_DATA_TYPE, &data_type);
			if (data_type==VAL_INTEGER)
			{
				return 0;
			}
			if ((change_switch==2)&&(ctrl==change_ctrl))
			{
				SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_INI,change_ini);
				SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_FIN,change_fin);
				SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_STEP,change_step);
				SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_NUMBER,change_number);
			}
			else
			{
				GetCtrlVal(panelM_B,ctrl,&ini);
				SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_INI,ini);
				SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_FIN,ini);
				SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_STEP,0.0);
				SetCtrlVal(panelCHNG,PANEL_CHNG_NUMERIC_NUMBER,0);
			}
			DisplayPanel (panelCHNG);
	}
	return 0;
}
