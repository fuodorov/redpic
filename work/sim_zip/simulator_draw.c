#include <cviauto.h>
#include "3DGraphCtrl.h"
#include <utility.h>
#include <userint.h>
#include <ansi_c.h>
#include "simulator_beam.h"
#include "simulator_vars.h"
#include "simulator_draw.h"

void make_image(struct particle *b,struct parametrs p, int check)
{
	int k,i,j,i0,j0, bitmap, bg_white, jet_map;
	unsigned char *bits;
	double *im, max=0, s, gamma;
	struct 
	{
	    double r,g,b;
	} c, cb;
	bits = calloc(4*height*width,sizeof(unsigned char));
	im =   calloc(  height*width,sizeof(double));
	GetCtrlVal(panelIM_SET,P_IM_SET_SIGMA,&s);
	GetCtrlVal(panelIM_SET,P_IM_SET_GAMMA,&gamma);
	GetCtrlVal(panelIM_SET,P_IM_SET_COLOR_R,&cb.r);
	GetCtrlVal(panelIM_SET,P_IM_SET_COLOR_G,&cb.g);
	GetCtrlVal(panelIM_SET,P_IM_SET_COLOR_B,&cb.b);
	for(k=0;k<p.t_number;k++)
	{
		if (b[k].status==on_screen || check == 0)  // check or not particles are only on screen
		{
			i0=width*(b[k].x*p.step/fl-p.x_min)/(p.x_max-p.x_min);
			j0=height-height*(b[k].y*p.step/fl-p.y_min)/(p.y_max-p.y_min);
			for (i=i0-4*s;i<i0+4*s;i++)
			{
				for (j=j0-4*s;j<j0+4*s;j++)
				{
					if ((i<width)&&(i>=0)&&(j<height)&&(j>=0))
					{
						im[i+j*width]+=exp(-((i-i0)*(i-i0)+(j-j0)*(j-j0))/(2*s*s));//G
						if (im[i+j*width]>max)	 max=im[i+j*width];
					}
				}
			}
		}
	}
	for (i=0;i<width*height;i++)
	{
		GetCtrlVal(panelIM_SET,P_IM_SET_GRID_3,&jet_map);
		if(jet_map)  // jet color map
		{
			if (im[i] < (0 + 0.25 * max)) 
			{
		      	c.r = 0;
      			c.g = 4 * (im[i] - 0) / max;
				c.b = 1;
   			} 
			else if (im[i] < (0 + 0.5 * max)) 
			{
		      	c.r = 0;
				c.g = 1;
      		  	c.b = 1 + 4 * (0 + 0.25 * max - im[i]) / max;
   			} 
			else if (im[i] < (0 + 0.75 * max)) 
			{
      		  	c.r = 4 * (im[i] - 0 - 0.5 * max) / max;
				c.g = 1;
      		  	c.b = 0;
   			} 
			else 
			{
	      		c.r = 1;
				c.g = 1 + 4 * (0 + 0.75 * max - im[i]) / max;
      			c.b = 0;
   			}
			bits[4*i+2] = 255*c.r; //Red
			bits[4*i+1] = 255*c.g; //Green
			bits[4*i+0] = 255*c.b; //Blue
		}
		else   // gray color map 
		{
			GetCtrlVal(panelIM_SET,P_IM_SET_GRID_2,&bg_white); 
			bits[4*i+2] = 255*bg_white - (2*bg_white-1)*255*pow(im[i]/max,gamma)*cb.r;//R
			bits[4*i+1] = 255*bg_white - (2*bg_white-1)*255*pow(im[i]/max,gamma)*cb.g;//G
			bits[4*i+0] = 255*bg_white - (2*bg_white-1)*255*pow(im[i]/max,gamma)*cb.b;//B
		}
	}
	NewBitmap (4*width*sizeof(unsigned char), 32, width, height, NULL, bits, NULL, &bitmap); free(bits); free(im);
	CanvasDrawBitmap (panelS, PANEL_S_CANVAS, bitmap, VAL_ENTIRE_OBJECT, MakeRect ( 0, 0, height, width));
}

void make_grid (struct parametrs p)
{
	int Px,Py,i, i_dx,i_dy;
	double gr_br, cx, cy, d_dx=1, d_dy=1;

	Px = floor(log10((p.x_max-p.x_min)/10)); i_dx=0.1*(p.x_max-p.x_min)/pow(10,Px);
	Py = floor(log10((p.y_max-p.y_min)/10)); i_dy=0.1*(p.y_max-p.y_min)/pow(10,Py);
	
	if ( i_dx<=1.5)				  d_dx =  1*pow(10,Px); if (i_dy<=1.5)				  d_dy =  1*pow(10,Py);
	if ((i_dx> 1.5)&&(i_dx<=3.5)) d_dx =  2*pow(10,Px); if ((i_dy> 1.5)&&(i_dy<=3.5)) d_dy =  2*pow(10,Py);
	if ((i_dx> 3.5)&&(i_dx<=7.5)) d_dx =  5*pow(10,Px); if ((i_dy> 3.5)&&(i_dy<=7.5)) d_dy =  5*pow(10,Py);
	if ( i_dx> 7.5)				  d_dx = 10*pow(10,Px); if (i_dy> 7.5)				  d_dy = 10*pow(10,Py);
	//d_dx = RoundRealToNearestInteger (  pow(2.154,RoundRealToNearestInteger (log(i_dx)/log(2.154)))  )*pow(10,Px);
	//d_dy = RoundRealToNearestInteger (  pow(2.154,RoundRealToNearestInteger (log(i_dy)/log(2.154)))  )*pow(10,Py);
									
	SetCtrlVal (panelS, PANEL_S_X_STEP, d_dx);
	SetCtrlVal (panelS, PANEL_S_Y_STEP, d_dy);
	GetCtrlVal(panelIM_SET,P_IM_SET_GRID_BR,&gr_br);
	SetCtrlAttribute (panelS, PANEL_S_CANVAS, ATTR_PEN_COLOR, (int)(0x19*gr_br)+(int)(0xaf*gr_br)*256+(int)(0x19*gr_br)*256*256);
	SetCtrlAttribute (panelS, PANEL_S_CANVAS, ATTR_PEN_STYLE, 2);
	SetCtrlAttribute (panelS, PANEL_S_CANVAS, ATTR_PEN_WIDTH, 1);
	
	cx = width /(p.x_max-p.x_min);
	cy = height/(p.y_max-p.y_min);
	
	for(i=0;i<((p.x_max-p.x_min)/d_dx);i++) 
	{ //  ==  draw grid lines  ==  //
		CanvasDrawLine (panelS, PANEL_S_CANVAS, 
			MakePoint (  (int)(cx*(i*d_dx+d_dx*ceil(p.x_min/d_dx)-p.x_min)), 0),
			  MakePoint ((int)(cx*(i*d_dx+d_dx*ceil(p.x_min/d_dx)-p.x_min)), height));
		CanvasDrawLine (panelS, PANEL_S_CANVAS, 
				MakePoint (0,height-(int)(cy*(i*d_dy+d_dy*ceil(p.y_min/d_dy)-p.y_min))), 
			MakePoint (width,height-(int)(cy*(i*d_dy+d_dy*ceil(p.y_min/d_dy)-p.y_min))));
	} //  ==  draw coordinate system  ==  //
	SetCtrlAttribute (panelS, PANEL_S_CANVAS, ATTR_PEN_COLOR, (int)(0xff*gr_br)+(int)(0xff*gr_br)*256+(int)(0xff*gr_br)*256*256);
	SetCtrlAttribute (panelS, PANEL_S_CANVAS, ATTR_PEN_STYLE, 2);
	SetCtrlAttribute (panelS, PANEL_S_CANVAS, ATTR_PEN_WIDTH, 1);
	CanvasDrawLine   (panelS, PANEL_S_CANVAS, MakePoint ((int)(cx*(-p.x_min)), 0),MakePoint ((int)(cx*(-p.x_min)), height));
	CanvasDrawLine   (panelS, PANEL_S_CANVAS, MakePoint (0,height-(int)(-cy*p.y_min)),MakePoint (width,height-(int)(-cy*p.y_min)));
}

void draw_3D (struct particle *b, struct parametrs p)
{
	CAObjHandle graphHandle, plotsHandle, plotHandle;
	int i;
	double *data1;
    double *data2;
    double *data3;
	int gAxis, axes, size;
	VARIANT xArray,yArray,zArray;
	
	GetObjHandleFromActiveXCtrl (panel_3D,PANEL_3D_GRAPH, &graphHandle);
    CW3DGraphLib__DCWGraph3DGetPlots (graphHandle, NULL, &plotsHandle);
    CW3DGraphLib_CWPlots3DAdd (plotsHandle, NULL,&plotHandle);
	CW3DGraphLib_CWPlots3DRemoveAll (plotsHandle, NULL);
	CW3DGraphLib__DCWGraph3DGetPlots (graphHandle, NULL, &plotsHandle);
    CW3DGraphLib_CWPlots3DAdd (plotsHandle, NULL,&plotHandle);
	CW3DGraphLib__DCWGraph3DGetAxes(graphHandle, NULL, &axes);
	CW3DGraphLib_CWAxes3DItem(axes, NULL, CA_VariantInt(0 + 1), &gAxis);
	CW3DGraphLib_CWAxis3DSetAutoScale (gAxis, NULL, VFALSE);
	CW3DGraphLib_CWAxis3DSetMaximum(gAxis, NULL, CA_VariantDouble( p.x_max*p.step/fl));
	CW3DGraphLib_CWAxis3DSetMinimum(gAxis, NULL, CA_VariantDouble( p.x_min*p.step/fl));
	CW3DGraphLib_CWAxis3DSetFormatString (gAxis, NULL, ".0");
	CW3DGraphLib_CWAxis3DSetInverted (gAxis, NULL, VTRUE);
	CW3DGraphLib_CWAxes3DItem(axes, NULL, CA_VariantInt(1 + 1), &gAxis);
	CW3DGraphLib_CWAxis3DSetAutoScale (gAxis, NULL, VFALSE);
	CW3DGraphLib_CWAxis3DSetMaximum(gAxis, NULL, CA_VariantDouble(0.1+p.screen*p.step/fl));
	CW3DGraphLib_CWAxis3DSetMinimum(gAxis, NULL, CA_VariantDouble(0.0));
	CW3DGraphLib_CWAxis3DSetFormatString (gAxis, NULL, ".");
	CW3DGraphLib_CWAxes3DItem(axes, NULL, CA_VariantInt(2 + 1), &gAxis);
	CW3DGraphLib_CWAxis3DSetAutoScale (gAxis, NULL, VFALSE);
	CW3DGraphLib_CWAxis3DSetMaximum(gAxis, NULL, CA_VariantDouble( p.y_max*p.step/fl));
	CW3DGraphLib_CWAxis3DSetMinimum(gAxis, NULL, CA_VariantDouble( p.y_min*p.step/fl));
	CW3DGraphLib_CWAxis3DSetFormatString (gAxis, NULL, ".0");
    
	size = p.t_number; //sizeof(b)/sizeof(struct particle);
	data1 = malloc(size*sizeof(double));
	data2 = malloc(size*sizeof(double));
	data3 = malloc(size*sizeof(double));
	for (i = 0; i<size; i++)
	{
		data1[i] = p.a_y*b[i].x*p.step/fl;
		data3[i] = p.a_y*b[i].y*p.step/fl;
		data2[i] = 		 b[i].z*p.step/fl;
	}
	CA_VariantSet1DArray (&xArray, CAVT_DOUBLE, size, data1);
	CA_VariantSet1DArray (&yArray, CAVT_DOUBLE, size, data2);
	CA_VariantSet1DArray (&zArray, CAVT_DOUBLE, size, data3);
	CW3DGraphLib_CWPlot3DSetStyle (plotHandle, NULL, CW3DGraphLibConst_cwPoint);
	//CW3DGraphLib_CWPlot3DSetFillStyle (plotHandle, NULL, CW3DGraphLibConst_cwFlat);
	CW3DGraphLib_CWPlot3DSetPointStyle (plotHandle, NULL, CW3DGraphLibConst_cwPoint3DSolidCircle);
	CW3DGraphLib_CWPlot3DSetPointSize (plotHandle, NULL, 2);
	CW3DGraphLib_CWPlot3DSetColorMapStyle (plotHandle, NULL, CW3DGraphLibConst_cwNone);
	CW3DGraphLib_CWPlot3DSetPointColor (plotHandle, NULL, VAL_GREEN);
    CW3DGraphLib_CWPlot3DPlot3DCurve (plotHandle, NULL, xArray, yArray, zArray, CA_DEFAULT_VAL); 
	CA_VariantClear (&xArray);
	CA_VariantClear (&yArray); 
	CA_VariantClear (&zArray);
	ProcessDrawEvents ();
	ProcessSystemEvents ();
}

int CVICALLBACK Refresh_image (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	int grid;
	switch (event)
	{
		case EVENT_VAL_CHANGED: //case EVENT_COMMIT:
		
			if (beam) 	make_image(beam,par,1); 
			else CanvasClear (panelS, PANEL_S_CANVAS, MakeRect ( 0, 0, height, width));
			GetCtrlVal(panelIM_SET,P_IM_SET_GRID,&grid);
			if (grid) 	make_grid(par);
	}
	return 0;
}

int CVICALLBACK Quit_IM_SET (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event) { case EVENT_COMMIT: HidePanel(panelIM_SET); }
	return 0;
}

int CVICALLBACK Plot_3D (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			GetCtrlVal (panelS, PANEL_S_PLOT_3D, &par.plot_3D);
			if (par.plot_3D) DisplayPanel (panel_3D); 
			else 			 HidePanel (panel_3D);
	}
	return 0;
}
