/**************************************************************************/
/* LabWindows/CVI User Interface Resource (UIR) Include File              */
/* Copyright (c) National Instruments 2005. All Rights Reserved.          */
/*                                                                        */
/* WARNING: Do not add to, delete from, or otherwise modify the contents  */
/*          of this include file.                                         */
/**************************************************************************/

#include <userint.h>

#ifdef __cplusplus
    extern "C" {
#endif

     /* Panels and Controls: */

#define  P_IM_SET                        1
#define  P_IM_SET_COLOR_R                2       /* callback function: Set_gamma */
#define  P_IM_SET_COLOR_G                3       /* callback function: Set_gamma */
#define  P_IM_SET_COLOR_B                4       /* callback function: Set_gamma */
#define  P_IM_SET_GAMMA                  5       /* callback function: Set_gamma */
#define  P_IM_SET_GRID_BR                6       /* callback function: Set_grid_brightness */
#define  P_IM_SET_QUITBUTTON             7       /* callback function: Quit_IM_SET */
#define  P_IM_SET_DECORATION             8
#define  P_IM_SET_TEXTMSG_2              9
#define  P_IM_SET_TEXTMSG                10
#define  P_IM_SET_AREAY                  11      /* callback function: Areay */
#define  P_IM_SET_AREAX                  12      /* callback function: Areax */
#define  P_IM_SET_GRID                   13      /* callback function: Make_grid */
#define  P_IM_SET_AUTO_SCALE             14
#define  P_IM_SET_BOTH                   15
#define  P_IM_SET_DECORATION_2           16

#define  PANEL_B                         2
#define  PANEL_B_BUNCH_N                 2
#define  PANEL_B_SIGMA_L                 3
#define  PANEL_B_BUNCH_Y                 4
#define  PANEL_B_SIGMA_R                 5
#define  PANEL_B_ENERGY_B                6
#define  PANEL_B_OKBUTTON                7       /* callback function: OkCallback_b */
#define  PANEL_B_RING_B                  8

#define  PANEL_G                         3
#define  PANEL_G_RADIUS                  2
#define  PANEL_G_RF_LENGTH               3
#define  PANEL_G_LENS_L                  4
#define  PANEL_G_SCAN_W                  5
#define  PANEL_G_SCAN_L                  6
#define  PANEL_G_SCAN                    7
#define  PANEL_G_BUNCH_W                 8
#define  PANEL_G_RF                      9
#define  PANEL_G_BUNCH                   10
#define  PANEL_G_LENS1                   11
#define  PANEL_G_SCREEN                  12
#define  PANEL_G_OKBUTTON                13      /* callback function: OkCallback_g */
#define  PANEL_G_COMMANDBUTTON           14      /* callback function: test */

#define  PANEL_O                         4
#define  PANEL_O_OKBUTTON                2       /* callback function: OkCallback_o */
#define  PANEL_O_TIME_S                  3
#define  PANEL_O_U                       4
#define  PANEL_O_B_MAX                   5
#define  PANEL_O_STEP                    6

#define  PANEL_S                         5
#define  PANEL_S_CANVAS                  2       /* callback function: Canvas_callback */
#define  PANEL_S_QUITBUTTON              3       /* callback function: QuitCallback */
#define  PANEL_S_BEAM                    4       /* callback function: Beam */
#define  PANEL_S_BUNCH                   5       /* callback function: Bunch */
#define  PANEL_S_START                   6       /* callback function: Start */
#define  PANEL_S_Y_MAX                   7
#define  PANEL_S_Y_MIN                   8
#define  PANEL_S_X_MAX                   9
#define  PANEL_S_X_MIN                   10
#define  PANEL_S_GEO                     11      /* callback function: Geo */
#define  PANEL_S_OPTION                  12      /* callback function: Option */
#define  PANEL_S_TEXTBOX                 13
#define  PANEL_S_IND                     14
#define  PANEL_S_REDRAW                  15      /* callback function: Redraw */

#define  PANEL_T                         6
#define  PANEL_T_ENERGY_T                2
#define  PANEL_T_TEMP                    3
#define  PANEL_T_NUM                     4
#define  PANEL_T_LENGTH                  5
#define  PANEL_T_SIZE                    6
#define  PANEL_T_OKBUTTON                7       /* callback function: OkCallback_t */
#define  PANEL_T_RING_T                  8


     /* Menu Bars, Menus, and Menu Items: */

          /* (no menu bars in the resource file) */


     /* Callback Prototypes: */ 

int  CVICALLBACK Areax(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Areay(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Beam(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Bunch(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Canvas_callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Geo(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Make_grid(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OkCallback_b(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OkCallback_g(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OkCallback_o(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OkCallback_t(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Option(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Quit_IM_SET(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK QuitCallback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Redraw(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Set_gamma(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Set_grid_brightness(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Start(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK test(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);


#ifdef __cplusplus
    }
#endif
