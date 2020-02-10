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

#define  PANEL                           1
#define  PANEL_A_R                       2       /* callback function: Refresh_Rho */
#define  PANEL_B_R                       3       /* callback function: Refresh_Rho */
#define  PANEL_C_R                       4       /* callback function: Refresh_Rho */
#define  PANEL_D_R                       5       /* callback function: Refresh_Rho */
#define  PANEL_O_R                       6       /* callback function: Refresh_Rho */
#define  PANEL_Fi_R                      7       /* callback function: Refresh_Rho */
#define  PANEL_E_R                       8       /* callback function: Refresh_Rho */
#define  PANEL_X0_R                      9       /* callback function: Refresh_Rho */
#define  PANEL_S_R                       10      /* callback function: Refresh_Rho */
#define  PANEL_F_R                       11      /* callback function: Refresh_Rho */
#define  PANEL_A_S                       12      /* callback function: Refresh_Sigma */
#define  PANEL_B_S                       13      /* callback function: Refresh_Sigma */
#define  PANEL_C_S                       14      /* callback function: Refresh_Sigma */
#define  PANEL_D_S                       15      /* callback function: Refresh_Sigma */
#define  PANEL_O_S                       16      /* callback function: Refresh_Sigma */
#define  PANEL_Fi_S                      17      /* callback function: Refresh_Sigma */
#define  PANEL_E_S                       18      /* callback function: Refresh_Sigma */
#define  PANEL_X0_S                      19      /* callback function: Refresh_Sigma */
#define  PANEL_S_S                       20      /* callback function: Refresh_Sigma */
#define  PANEL_F_S                       21      /* callback function: Refresh_Sigma */
#define  PANEL_A_D                       22      /* callback function: Refresh_Delta */
#define  PANEL_B_D                       23      /* callback function: Refresh_Delta */
#define  PANEL_C_D                       24      /* callback function: Refresh_Delta */
#define  PANEL_D_D                       25      /* callback function: Refresh_Delta */
#define  PANEL_O_D                       26      /* callback function: Refresh_Delta */
#define  PANEL_Fi_D                      27      /* callback function: Refresh_Delta */
#define  PANEL_E_D                       28      /* callback function: Refresh_Delta */
#define  PANEL_X0_D                      29      /* callback function: Refresh_Delta */
#define  PANEL_S_D                       30      /* callback function: Refresh_Delta */
#define  PANEL_F_D                       31      /* callback function: Refresh_Delta */
#define  PANEL_XMIN                      32      /* callback function: Refresh */
#define  PANEL_XMAX                      33      /* callback function: Refresh */
#define  PANEL_NUM                       34      /* callback function: Refresh */
#define  PANEL_SAVE_FILE                 35      /* callback function: Save_file */
#define  PANEL_REFRESH_                  36      /* callback function: Refresh_M_B */
#define  PANEL_EXIT                      37      /* callback function: Exit */
#define  PANEL_GRAPH_2                   38
#define  PANEL_GRAPH                     39
#define  PANEL_CANVAS                    40
#define  PANEL_AUTO_REFRESH              41      /* callback function: Auto_refresh */
#define  PANEL_Y_Z                       42      /* callback function: y_z */
#define  PANEL_TEXTMSG                   43
#define  PANEL_TEXTMSG_4                 44
#define  PANEL_TEXTMSG_3                 45
#define  PANEL_TEXTMSG_8                 46
#define  PANEL_TEXTMSG_9                 47
#define  PANEL_TEXTMSG_7                 48
#define  PANEL_TEXTMSG_10                49
#define  PANEL_TEXTMSG_11                50
#define  PANEL_TEXTMSG_12                51
#define  PANEL_TEXTMSG_13                52
#define  PANEL_TEXTMSG_6                 53
#define  PANEL_TEXTMSG_14                54
#define  PANEL_TEXTMSG_15                55
#define  PANEL_TEXTMSG_16                56
#define  PANEL_TEXTMSG_17                57
#define  PANEL_TEXTMSG_2                 58
#define  PANEL_TEXTMSG_19                59
#define  PANEL_TEXTMSG_20                60
#define  PANEL_TEXTMSG_21                61
#define  PANEL_TEXTMSG_22                62
#define  PANEL_TEXTMSG_23                63
#define  PANEL_TEXTMSG_24                64
#define  PANEL_TEXTMSG_25                65
#define  PANEL_TEXTMSG_26                66
#define  PANEL_TEXTMSG_27                67
#define  PANEL_TEXTMSG_28                68
#define  PANEL_TEXTMSG_29                69
#define  PANEL_TEXTMSG_30                70
#define  PANEL_TEXTMSG_31                71
#define  PANEL_TEXTMSG_32                72
#define  PANEL_TEXTMSG_33                73
#define  PANEL_TEXTMSG_34                74
#define  PANEL_TEXTMSG_36                75
#define  PANEL_TEXTMSG_37                76
#define  PANEL_TEXTMSG_38                77
#define  PANEL_TEXTMSG_39                78
#define  PANEL_TEXTMSG_40                79
#define  PANEL_TEXTMSG_41                80
#define  PANEL_TEXTMSG_42                81
#define  PANEL_TEXTMSG_43                82
#define  PANEL_TEXTMSG_44                83
#define  PANEL_TEXTMSG_45                84
#define  PANEL_TEXTMSG_46                85
#define  PANEL_TEXTMSG_47                86
#define  PANEL_TEXTMSG_48                87
#define  PANEL_TEXTMSG_49                88
#define  PANEL_TEXTMSG_50                89
#define  PANEL_TEXTMSG_51                90
#define  PANEL_TEXTMSG_5                 91
#define  PANEL_TEXTMSG_18                92
#define  PANEL_TEXTMSG_60                93
#define  PANEL_TEXTMSG_61                94
#define  PANEL_TEXTMSG_SIGMA             95
#define  PANEL_TEXTMSG_DELTA             96
#define  PANEL_TEXTMSG_35                97


     /* Menu Bars, Menus, and Menu Items: */

          /* (no menu bars in the resource file) */


     /* Callback Prototypes: */ 

int  CVICALLBACK Auto_refresh(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Exit(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Refresh(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Refresh_Delta(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Refresh_M_B(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Refresh_Rho(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Refresh_Sigma(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Save_file(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK y_z(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);


#ifdef __cplusplus
    }
#endif
