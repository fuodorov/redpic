/**************************************************************************/
/* LabWindows/CVI User Interface Resource (UIR) Include File              */
/* Copyright (c) National Instruments 2016. All Rights Reserved.          */
/*                                                                        */
/* WARNING: Do not add to, delete from, or otherwise modify the contents  */
/*          of this include file.                                         */
/**************************************************************************/

#include <userint.h>

#ifdef __cplusplus
    extern "C" {
#endif

     /* Panels and Controls: */

#define  P_IM_SET                         1
#define  P_IM_SET_COLOR_R                 2       /* callback function: Refresh_image */
#define  P_IM_SET_COLOR_G                 3       /* callback function: Refresh_image */
#define  P_IM_SET_COLOR_B                 4       /* callback function: Refresh_image */
#define  P_IM_SET_GAMMA                   5       /* callback function: Refresh_image */
#define  P_IM_SET_GRID_BR                 6       /* callback function: Refresh_image */
#define  P_IM_SET_QUITBUTTON              7       /* callback function: Quit_IM_SET */
#define  P_IM_SET_DECORATION              8
#define  P_IM_SET_TEXTMSG_2               9
#define  P_IM_SET_TEXTMSG                 10
#define  P_IM_SET_SIGMA                   11      /* callback function: Refresh_image */
#define  P_IM_SET_GRID_2                  12      /* callback function: Refresh_image */
#define  P_IM_SET_GRID_3                  13      /* callback function: Refresh_image */
#define  P_IM_SET_GRID                    14      /* callback function: Refresh_image */

#define  PANEL                            2
#define  PANEL_RING_Q                     2
#define  PANEL_QUADS_AMOUNT               3
#define  PANEL_QUAD                       4
#define  PANEL_QUAD_L_2                   5
#define  PANEL_QUAD_L                     6

#define  PANEL_3D                         3
#define  PANEL_3D_GRAPH                   2

#define  PANEL_A                          4
#define  PANEL_A_QUITBUTTON               2       /* callback function: QuitCallback_A */
#define  PANEL_A_GRAPH                    3
#define  PANEL_A_ANALYSE_S                4       /* callback function: Save_a */
#define  PANEL_A_POLYFIT                  5       /* callback function: Polyfit */
#define  PANEL_A_TEST                     6       /* callback function: Test */
#define  PANEL_A_ANALYSE                  7       /* callback function: Analyse */
#define  PANEL_A_XMAX                     8
#define  PANEL_A_XMIN                     9
#define  PANEL_A_N                        10
#define  PANEL_A_D                        11
#define  PANEL_A_K                        12
#define  PANEL_A_ORDER                    13

#define  PANEL_CHNG                       5
#define  PANEL_CHNG_NUMERIC_NUMBER        2       /* callback function: Set_number */
#define  PANEL_CHNG_NUMERIC_STEP          3       /* callback function: Set_step */
#define  PANEL_CHNG_NUMERIC_FIN           4       /* callback function: Set_range */
#define  PANEL_CHNG_NUMERIC_INI           5       /* callback function: Set_range */
#define  PANEL_CHNG_OKBUTTON              6       /* callback function: OkCallbackCHNG */

#define  PANEL_G                          6
#define  PANEL_G_GRAPH                    2
#define  PANEL_G_QUITBUTTON               3       /* callback function: QuitCallback_graph */
#define  PANEL_G_RING_FIELDS              4       /* callback function: Plot_graph */
#define  PANEL_G_FIELD_X                  5       /* callback function: Plot_graph */
#define  PANEL_G_FIELD_Y                  6       /* callback function: Plot_graph */

#define  PANEL_M_B                        7       /* callback function: Panel_make_bunch */
#define  PANEL_M_B_A_R                    2       /* callback function: Refresh_Rho */
#define  PANEL_M_B_B_R                    3       /* callback function: Refresh_Rho */
#define  PANEL_M_B_C_R                    4       /* callback function: Refresh_Rho */
#define  PANEL_M_B_D_R                    5       /* callback function: Refresh_Rho */
#define  PANEL_M_B_O_R                    6       /* callback function: Refresh_Rho */
#define  PANEL_M_B_Fi_R                   7       /* callback function: Refresh_Rho */
#define  PANEL_M_B_E_R                    8       /* callback function: Refresh_Rho */
#define  PANEL_M_B_X0_R                   9       /* callback function: Refresh_Rho */
#define  PANEL_M_B_S_R                    10      /* callback function: Refresh_Rho */
#define  PANEL_M_B_F_R                    11      /* callback function: Refresh_Rho */
#define  PANEL_M_B_A_S                    12      /* callback function: Refresh_Sigma */
#define  PANEL_M_B_B_S                    13      /* callback function: Refresh_Sigma */
#define  PANEL_M_B_C_S                    14      /* callback function: Refresh_Sigma */
#define  PANEL_M_B_D_S                    15      /* callback function: Refresh_Sigma */
#define  PANEL_M_B_O_S                    16      /* callback function: Refresh_Sigma */
#define  PANEL_M_B_Fi_S                   17      /* callback function: Refresh_Sigma */
#define  PANEL_M_B_E_S                    18      /* callback function: Refresh_Sigma */
#define  PANEL_M_B_X0_S                   19      /* callback function: Refresh_Sigma */
#define  PANEL_M_B_S_S                    20      /* callback function: Refresh_Sigma */
#define  PANEL_M_B_F_S                    21      /* callback function: Refresh_Sigma */
#define  PANEL_M_B_A_D                    22      /* callback function: Refresh_Delta */
#define  PANEL_M_B_B_D                    23      /* callback function: Refresh_Delta */
#define  PANEL_M_B_C_D                    24      /* callback function: Refresh_Delta */
#define  PANEL_M_B_D_D                    25      /* callback function: Refresh_Delta */
#define  PANEL_M_B_O_D                    26      /* callback function: Refresh_Delta */
#define  PANEL_M_B_Fi_D                   27      /* callback function: Refresh_Delta */
#define  PANEL_M_B_E_D                    28      /* callback function: Refresh_Delta */
#define  PANEL_M_B_X0_D                   29      /* callback function: Refresh_Delta */
#define  PANEL_M_B_S_D                    30      /* callback function: Refresh_Delta */
#define  PANEL_M_B_F_D                    31      /* callback function: Refresh_Delta */
#define  PANEL_M_B_XMIN                   32      /* callback function: RefreshM_B */
#define  PANEL_M_B_XMAX                   33      /* callback function: RefreshM_B */
#define  PANEL_M_B_NUM                    34      /* callback function: RefreshM_B */
#define  PANEL_M_B_CANVAS                 35
#define  PANEL_M_B_GRAPH_2                36
#define  PANEL_M_B_GRAPH                  37
#define  PANEL_M_B_REFRESH_M_B            38      /* callback function: RefreshM_B */
#define  PANEL_M_B_AUTO_REFRESH           39      /* callback function: Auto_refresh */
#define  PANEL_M_B_TEXTMSG                40
#define  PANEL_M_B_TEXTMSG_4              41
#define  PANEL_M_B_TEXTMSG_3              42
#define  PANEL_M_B_TEXTMSG_8              43
#define  PANEL_M_B_TEXTMSG_9              44
#define  PANEL_M_B_TEXTMSG_7              45
#define  PANEL_M_B_TEXTMSG_10             46
#define  PANEL_M_B_TEXTMSG_11             47
#define  PANEL_M_B_TEXTMSG_12             48
#define  PANEL_M_B_TEXTMSG_13             49
#define  PANEL_M_B_TEXTMSG_6              50
#define  PANEL_M_B_TEXTMSG_14             51
#define  PANEL_M_B_TEXTMSG_15             52
#define  PANEL_M_B_TEXTMSG_16             53
#define  PANEL_M_B_TEXTMSG_17             54
#define  PANEL_M_B_TEXTMSG_2              55
#define  PANEL_M_B_TEXTMSG_19             56
#define  PANEL_M_B_TEXTMSG_20             57
#define  PANEL_M_B_TEXTMSG_21             58
#define  PANEL_M_B_TEXTMSG_22             59
#define  PANEL_M_B_TEXTMSG_23             60
#define  PANEL_M_B_TEXTMSG_24             61
#define  PANEL_M_B_TEXTMSG_25             62
#define  PANEL_M_B_TEXTMSG_26             63
#define  PANEL_M_B_TEXTMSG_27             64
#define  PANEL_M_B_TEXTMSG_28             65
#define  PANEL_M_B_TEXTMSG_29             66
#define  PANEL_M_B_TEXTMSG_30             67
#define  PANEL_M_B_TEXTMSG_31             68
#define  PANEL_M_B_TEXTMSG_32             69
#define  PANEL_M_B_TEXTMSG_33             70
#define  PANEL_M_B_TEXTMSG_54             71
#define  PANEL_M_B_TEXTMSG_53             72
#define  PANEL_M_B_TEXTMSG_52             73
#define  PANEL_M_B_TEXTMSG_34             74
#define  PANEL_M_B_TEXTMSG_36             75
#define  PANEL_M_B_TEXTMSG_37             76
#define  PANEL_M_B_TEXTMSG_38             77
#define  PANEL_M_B_TEXTMSG_39             78
#define  PANEL_M_B_TEXTMSG_40             79
#define  PANEL_M_B_TEXTMSG_41             80
#define  PANEL_M_B_TEXTMSG_42             81
#define  PANEL_M_B_TEXTMSG_51             82
#define  PANEL_M_B_TEXTMSG_43             83
#define  PANEL_M_B_TEXTMSG_5              84
#define  PANEL_M_B_TEXTMSG_44             85
#define  PANEL_M_B_TEXTMSG_18             86
#define  PANEL_M_B_TEXTMSG_45             87
#define  PANEL_M_B_TEXTMSG_46             88
#define  PANEL_M_B_TEXTMSG_47             89
#define  PANEL_M_B_TEXTMSG_SIGMA          90
#define  PANEL_M_B_TEXTMSG_48             91
#define  PANEL_M_B_TEXTMSG_DELTA          92
#define  PANEL_M_B_TEXTMSG_49             93
#define  PANEL_M_B_TEXTMSG_35             94
#define  PANEL_M_B_TEXTMSG_50             95
#define  PANEL_M_B_Ok                     96      /* callback function: Ok_M_B */
#define  PANEL_M_B_Y_Z                    97      /* callback function: y_z */

#define  PANEL_PHSP                       8
#define  PANEL_PHSP_QUITBUTTON            2       /* callback function: QuitCallback_PHSP */
#define  PANEL_PHSP_GRAPH_Z               3
#define  PANEL_PHSP_GRAPH_Y               4
#define  PANEL_PHSP_GRAPH_X               5

#define  PANEL_S                          9       /* callback function: Refresh */
#define  PANEL_S_CANVAS                   2       /* callback function: Canvas_callback */
#define  PANEL_S_QUITBUTTON               3       /* callback function: QuitCallback */
#define  PANEL_S_PHASE_SPACE              4       /* callback function: Phase_space */
#define  PANEL_S_ANALYSIS                 5       /* callback function: Analysis */
#define  PANEL_S_SAVE                     6       /* callback function: Save */
#define  PANEL_S_START                    7       /* callback function: Start */
#define  PANEL_S_Y_MAX                    8
#define  PANEL_S_Y_STEP                   9
#define  PANEL_S_Y_MIN                    10
#define  PANEL_S_X_MAX                    11
#define  PANEL_S_TIKALKA                  12
#define  PANEL_S_X_STEP                   13
#define  PANEL_S_X_MIN                    14
#define  PANEL_S_TEXTBOX                  15
#define  PANEL_S_REDRAW                   16      /* callback function: Redraw */
#define  PANEL_S_SETTINGS                 17      /* callback function: Settings */
#define  PANEL_S_CANVAS_2                 18
#define  PANEL_S_GRAPH_MODE               19      /* callback function: Graph_mode */
#define  PANEL_S_A_Y                      20      /* callback function: Graph_mode */
#define  PANEL_S_PLOT_2D                  21
#define  PANEL_S_PLOT_3D                  22      /* callback function: Plot_3D */
#define  PANEL_S_CLEAR_GRAPH              23
#define  PANEL_S_SV                       24
#define  PANEL_S_NP                       25

#define  PANEL_SET                        10      /* callback function: Panel_set */
#define  PANEL_SET_LENGTH_G               2
#define  PANEL_SET_SIZE                   3
#define  PANEL_SET_LENGTH                 4
#define  PANEL_SET_TEMP                   5
#define  PANEL_SET_WEIGHT_T               6
#define  PANEL_SET_NUM                    7
#define  PANEL_SET_CURRENT                8
#define  PANEL_SET_SPREAD                 9
#define  PANEL_SET_ENERGY_T               10
#define  PANEL_SET_RING_T                 11
#define  PANEL_SET_TEXTMSG                12
#define  PANEL_SET_RING_B_LN_D            13
#define  PANEL_SET_RING_B_TRD_2           14
#define  PANEL_SET_RING_B_TRD_1           15      /* callback function: Bunch_transverse_ring */
#define  PANEL_SET_BUNCH_TYPE             16
#define  PANEL_SET_BUNCH_N_2              17      /* callback function: Charge_2 */
#define  PANEL_SET_BUNCH_N                18      /* callback function: Charge_1 */
#define  PANEL_SET_BUNCH_X_2              19
#define  PANEL_SET_BUNCH_Z_2              20
#define  PANEL_SET_BUNCH_Y_2              21
#define  PANEL_SET_SIGMA_R_2              22
#define  PANEL_SET_BUNCH_ENERGY           23
#define  PANEL_SET_QUAD_X_4               24
#define  PANEL_SET_QUAD_Y_4               25
#define  PANEL_SET_QUAD_X_3               26
#define  PANEL_SET_QUAD_Y_3               27
#define  PANEL_SET_QUAD_X_1               28
#define  PANEL_SET_QUAD_Y_1               29
#define  PANEL_SET_BUNCH_X_1              30
#define  PANEL_SET_BUNCH_Z_1              31
#define  PANEL_SET_BUNCH_Y_1              32
#define  PANEL_SET_SIGMA_R_1              33
#define  PANEL_SET_SIGMA_L_2              34
#define  PANEL_SET_SIGMA_L_1              35
#define  PANEL_SET_TEXTMSG_2              36
#define  PANEL_SET_RADIUS                 37
#define  PANEL_SET_QUAD_Z_4               38
#define  PANEL_SET_QUAD_L_4               39
#define  PANEL_SET_QUAD_Z_3               40
#define  PANEL_SET_QUAD_L_3               41
#define  PANEL_SET_QUAD_Z_2               42
#define  PANEL_SET_QUAD_L_2               43
#define  PANEL_SET_QUAD_Z_1               44
#define  PANEL_SET_QUAD_L_1               45
#define  PANEL_SET_LENS                   46
#define  PANEL_SET_LENS_L                 47
#define  PANEL_SET_CR                     48
#define  PANEL_SET_CR_L                   49
#define  PANEL_SET_RF                     50
#define  PANEL_SET_RF_L                   51
#define  PANEL_SET_SCAN                   52
#define  PANEL_SET_SCAN_L                 53
#define  PANEL_SET_SCAN_W                 54
#define  PANEL_SET_ACC                    55
#define  PANEL_SET_ACC_G                  56
#define  PANEL_SET_ACC_L                  57
#define  PANEL_SET_BUNCH                  58
#define  PANEL_SET_BUNCH_RADIUS           59
#define  PANEL_SET_QUAD_G_4               60
#define  PANEL_SET_STOP                   61
#define  PANEL_SET_SCREEN                 62
#define  PANEL_SET_DIAPH                  63
#define  PANEL_SET_QUAD_G_3               64
#define  PANEL_SET_DIAPH_R                65
#define  PANEL_SET_CR_BX                  66
#define  PANEL_SET_CR_BY                  67
#define  PANEL_SET_STEP                   68
#define  PANEL_SET_QUAD_PHI_4             69
#define  PANEL_SET_QUAD_PHI_2             70
#define  PANEL_SET_QUAD_G_2               71
#define  PANEL_SET_U_RF                   72
#define  PANEL_SET_QUAD_PHI_3             73
#define  PANEL_SET_LAMBDA                 74
#define  PANEL_SET_QUAD_PHI_1             75
#define  PANEL_SET_ACC_PERIOD             76
#define  PANEL_SET_PERIOD                 77
#define  PANEL_SET_QUAD_G_1               78
#define  PANEL_SET_B_MAX                  79
#define  PANEL_SET_UX_SCAN                80
#define  PANEL_SET_UY_SCAN                81
#define  PANEL_SET_TIME_S                 82
#define  PANEL_SET_TEXTMSG_3              83
#define  PANEL_SET_RING_QUAD_4            84
#define  PANEL_SET_RING_RF                85
#define  PANEL_SET_RING_QUAD_3            86
#define  PANEL_SET_RING_QUAD_2            87
#define  PANEL_SET_RING_QUAD_1            88
#define  PANEL_SET_RING_SCAN              89
#define  PANEL_SET_PLOT_GRAPH_3           90      /* callback function: Load_config2 */
#define  PANEL_SET_PLOT_GRAPH_2           91      /* callback function: Save_config */
#define  PANEL_SET_PLOT_GRAPH             92      /* callback function: Plot_graph */
#define  PANEL_SET_OKBUTTON               93      /* callback function: OkCallback */
#define  PANEL_SET_SPLITTER               94
#define  PANEL_SET_TUNING                 95      /* callback function: Tuning */
#define  PANEL_SET_RING_T_S               96      /* callback function: Ring_beam */
#define  PANEL_SET_RING_B_S               97      /* callback function: Ring_bunch */
#define  PANEL_SET_SPACE_CHARGE           98
#define  PANEL_SET_SPLITTER_2             99
#define  PANEL_SET_NUMERICSLIDE_1         100     /* callback function: Slide */
#define  PANEL_SET_NUMERICSLIDE_2         101     /* callback function: Slide */
#define  PANEL_SET_DECORATION             102
#define  PANEL_SET_TEXTMSG_4              103
#define  PANEL_SET_N_FLAT                 104


     /* Menu Bars, Menus, and Menu Items: */

          /* (no menu bars in the resource file) */


     /* Callback Prototypes: */

int  CVICALLBACK Analyse(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Analysis(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Auto_refresh(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Bunch_transverse_ring(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Canvas_callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Charge_1(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Charge_2(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Graph_mode(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Load_config2(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Ok_M_B(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OkCallback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OkCallbackCHNG(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Panel_make_bunch(int panel, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Panel_set(int panel, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Phase_space(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Plot_3D(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Plot_graph(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Polyfit(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Quit_IM_SET(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK QuitCallback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK QuitCallback_A(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK QuitCallback_graph(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK QuitCallback_PHSP(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Redraw(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Refresh(int panel, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Refresh_Delta(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Refresh_image(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Refresh_Rho(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Refresh_Sigma(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK RefreshM_B(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Ring_beam(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Ring_bunch(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Save(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Save_a(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Save_config(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Set_number(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Set_range(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Set_step(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Settings(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Slide(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Start(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Test(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Tuning(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK y_z(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);


#ifdef __cplusplus
    }
#endif
