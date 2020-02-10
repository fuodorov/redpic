#include <ansi_c.h>
#include <userint.h>
#include "make_bunch.h"
#include "make_bunch_my.h"

int CVICALLBACK Save_file (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
int status_file;
char file_name[MAX_PATHNAME_LEN];
FILE* fo;
	switch (event) {
		case EVENT_COMMIT:
			Refresh_rho();
			Refresh_sigma();
			Refresh_delta();
			status_file = FileSelectPopup ("Save", "*.dat", "", "", VAL_SAVE_BUTTON,0, 0, 1, 1, file_name);
			fo=fopen (file_name, "w");
			fprintf(fo,"%d\n",NUM);
			for(i=0;i<NUM;i++)
			{
				fprintf(fo,"%f  %f  %f  %f  %f  %f\n",X[i],RHO[i],SIGMA[0][i],DELTA[0][i],SIGMA[1][i],DELTA[1][i]);
			}
			fclose(fo);

			break;
	}
	return 0;
}
