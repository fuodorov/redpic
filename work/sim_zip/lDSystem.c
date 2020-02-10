#include <ansi_c.h>
#include "lDSystem.h"


static double xyz[3][10000];
static double field[4][10000];

const int m_nCommentStrCount = 7;
const int m_nElementCount = 1851;
const int m_nCountElementsInOrder = 39;
const int m_nCountCommentInitial = 8;
const int m_nCountCommentsStrInTable = 3;
const int m_nNodesCount = 8;
const int m_nElementNumberToRead = 1;
const double dz = 0.009;
const double normU = 0.0001;
const char *m_sGeometryFile = "data\\geometry.dat";
const char *m_sFieldFile = "data\\efield.dat";

void DSReadLn (FILE *stream, char str[], int max)
{
	char ch;
	int i;
	
	for (i=0; i<max; i++) str[i] = 32;
	/* Read in first 80 characters and place them in "buffer": */
	ch = fgetc( stream );
	for( i=0; (i < max ) && ( feof( stream ) == 0 ) && (ch != 10); i++ )
	{
		str[i] = (char)ch;
		ch = fgetc( stream );
	}
}

void DSReadGeometryStr (FILE *stream, int *num, double *gx, double *gy, double *gz)
{
	char ch, snum[5]={0}, sx[13]={0}, sy[13]={0}, sz[13]={0};
	int i=0;
	
	ch = fgetc( stream );
	while ((ch == 32) || (ch == 10)) ch = fgetc( stream );
	while (ch != 32)
	{
		snum[i] = ch;
		i++;
		ch = fgetc( stream );
	}
	*num = atoi(snum);
	fread( sx, sizeof( char ), 12, stream );
	fread( sy, sizeof( char ), 12, stream );
	fread( sz, sizeof( char ), 12, stream );
	
	*gx = atof(sx);
	*gy = atof(sy);
	*gz = atof(sz);
	ch = fgetc( stream ); // Чтение перевода строки ASCI #10
}

void DSReadFieldStr (FILE *stream, int *num, double *gx, double *gy, double *gz, double *gsum)
{
	char ch, snum[5]={0}, sx[13]={0}, sy[13]={0}, sz[13]={0}, ssum[13]={0};
	int i=0;
	
	ch = fgetc( stream );
	while ((ch == 32) || (ch == 10)) ch = fgetc( stream );
	while (ch != 32)
	{
		snum[i] = ch;
		i++;
		ch = fgetc( stream );
	}
	*num = atoi(snum);
	fread( sx, sizeof( char ), 12, stream );
	fread( sy, sizeof( char ), 12, stream );
	fread( sz, sizeof( char ), 12, stream );
	fread( ssum, sizeof( char ), 12, stream );
	*gx = atof(sx);
	*gy = atof(sy);
	*gz = atof(sz);
	*gsum = atof(ssum);
	ch = fgetc( stream ); // Чтение перевода строки ASCI #10
}

int DSReadAnsysData()
{	
	FILE *stream;
	char str[100];
	int i, j, num=0, cnt=0;
	double node[4][100]={0};

	// Чтение файла геометрии
	if( (stream = fopen( m_sGeometryFile, "r-b" )) != NULL )
	{
		DSReadLn(stream, str, 100);
		DSReadLn(stream, str, 100);
		while ((feof( stream ) == 0) && (cnt < m_nElementCount))
		{
			for (i=0; i<m_nCommentStrCount; i++) DSReadLn(stream, str, 100);
			for (i=0; (i<m_nCountElementsInOrder) && (cnt < m_nElementCount); i++)
			{
				DSReadGeometryStr(stream, &num,&xyz[0][cnt],&xyz[1][cnt],&xyz[2][cnt]);
				cnt++;
			}
		}
		fclose(stream);
	}
	else return -1;	// Exit with error code (-1)

	cnt = 0;
	// Чтение файла данных
	if( (stream = fopen( m_sFieldFile, "r-b" )) != NULL )
	{
		DSReadLn(stream, str, 100);
		DSReadLn(stream, str, 100);
		while ((feof( stream ) == 0) && (cnt < m_nElementCount))
		{
			for (i=0; i<m_nCountCommentInitial; i++) DSReadLn(stream, str, 100);
			for (j=0; (j<3) && (cnt < m_nElementCount); j++)
			{
				for (i=0; i<m_nCountCommentsStrInTable; i++) DSReadLn(stream, str, 100);
				for (i=0; i<m_nNodesCount; i++) DSReadFieldStr(stream, &num,&node[0][i],&node[1][i],&node[2][i], &node[3][i]);
			
				field[0][cnt] = node[0][m_nElementNumberToRead-1];
				field[1][cnt] = node[1][m_nElementNumberToRead-1];
				field[2][cnt] = node[2][m_nElementNumberToRead-1];
				field[3][cnt] = node[3][m_nElementNumberToRead-1];

				cnt++;
			}
		}
		fclose(stream);
	}
	else return -1;	// Exit with error code (-1)
	return 0;	// Reading success
}

void DSfield(double U, double rx, double ry, double rz, double *Ex, double *Ey, double *Ez)
{// Ox и Oz меняем местами
	int j;
	double r1;
	double Ampl = U*normU;
	int nearN = 0;
	double nearR = 100.0;
	
	rx += dz;
	r1 = 100;
	for (j=0; j<m_nElementCount; j++)
	{
		r1 = (((xyz[0][j]-rz)*(xyz[0][j]-rz)) + ((xyz[1][j]-ry)*(xyz[1][j]-ry)) + ((xyz[2][j]-rx)*(xyz[2][j]-rx)));
		if (r1<nearR)
		{
			nearN=j;
			nearR=r1;
		}
	}
	*Ex = Ampl*field[2][nearN];
	*Ey = Ampl*field[1][nearN];
	*Ez = Ampl*field[0][nearN];
}
