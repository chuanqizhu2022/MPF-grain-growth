// Adapted from Prof. Koyama's MPF code in his textbook
// Author: Chuanqi Zhu
// Created on: 2022/5/19

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "CImg.h" //CImg ライブラリ（描画用）使用のためのヘッダ

using namespace cimg_library;

#define NDX 64
#define NDY 64
#define NDZ 64

#define N 2

int ndx = NDX;
int ndy = NDY;
int ndz = NDZ;
int ndmx = NDX - 1;
int ndmy = NDY - 1;
int ndmz = NDZ - 1;
int nm = N - 1;
double PI = 3.141592;
double RR = 8.3145;

double ****phi, ****phi2;
int ***phiNum, ****phiIdx;
double **aij, **wij, **mij, **fij;
double **anij, **thij, **vpij, **etaij;

CImg<unsigned char> phi_fldxy(NDX, NDY, 1, 3);
char outFilePhi_xy[64];

int phinum;
int i, j, k, ii, jj, kk;
int ni;
int ip, im, jp, jm, kp, km;
int n1, n2, n3;

int istep, nstep, pstep;
double dtime, L, dx;
double M0;				 //粒界の易動度
double W0;				 //ペナルティー項の係数
double A0;				 //勾配エネルギー係数
double F0;				 //粒界移動の駆動力
double temp;			 //温度
double sum1, sum2, sum3; //各種の和の作業変数
double pddtt;			 //フェーズフィールドの時間変化率

double gamma0; //粒界エネルギ密度
double delta;  //粒界幅（差分ブロック数にて表現）
double mobi;   //粒界の易動度
double vm0;	   //モル体積

double epsilon0;
double termiikk, termjjkk;

double phidx, phidy, phidz;
double phidxx, phidyy, phidzz;

//******* メインプログラム ******************************************
int main()
{
	phi = new double ***[N];
	phi2 = new double ***[N];
	for (ni = 0; ni <= nm; ni++)
	{
		phi[ni] = new double **[NDX];
		phi2[ni] = new double **[NDX];
		for (i = 0; i <= ndmx; i++)
		{
			phi[ni][i] = new double *[NDY];
			phi2[ni][i] = new double *[NDY];
			for (j = 0; j <= ndmy; j++)
			{
				phi[ni][i][j] = new double[NDZ];
				phi2[ni][i][j] = new double[NDZ];
			}
		}
	}

	phiIdx = new int ***[N + 1];
	for (ni = 0; ni <= N; ni++)
	{
		phiIdx[ni] = new int **[NDX];
		for (i = 0; i <= ndmx; i++)
		{
			phiIdx[ni][i] = new int *[NDY];
			for (j = 0; j <= ndmy; j++)
			{
				phiIdx[ni][i][j] = new int[NDZ];
			}
		}
	}

	phiNum = new int **[NDX];
	for (i = 0; i <= ndmx; i++)
	{
		phiNum[i] = new int *[NDY];

		for (j = 0; j <= ndmy; j++)
		{
			phiNum[i][j] = new int[NDZ];
		}
	}

	aij = new double *[N];
	wij = new double *[N];
	mij = new double *[N];
	fij = new double *[N];
	anij = new double *[N];
	thij = new double *[N];
	vpij = new double *[N];
	etaij = new double *[N];
	for (ni = 0; ni <= nm; ni++)
	{
		aij[ni] = new double[N];
		wij[ni] = new double[N];
		mij[ni] = new double[N];
		fij[ni] = new double[N];
		anij[ni] = new double[N];
		thij[ni] = new double[N];
		vpij[ni] = new double[N];
		etaij[ni] = new double[N];
	}

	nstep = 2001;
	pstep = 200;
	dtime = 5.0;
	temp = 1000.0;
	L = 2000.0;
	vm0 = 7.0e-6;
	delta = 7.0;
	mobi = 1.0;

	dx = L / 100.0 * 1.0e-9;			 //差分プロック１辺の長さ(m)
	gamma0 = 0.5 * vm0 / RR / temp / dx; //粒界エネルギ密度（0.5J/m^2）を無次元化x
	A0 = 8.0 * delta * gamma0 / PI / PI; //勾配エネルギー係数[式(4.40)]
	W0 = 4.0 * gamma0 / delta;			 //ペナルティー項の係数[式(4.40)]
	M0 = mobi * PI * PI / (8.0 * delta); //粒界の易動度[式(4.40)]
	F0 = 100.0 / RR / temp;				 //粒界移動の駆動力

	for (ii = 0; ii <= nm; ii++)
	{
		for (jj = 0; jj <= nm; jj++)
		{
			wij[ii][jj] = W0;
			aij[ii][jj] = A0;
			mij[ii][jj] = M0;
			fij[ii][jj] = 0.0;
			anij[ii][jj] = 0;
			thij[ii][jj] = 0.0;
			vpij[ii][jj] = 0.0;
			etaij[ii][jj] = 0.0;
			if ((ii == 0) || (jj == 0))
			{
				fij[ii][jj] = F0;
				anij[ii][jj] = 1;
			}
			if (ii < jj)
			{
				fij[ii][jj] = -fij[ii][jj];
			}
			if (ii == jj)
			{
				wij[ii][jj] = 0.0;
				aij[ii][jj] = 0.0;
				mij[ii][jj] = 0.0;
				fij[ii][jj] = 0.0;
				anij[ii][jj] = 0;
			}
		}
	}

	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			for (k = 0; k <= ndmz; k++)
			{
				if ((i - NDX / 2) * (i - NDX / 2) + (j - NDY / 2) * (j - NDY / 2) + (k - NDZ / 2) * (k - NDZ / 2) < 36)
				{
					phi[1][i][j][k] = 1.0;
					phi[0][i][j][k] = 0.0;
				}
				else
				{
					phi[1][i][j][k] = 0.0;
					phi[0][i][j][k] = 1.0;
				}
			}
		}
	}

start:;

	if (istep % pstep == 0)
	{
		cimg_forXY(phi_fldxy, x, y)
		{
			phi_fldxy(x, y, 0) = 255. * (phi[1][x][y][NDZ / 2]); // red
			phi_fldxy(x, y, 1) = 255. * (phi[1][x][y][NDZ / 2]); // green
			phi_fldxy(x, y, 2) = 255. * (phi[1][x][y][NDZ / 2]); // blue
		}
		sprintf(outFilePhi_xy, "figures/phi/2dxy%d.png", istep);
		phi_fldxy.save_jpeg(outFilePhi_xy);

		FILE *stream;
		char buffer[30];
		sprintf(buffer, "data/phi/3d%d.vtk", istep);
		stream = fopen(buffer, "a");

		fprintf(stream, "# vtk DataFile Version 1.0\n");
		fprintf(stream, "phi_%d.vtk\n", istep);
		fprintf(stream, "ASCII\n");
		fprintf(stream, "DATASET STRUCTURED_POINTS\n");
		fprintf(stream, "DIMENSIONS %d %d %d\n", NDX, NDY, NDZ);
		fprintf(stream, "ORIGIN 0.0 0.0 0.0\n");
		fprintf(stream, "ASPECT_RATIO 1.0 1.0 1.0\n");
		fprintf(stream, "\n");
		fprintf(stream, "POINT_DATA %d\n", NDX * NDY * NDZ);
		fprintf(stream, "SCALARS scalars float\n");
		fprintf(stream, "LOOKUP_TABLE default\n");

		for (i = 0; i <= ndmx; i++)
		{
			for (j = 0; j <= ndmy; j++)
			{
				for (k = 0; k <= ndmz; k++)
				{
					fprintf(stream, "%e\n", phi[1][i][j][k]);
				}
			}
		}
		fclose(stream);
	}

	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			for (k = 0; k <= ndmz; k++)
			{
				ip = i + 1;
				im = i - 1;
				jp = j + 1;
				jm = j - 1;
				kp = k + 1;
				km = k - 1;
				if (i == ndmx)
				{
					ip = 0;
				}
				if (i == 0)
				{
					im = ndmx;
				}
				if (j == ndmy)
				{
					jp = 0;
				}
				if (j == 0)
				{
					jm = ndmy;
				}
				if (k == ndmz)
				{
					kp = 0;
				}
				if (k == 0)
				{
					km = ndmz;
				}

				phinum = 0;
				for (ii = 0; ii <= nm; ii++)
				{
					if ((phi[ii][i][j][k] > 0.0) ||
						((phi[ii][i][j][k] == 0.0) && (phi[ii][ip][j][k] > 0.0) ||
						 (phi[ii][im][j][k] > 0.0) ||
						 (phi[ii][i][jp][k] > 0.0) ||
						 (phi[ii][i][jm][k] > 0.0) ||
						 (phi[ii][i][j][kp] > 0.0) ||
						 (phi[ii][i][j][km] > 0.0)))
					{
						phinum++;
						phiIdx[phinum][i][j][k] = ii;
					}
				}
				phiNum[i][j][k] = phinum;
			}
		}
	}

	// Evolution Equations
	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			for (k = 0; k <= ndmz; k++)
			{
				ip = i + 1;
				im = i - 1;
				jp = j + 1;
				jm = j - 1;
				kp = k + 1;
				km = k - 1;
				if (i == ndmx)
				{
					ip = 0;
				}
				if (i == 0)
				{
					im = ndmx;
				} //周期的境界条件
				if (j == ndmy)
				{
					jp = 0;
				}
				if (j == 0)
				{
					jm = ndmy;
				}
				if (k == ndmz)
				{
					kp = 0;
				}
				if (k == 0)
				{
					km = ndmz;
				}

				for (n1 = 1; n1 <= phiNum[i][j][k]; n1++)
				{
					ii = phiIdx[n1][i][j][k];
					pddtt = 0.0;
					for (n2 = 1; n2 <= phiNum[i][j][k]; n2++)
					{
						jj = phiIdx[n2][i][j][k];
						sum1 = 0.0;
						for (n3 = 1; n3 <= phiNum[i][j][k]; n3++)
						{
							kk = phiIdx[n3][i][j][k];

							// calculate the interface normal and deirivatives of the phase field
							phidx = (phi[kk][ip][j][k] - phi[kk][im][j][k]) / 2.0;
							phidy = (phi[kk][i][jp][k] - phi[kk][i][jm][k]) / 2.0;
							phidz = (phi[kk][i][j][kp] - phi[kk][i][j][km]) / 2.0;

							phidxx = (phi[kk][ip][j][k] + phi[kk][im][j][k] - 2.0 * phi[kk][i][j][k]); //フェーズフィールドの空間２階微分
							phidyy = (phi[kk][i][jp][k] + phi[kk][i][jm][k] - 2.0 * phi[kk][i][j][k]);
							phidzz = (phi[kk][i][j][kp] + phi[kk][i][j][km] - 2.0 * phi[kk][i][j][k]);

							termiikk = aij[ii][kk] * (phidxx + phidyy + phidzz);

							termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz);

							sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j][k];
						}
						pddtt += -2.0 * mij[ii][jj] / (double)phiNum[i][j][k] * (sum1 - 8.0 / PI * fij[ii][jj] * sqrt(phi[ii][i][j][k] * phi[jj][i][j][k]));
						//フェーズフィールドの発展方程式[式(4.31)]
					}
					phi2[ii][i][j][k] = phi[ii][i][j][k] + pddtt * dtime; //フェーズフィールドの時間発展（陽解法）
					if (phi2[ii][i][j][k] >= 1.0)
					{
						phi2[ii][i][j][k] = 1.0;
					} //フェーズフィールドの変域補正
					if (phi2[ii][i][j][k] <= 0.0)
					{
						phi2[ii][i][j][k] = 0.0;
					}
				}
			} // j
		}	  // i
	}

	for (ni = 0; ni <= nm; ni++)
	{
		for (i = 0; i <= ndmx; i++)
		{
			for (j = 0; j <= ndmy; j++)
			{
				for (k = 0; k <= ndmz; k++)
				{
					phi[ni][i][j][k] = phi2[ni][i][j][k];
				}
			}
		}
	}

	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			for (k = 0; k <= ndmz; k++)
			{
				sum1 = 0.0;
				for (ni = 0; ni <= nm; ni++)
				{
					sum1 += phi[ni][i][j][k];
				}
				for (ni = 0; ni <= nm; ni++)
				{
					phi[ni][i][j][k] = phi[ni][i][j][k] / sum1;
				}
			}
		}
	}

	istep = istep + 1;
	if (istep < nstep)
	{
		goto start;
	}

end:;
	return 0;
}
