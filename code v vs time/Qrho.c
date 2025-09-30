/* Binary Fluid - coupled to hydrodynamics and P - */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* useful constants */
#define Pi 3.141592653589793
#define Piover2 1.570796326794897
#define threePiover2 4.71238898
#define twoPi 6.283185307
#define fivePiover2 7.853981634

/*
 *	Lattice Boltzmann
 */

/* program parameters */
#define Lx 64
#define Ly 64 // 100 //96
#define Lz 1 // 100 //96				/*32*/ /*240*/ /*100*/

#define Nmax 1000000			/*2000000*/ /* total number of iterations */
#define stepskip 500			/*2500*/ /* output graphics every stepskip steps */
#define tau1 2.5			/* 2.5 proportional to viscosity */
#define dt 1.0
#define densityinit 2.0			/* initial value for the density everywhere */

int iupa[Lx],idwna[Lx];
int jupa[Ly],jdwna[Ly];
int kupa[Lz],kdwna[Lz];

/* physical parameters */
#define temperature 0.3333333 //0.5
double bodyforce = 0.0000;			/* body force to set Pouiseuille flow */

double friction;
double amplitude;			/* amplitude for perturbation */

double vwtp = 0.0;
double vwbt = 0.0;

/* function declarations */
void observables(double [Lx][Ly][Lz][15]);
void equilibriumdistQ(void);

void initialize(void);

void collision  (double [Lx][Ly][Lz][15], double [Lx][Ly][Lz][15]);
void collisionpr(double [Lx][Ly][Lz][15], double [Lx][Ly][Lz][15]);

void update0(double [Lx][Ly][Lz][15], double [Lx][Ly][Lz][15]);
void update (double [Lx][Ly][Lz][15], double [Lx][Ly][Lz][15], double [Lx][Ly][Lz][15]);

void computedensityfluctuations(int boxsize);

void plot(void);
void plotQ(void);
void multipleplot(int n);
void plotdensityfluctuations(void);

void plotbifurcation(void);
void plotV(int n);

void correction(int);
double globalP[3];
double densitycorrection;

#define n 3
void jacobi(double (*a)[n], double d[], double (*v)[n], int *nrot);
#undef n

/* lattice Boltzmann variables */
double C[Lx][Ly][Lz][15];
double Cpr[Lx][Ly][Lz][15];
double feq[Lx][Ly][Lz][15];

double f[Lx][Ly][Lz][15];
double fold[Lx][Ly][Lz][15];
double fpr[Lx][Ly][Lz][15];
double fnew[Lx][Ly][Lz][15];

double flowparameter[Lx][Ly][Lz];
double compressibility[Lx][Ly][Lz];

int e[15][3];

/* observables */
double density[Lx][Ly][Lz];
double ux[Lx][Ly][Lz];
double uy[Lx][Ly][Lz];
double uz[Lx][Ly][Lz];
double umax;

void plotumax(int n);

/* switch for boundary conditions*/
#define BC 0
#define FIXEDBC 0

double wallamp;

double angztop=90.0; // anchoring at top wall: 0.0 is homeotropic/normal
double angzbot=90.0; // anchoring at bottom wall

double angxytop=90.0;
double angxybot=90.0;

/*
 *	Binary Fluid
 */

/* binary fluid parameters */
#define alpha2phi 0.01 //0.07
#define alpha4phi 0.0 //0.07
#define Kphi 0.00 // 0.01 // 0.14

double M = 1.0; //0.03;				/*0.05*/ /* need to tweak */
double radius = 40;

double noisestrength=0.0;

double phi0=1.0; /*initial value of the concentration*/
 
/* switch for initial conditions*/
#define DROPLET 0
#define SPINODAL 0
#define TWIST 0

/* variables */
double phi[Lx][Ly][Lz];
double phiold[Lx][Ly][Lz];

int ndatatot;
double phiav[Ly],phivar[Ly];

double mu[Lx][Ly][Lz];		/* chemical potential */

double sxx[Lx][Ly][Lz];		/* dxphidxphi */
double sxy[Lx][Ly][Lz];
double sxz[Lx][Ly][Lz];
double syy[Lx][Ly][Lz];
double syz[Lx][Ly][Lz];
double szz[Lx][Ly][Lz];

double Ja[Lx][Ly][Lz][3];	/* advective current */
double Jd[Lx][Ly][Lz][3];	/* diffusive current */

double h_phi[Lx][Ly][Lz];
double h_phiold[Lx][Ly][Lz];

double Vx = 0.0;
double Vy = 0.0;
double Vz = 0.0;

int R1;

double phi_average;
double phi_average0;

/* functions */
void phi_function(void);

void initializephi(void);

void updatephi0(double phipr[Lx][Ly][Lz], double phiold[Lx][Ly][Lz]);
void updatephi(double phinew[Lx][Ly][Lz], double phiold[Lx][Ly][Lz]);

/* mathematical functions */
double max(double, double);
double min(double, double);
double dx_(double phi[Lx][Ly][Lz], int i, int j, int k);
double dy_(double phi[Lx][Ly][Lz], int i, int j, int k);
double dz_(double phi[Lx][Ly][Lz], int i, int j, int k);
double laplacian_(double phi[Lx][Ly][Lz], int i , int j , int k);

/*
 *	Leslie-Ericksen
 */

/* Leslie-Ericksen parameters */
#define alpha2 0.1
#define alpha4 0.1
#define K 0.04

double w1 = 0.0;			/* self advection */
double zeta = 0.01;			/*0.0*/ /* activity */
double xi = 0.7; 		        /* determines whether flow-aligning or flow-tumbling, keep fixed to begin with; xi>0: rod-like, xi<0: disc-like */
double G = 1.0;				/* rotational diffusion constant */

/* function declarations */
void P_function(void);
void Q_function(void);

void initializeQ(void);

void updateP0(double Pxpr[Lx][Ly][Lz],  double Pypr[Lx][Ly][Lz],  double Pzpr[Lx][Ly][Lz],
	      double Pxold[Lx][Ly][Lz], double Pyold[Lx][Ly][Lz], double Pzold[Lx][Ly][Lz]);
void updateP (double Pxnew[Lx][Ly][Lz], double Pynew[Lx][Ly][Lz], double Pznew[Lx][Ly][Lz],
	      double Pxold[Lx][Ly][Lz], double Pyold[Lx][Ly][Lz], double Pzold[Lx][Ly][Lz]);

void updateQ0(double Qxxpr[Lx][Ly][Lz],  double Qxypr[Lx][Ly][Lz],  double Qxzpr[Lx][Ly][Lz], double Qyypr[Lx][Ly][Lz],  double Qyzpr[Lx][Ly][Lz],
	      double Qxxold[Lx][Ly][Lz], double Qxyold[Lx][Ly][Lz], double Qxzold[Lx][Ly][Lz], double Qyyold[Lx][Ly][Lz], double Qyzold[Lx][Ly][Lz]);
void updateQ(double Qxxnew[Lx][Ly][Lz], double Qxynew[Lx][Ly][Lz], double Qxznew[Lx][Ly][Lz], double Qyynew[Lx][Ly][Lz], double Qyznew[Lx][Ly][Lz],
	     double Qxxold[Lx][Ly][Lz], double Qxyold[Lx][Ly][Lz], double Qxzold[Lx][Ly][Lz], double Qyyold[Lx][Ly][Lz], double Qyzold[Lx][Ly][Lz]);

void add_perturbation(void);
void add_Qperturbation(void);

/* Leslie-Ericksen variables */
double Px[Lx][Ly][Lz], Py[Lx][Ly][Lz], Pz[Lx][Ly][Lz];
double Pxold[Lx][Ly][Lz], Pyold[Lx][Ly][Lz], Pzold[Lx][Ly][Lz];

double hmol[Lx][Ly][Lz][3];
double h[Lx][Ly][Lz][3];
double Fh[Lx][Ly][Lz][15];

double hold[Lx][Ly][Lz][3];
double hnew[Lx][Ly][Lz][3];

double tauxy[Lx][Ly][Lz], tauxz[Lx][Ly][Lz];
double tauyz[Lx][Ly][Lz];

/*Q tensor variables*/
double Qxx[Lx][Ly][Lz],Qxy[Lx][Ly][Lz],Qxz[Lx][Ly][Lz],Qyy[Lx][Ly][Lz],Qyz[Lx][Ly][Lz];
double Qxxold[Lx][Ly][Lz],Qxyold[Lx][Ly][Lz],Qxzold[Lx][Ly][Lz],Qyyold[Lx][Ly][Lz],Qyzold[Lx][Ly][Lz];
double Qgradphix[Lx][Ly][Lz],Qgradphiy[Lx][Ly][Lz],Qgradphiz[Lx][Ly][Lz];
double DEHxx[Lx][Ly][Lz],DEHxy[Lx][Ly][Lz],DEHxz[Lx][Ly][Lz],DEHyy[Lx][Ly][Lz],DEHyz[Lx][Ly][Lz];
double DEH2xx[Lx][Ly][Lz],DEH2xy[Lx][Ly][Lz],DEH2xz[Lx][Ly][Lz],DEH2yy[Lx][Ly][Lz],DEH2yz[Lx][Ly][Lz];
double DEHxxold[Lx][Ly][Lz],DEHxyold[Lx][Ly][Lz],DEHxzold[Lx][Ly][Lz],DEHyyold[Lx][Ly][Lz],DEHyzold[Lx][Ly][Lz];
double txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz;
double laplacianphi[Lx][Ly][Lz];

double DG2xx[Lx][Ly][Lz],DG2yy[Lx][Ly][Lz],DG2xy[Lx][Ly][Lz];
double DG2xz[Lx][Ly][Lz],DG2yz[Lx][Ly][Lz],DG2zz[Lx][Ly][Lz],pG[Lx][Ly][Lz];
double Abulk = 0.000; //1.0; /*determines bulk liquid crystal free energy*/
double L1 = 0.01; //0.001 /*liquid crystal elastic constant*/
double W = -0.00; // -0.01 /*anchoring term at interfaces*/
double phivr = 1.0;
double gammac;

#define gamma 2.0 // 1.0 /*baseline value of gamma - N/I transition at 2.7*/
#define delta 1.0 // 0.5

double energy;

/* Switch */
#define PERTURBED 0

/*
 *	coupling parameters
 */

double beta = 0.0;		/* 0.05, not 0.1 */

FILE *output1;
char filename1[20];
FILE *output2;
char filename2[20];
FILE *output3;
char filename3[20];

int main(int argc, char** argv)
{
	int i, j, k, l;
	int n;

	double Qsqxx,Qsqxy,Qsqxz,Qsqyy,Qsqyz,Qsqzz;
	double Qxxl,Qxyl,Qxzl,Qyyl,Qyzl,Qzzl;
	double TrQ2,qorder;

	double deltaQ;


	w1 = 0.0;
	friction = 0.0; //0.01;

	double tol = 1.0E-12;
	double tolerance, Vybefore, Vzbefore;

    scanf("%lf %lf",&zeta,&bodyforce);
	initialize();
	initializephi();
	initializeQ();

	phi_function();
	Q_function();
	observables(f);

	phi_average0 = phi_average;

	/*********** plot initial configuration ************/
	plotQ();
	/***************************************************/

	angztop=0.0;
	angzbot=0.0;

	ndatatot =0;
	for(j = 0; j<Ly; j++){
	  phiav[j]=0.0;
	  phivar[j]=0.0;
	}

	//main update loop starts on next row
	for (n = 1; n <= Nmax; n++) {

		/************** add noise ***************/

		/*if (n == 5000)  {add_Qperturbation();}
		if (n == 20000)  {add_Qperturbation();}
		if (n == 50000)  {add_Qperturbation();}*/

		//if (n == 100000) {add_perturbation();}
		//if (n == 200000) {add_perturbation();}
		//if (n == 300000) {add_perturbation();friction = 0.0;}
		/************** add noise ***************/

		/************** Cahn-Hilliard equation ****************************************/
		for (i = 0; i < Lx; i++) {
			for (j = 0; j < Ly; j++) {
				for (k = 0; k < Lz; k++) {
					phiold[i][j][k] = phi[i][j][k];
				}
			}
		}

		/* predictor */
		phi_function();					/* find mu */
		updatephi0(phi, phiold);			/* find phi primed */

		/* corrector */
		//phi_function();				/* find mu primed */
		//updatephi(phi, phiold);			/* find phi */
		/************** Cahn-Hiliard equation ****************************************/

		phi_function();
		Q_function();


		/************************ update Q tensor ******************************/
		for (i = 0; i < Lx; i++) {
			for (j = 0; j < Ly; j++) {
				for (k = 0; k < Lz; k++) {
					Qxxold[i][j][k] = Qxx[i][j][k];
					Qxyold[i][j][k] = Qxy[i][j][k];
					Qxzold[i][j][k] = Qxz[i][j][k];
					Qyyold[i][j][k] = Qyy[i][j][k];
					Qyzold[i][j][k] = Qyz[i][j][k];
				}
			}
		}
		/* predictor */
		Q_function();							/* find molecular field */
		updateQ0(Qxx, Qxy, Qxz, Qyy, Qyz, Qxxold, Qxyold, Qxzold, Qyyold, Qyzold);			/* find Q primed */

		/* below to ensure 2D */
		for (i = 0; i < Lx; i++) {
		  for (j = 0; j < Ly; j++) {
		    for (k = 0; k < Lz; k++) {
		      Qxxl=Qxx[i][j][k];
		      Qxyl=Qxy[i][j][k];
		      Qyyl=Qyy[i][j][k];
		      Qxzl=Qxz[i][j][k];
		      Qyzl=Qyz[i][j][k];
	      
		      Qzzl=-Qxx[i][j][k]-Qyy[i][j][k];
		      
		      Qsqxx=Qxxl+Qyyl;
		      Qsqzz=Qsqxx*Qsqxx+Qxzl*Qxzl+Qyzl*Qyzl;
		      Qsqxy=Qxyl*Qsqxx+Qxzl*Qyzl;
		      Qsqxz=Qxxl*Qxzl-Qxzl*Qsqxx+Qxyl*Qyzl;
		      Qsqyz=Qxyl*Qxzl-Qyzl*Qsqxx+Qyyl*Qyzl;
		      Qsqxx=Qxxl*Qxxl+Qxyl*Qxyl+Qxzl*Qxzl;
		      Qsqyy=Qyyl*Qyyl+Qxyl*Qxyl+Qyzl*Qyzl;
		      TrQ2=Qsqxx+Qsqyy+Qsqzz;
		      
		      qorder=sqrt(TrQ2*3.0/2.0);
				
		      Qxz[i][j][k]=0.0;
		      Qyz[i][j][k]=0.0;
		      deltaQ=-Qxx[i][j][k]-Qyy[i][j][k]+qorder*1.0/3.0;
		      Qxx[i][j][k]+=deltaQ/2.0;
		      Qyy[i][j][k]+=deltaQ/2.0;
		
		      /*Qxz[i][j][k]=0.0;
		      Qxy[i][j][k]=0.0;
		      deltaQ=Qxx[i][j][k]+qorder*1.0/3.0;
		      Qxx[i][j][k]-=deltaQ;
		      Qyy[i][j][k]+=deltaQ/2.0;*/

		      /*Qxy[i][j][k]=0.0;
		      Qyz[i][j][k]=0.0;
		      deltaQ=Qyy[i][j][k]+qorder*1.0/3.0;
		      Qxx[i][j][k]+=deltaQ/2.0;
		      Qyy[i][j][k]-=deltaQ;*/


		    }
		  }
		}
 
		
		
		/* corrector */
		//P_function();							/* find h primed */
		//updateP(Px, Py, Pz, Pxold, Pyold, Pzold);			/* find P */

		/************************ update director field ******************************/

		phi_function();
		Q_function();

		/********************************* LB *******************************************/
		for (i = 0; i < Lx; i++) {
			for (j = 0; j < Ly; j++) {
				for (k = 0; k < Lz; k++) {
					for (l = 0; l < 15; l++) {
						fold[i][j][k][l] = f[i][j][k][l];
					}
				}
			}
		}

		/* predictor */
		observables(fold);
		equilibriumdistQ();				/* find feq */
		collision(feq,fold);				/* find C */

		update0(fpr,fold);				/* find f primed */

		/* corrector */
		observables(fpr);
		equilibriumdistQ();
		collisionpr(feq,fpr);				/* find C primed */

		update(f,fpr,fold);				/* find f */
		/********************************* LB *******************************************/

		observables(f);
		correction(n);

		/***************** plot ***************************************/
		if (n % stepskip == 0) {
		  plotQ();
		  if(n>=10000) plotdensityfluctuations();
		  if(n>=10000) computeaveragevelocity();
            
			tolerance = (Vy - Vybefore)/Vybefore;
			Vybefore = Vy;

           		printf("%d \t %E \t %E \t %E \n", n, phi_average, Vy, tolerance);

			/*			if (tolerance > -tol && tolerance < tol) {
				break;
				}*/
		}
		/**************************************************************/
	}

	printf("final Vy = %E \n \n", Vy);
}

/*
 * 	Lattice Boltzmann
 */

/* initialization */
void initialize(void)
{
	int i, j, k, l;

	e[0][0]= 0;
	e[0][1]= 0;
	e[0][2]= 0;

	e[1][0]= 1;
	e[1][1]= 0;
	e[1][2]= 0;

	e[2][0]= 0;
	e[2][1]= 1;
	e[2][2]= 0;

	e[3][0]= -1;
	e[3][1]= 0;
	e[3][2]= 0;

	e[4][0]= 0;
	e[4][1]= -1;
	e[4][2]= 0;

	e[5][0]= 0;
	e[5][1]= 0;
	e[5][2]= 1;

	e[6][0]= 0;
	e[6][1]= 0;
	e[6][2]= -1;

	e[7][0]= 1;
	e[7][1]= 1;
	e[7][2]= 1;

	e[8][0]= -1;
	e[8][1]= 1;
	e[8][2]= 1;

	e[9][0]= -1;
	e[9][1]= -1;
	e[9][2]= 1;

	e[10][0]= 1;
	e[10][1]= -1;
	e[10][2]= 1;

	e[11][0]= 1;
	e[11][1]= 1;
	e[11][2]= -1;

	e[12][0]= -1;
	e[12][1]= 1;
	e[12][2]= -1;

	e[13][0]= -1;
	e[13][1]= -1;
	e[13][2]= -1;

	e[14][0]= 1;
	e[14][1]= -1;
	e[14][2]= -1;

	for (i = 0; i < Lx; i++) {
		if (i == Lx-1) {iupa[i]  = 0;}		else {iupa[i]  = i+1;}
		if (i == 0)    {idwna[i] = Lx-1;}	else {idwna[i] = i-1;}
		for (j = 0; j < Ly; j++) {
			if (j == Ly-1) {jupa[j]  = 0;}		else {jupa[j]  = j+1;}
			if (j == 0)    {jdwna[j] = Ly-1;}	else {jdwna[j] = j-1;}
			for (k = 0; k < Lz; k++) {
				if (k == Lz-1) {kupa[k]  = 0;}		else {kupa[k]  = k+1;}
				if (k == 0)    {kdwna[k] = Lz-1;}	else {kdwna[k] = k-1;}

				density[i][j][k] = densityinit;

				/* initialise f */
				for (l = 0; l < 15; l++) {
					f[i][j][k][l] = density[i][j][k]/15.0;
				}
			}
		}
	}
}

/* predictor-corrector of the LB step */
void update0(double fpr[Lx][Ly][Lz][15], double fold[Lx][Ly][Lz][15])
{
	int i, j, k, l, imod, jmod, kmod;
	double rb;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				/* main LB equation */
				for (l = 0; l < 15; l++) {
					imod = (i - e[l][0] + Lx) % Lx;
					jmod = (j - e[l][1] + Ly) % Ly;
					kmod = (k - e[l][2] + Lz) % Lz;

					fpr[i][j][k][l] = fold[imod][jmod][kmod][l] + dt*C[imod][jmod][kmod][l];
				}
	/* moving wall BC */
#if BC
	if (k == 0) {
	  fpr[i][j][k][5] = fpr[i][j][k][6];

	  rb = fpr[i][j][k][0] + fpr[i][j][k][1] + fpr[i][j][k][2] + fpr[i][j][k][3] + fpr[i][j][k][4] +
	    2.0 * (fpr[i][j][k][6] + fpr[i][j][k][11] + fpr[i][j][k][12] + fpr[i][j][0][13] + fpr[i][j][0][14]);

	  fpr[i][j][k][7] = 0.25 * (-fpr[i][j][k][1] - fpr[i][j][k][2] + fpr[i][j][k][3] + fpr[i][j][k][4]) +
	    0.25 * (-fpr[i][j][k][11] + fpr[i][j][k][12] + fpr[i][j][k][14] + rb * vwbt) + 0.75 * fpr[i][j][k][13];
	  
	  fpr[i][j][k][8] = 0.25 * (fpr[i][j][k][1] - fpr[i][j][k][2] - fpr[i][j][k][3] + fpr[i][j][k][4]) +
	    0.25 * (fpr[i][j][k][11] - fpr[i][j][k][12] + fpr[i][j][k][13] + rb * vwbt) + 0.75 * fpr[i][j][k][14];

	  fpr[i][j][k][9] = 0.25 * (fpr[i][j][k][1] + fpr[i][j][k][2] - fpr[i][j][k][3] - fpr[i][j][k][4]) +
	    0.25 * (-fpr[i][j][k][13] + fpr[i][j][k][14] + fpr[i][j][k][12] - rb * vwbt) + 0.75 * fpr[i][j][k][11];
	  fpr[i][j][k][10] = 0.25 * (-fpr[i][j][k][1] + fpr[i][j][k][2] + fpr[i][j][k][3] - fpr[i][j][0][4]) +
	    0.25 * (fpr[i][j][k][11] + fpr[i][j][k][13] - fpr[i][j][k][14] - rb * vwbt) +
	    0.75 * fpr[i][j][k][12];
	}
	else if (k == Lz - 1) {
	  fpr[i][j][k][6] = fpr[i][j][k][5];

                        rb = fpr[i][j][k][0] + fpr[i][j][k][1] + fpr[i][j][k][2] + fpr[i][j][k][3] + fpr[i][j][k][4] +
			  2.0 * (fpr[i][j][k][5] + fpr[i][j][k][7] + fpr[i][j][k][8] + fpr[i][j][k][9] + fpr[i][j][k][10]);

                        fpr[i][j][k][11] = 0.25 * (-fpr[i][j][k][1] - fpr[i][j][k][2] + fpr[i][j][k][3] + fpr[i][j][k][4]) +
			  0.25 * (-fpr[i][j][k][7] + fpr[i][j][k][8] + fpr[i][j][k][10] + rb * vwtp) +
			  0.75 * fpr[i][j][k][9];

                        fpr[i][j][k][12] = 0.25 * (fpr[i][j][k][1] - fpr[i][j][k][2] - fpr[i][j][k][3] + fpr[i][j][k][4]) +
			  0.25 * (fpr[i][j][k][7] - fpr[i][j][k][8] + fpr[i][j][k][9] + rb * vwtp) +
			  0.75 * fpr[i][j][k][10];

                        fpr[i][j][k][13] = 0.25 * (fpr[i][j][k][1] + fpr[i][j][k][2] - fpr[i][j][k][3] - fpr[i][j][k][4]) +
			  0.25 * (-fpr[i][j][k][9] + fpr[i][j][k][10] + fpr[i][j][k][8] - rb * vwtp) +
			  0.75 * fpr[i][j][k][7];

                        fpr[i][j][k][14] = 0.25 * (-fpr[i][j][k][1] + fpr[i][j][k][2] + fpr[i][j][k][3] - fpr[i][j][k][4]) +
			  0.25 * (fpr[i][j][k][7] + fpr[i][j][k][9] - fpr[i][j][k][10] - rb * vwtp) +
			  0.75 * fpr[i][j][k][8];
	}             
#endif
			}
		}
	}

}
void update(double fnew[Lx][Ly][Lz][15], double fpr[Lx][Ly][Lz][15], double fold[Lx][Ly][Lz][15])
{
	int i, j, k, l, imod, jmod, kmod;
	double rb;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				/* main LB equation */
				for (l = 0; l < 15; l++) {
					imod = (i - e[l][0] + Lx) % Lx;
					jmod = (j - e[l][1] + Ly) % Ly;
					kmod = (k - e[l][2] + Lz) % Lz;

					fnew[i][j][k][l] = fold[imod][jmod][kmod][l] + 0.5*dt*(C[imod][jmod][kmod][l] + Cpr[i][j][k][l]);
				}
				/* moving wall BC */
#if BC
				if (k == 0) {
				  fnew[i][j][k][5] = fnew[i][j][0][6];

                    rb = fnew[i][j][k][0] + fnew[i][j][k][1] + fnew[i][j][k][2] + fnew[i][j][k][3] + fnew[i][j][k][4] +
		      2.0 * (fnew[i][j][k][6] + fnew[i][j][k][11] + fnew[i][j][k][12] + fnew[i][j][k][13] + fnew[i][j][k][14]);

                    fnew[i][j][k][7] = 0.25 * (-fnew[i][j][k][1] - fnew[i][j][k][2] + fnew[i][j][k][3] + fnew[i][j][k][4]) +
		      0.25 * (-fnew[i][j][k][11] + fnew[i][j][k][12] + fnew[i][j][k][14] + rb * vwbt) +
		      0.75 * fnew[i][j][k][13];

                    fnew[i][j][k][8] = 0.25 * (fnew[i][j][k][1] - fnew[i][j][k][2] - fnew[i][j][k][3] + fnew[i][j][k][4]) +
		      0.25 * (fnew[i][j][k][11] - fnew[i][j][k][12] + fnew[i][j][k][13] + rb * vwbt) +
		      0.75 * fnew[i][j][k][14];

                    fnew[i][j][k][9] = 0.25 * (fnew[i][j][k][1] + fnew[i][j][k][2] - fnew[i][j][k][3] - fnew[i][j][k][4]) +
		      0.25 * (-fnew[i][j][k][13] + fnew[i][j][k][14] + fnew[i][j][k][12] - rb * vwbt) +
		      0.75 * fnew[i][j][k][11];

                    fnew[i][j][k][10] = 0.25 * (-fnew[i][j][k][1] + fnew[i][j][k][2] + fnew[i][j][k][3] - fnew[i][j][k][4]) +
		      0.25 * (fnew[i][j][k][11] + fnew[i][j][k][13] - fnew[i][j][k][14] - rb * vwbt) +
		      0.75 * fnew[i][j][k][12];
				} else if (k == Lz - 1) {
				  fnew[i][j][k][6] = fnew[i][j][k][5];

                    rb = fnew[i][j][k][0] + fnew[i][j][k][1] + fnew[i][j][k][2] + fnew[i][j][k][3] + fnew[i][j][k][4] +
		      2.0 * (fnew[i][j][k][5] + fnew[i][j][k][7] + fnew[i][j][k][8] + fnew[i][j][k][9] + fnew[i][j][k][10]);

                    fnew[i][j][k][11] = 0.25 * (-fnew[i][j][k][1] - fnew[i][j][k][2] + fnew[i][j][k][3] + fnew[i][j][k][4]) +
		      0.25 * (-fnew[i][j][k][7] + fnew[i][j][k][8] + fnew[i][j][k][10] + rb * vwtp) +
		      0.75 * fnew[i][j][k][9];

                    fnew[i][j][k][12] = 0.25 * (fnew[i][j][k][1] - fnew[i][j][k][2] - fnew[i][j][k][3] + fnew[i][j][k][4]) +
		      0.25 * (fnew[i][j][k][7] - fnew[i][j][k][8] + fnew[i][j][k][9] + rb * vwtp) +
		      0.75 * fnew[i][j][k][10];

                    fnew[i][j][k][13] = 0.25 * (fnew[i][j][k][1] + fnew[i][j][k][2] - fnew[i][j][k][3] - fnew[i][j][k][4]) +
		      0.25 * (-fnew[i][j][k][9] + fnew[i][j][k][10] + fnew[i][j][k][8] - rb * vwtp) +
		      0.75 * fnew[i][j][k][7];

                    fnew[i][j][k][14] = 0.25 * (-fnew[i][j][k][1] + fnew[i][j][k][2] + fnew[i][j][k][3] - fnew[i][j][k][4]) +
		      0.25 * (fnew[i][j][k][7] + fnew[i][j][k][9] - fnew[i][j][k][10] - rb * vwtp) +
		      0.75 * fnew[i][j][k][8];
				}
#endif
			}
		}
	}
}

/* collision operator */
void collision(double feq[Lx][Ly][Lz][15], double f[Lx][Ly][Lz][15])
{
	int i, j, k, l;

	for (i = 0; i < Lx; i++)
		for (j = 0; j < Ly; j++)
			for (k = 0; k < Lz; k++)
				for (l = 0; l < 15; l++) {
					C[i][j][k][l] = (feq[i][j][k][l] - f[i][j][k][l])/tau1;
				}
}
void collisionpr(double feqpr[Lx][Ly][Lz][15], double fpr[Lx][Ly][Lz][15])
{
	int i, j, k, l;

	for (i = 0; i < Lx; i++)
		for (j = 0; j < Ly; j++)
			for (k = 0; k < Lz; k++)
				for (l = 0; l < 15; l++) {
					Cpr[i][j][k][l] = (feqpr[i][j][k][l] - fpr[i][j][k][l])/tau1;
				}
}

/* calculate equilibrium distribution for Q tensor*/
void equilibriumdistQ(void)
{
	double A0,A1,A2,B1,B2,C0,C1,C2,D1,D2;
	double G1xx,G2xx,G2xy,G2xz,G2yz,G1yy,G2yy,G1zz,G2zz;

	double rho,phil,Qxxl,Qxyl,Qyyl,Qxzl,Qyzl,Qzzl,usq,udote,omdote;
	double nnxxl,nnyyl;
	double Hxx,Hyy,Hxy,Hxz,Hyz,Qsqxx,Qsqxy,Qsqyy,Qsqzz,Qsqxz,Qsqyz,TrQ2;

	double sigxx,sigyy,sigxy,sigxz,sigyz,sigzz;

   	double Force[3];

	double dbdtauxb,dbdtauyb,dbdtauzb;
	double Pxl,Pyl,Pzl,Psquare,Pdoth;
	double hx,hy,hz;


	int i, j, k, l;
	int iup, jup, kup;
   	int iup2, jup2, kup2;
    	int idwn, jdwn, kdwn;
    	int idwn2, jdwn2, kdwn2;

	for (i = 0; i < Lx; i++) {
	  for (j = 0; j < Ly; j++) {
	    for (k = 0; k < Lz; k++) {
	      int iup = iupa[i];
	      int jup = jupa[j];
	      int kup = kupa[k];
	      
	      int iup2 = iupa[iupa[i]];
	      int jup2 = jupa[jupa[j]];
	      int kup2 = kupa[kupa[k]];
	      
	      int idwn = idwna[i];
	      int jdwn = jdwna[j];
	      int kdwn = kdwna[k];
	      
	      int idwn2 = idwna[idwna[i]];
	      int jdwn2 = jdwna[jdwna[j]];
	      int kdwn2 = kdwna[kdwna[k]];
	      
	      rho = density[i][j][k];
	      phil = phi[i][j][k];

	      Qxxl=Qxx[i][j][k];
	      Qxyl=Qxy[i][j][k];
	      Qyyl=Qyy[i][j][k];
	      Qxzl=Qxz[i][j][k];
	      Qyzl=Qyz[i][j][k];
	      
	      Qzzl=-Qxx[i][j][k]-Qyy[i][j][k];

	      nnxxl= Qxxl+1.0/3.0;
	      nnyyl= Qyyl+1.0/3.0;

	      Qsqxx=Qxxl+Qyyl;
	      Qsqzz=Qsqxx*Qsqxx+Qxzl*Qxzl+Qyzl*Qyzl;
	      Qsqxy=Qxyl*Qsqxx+Qxzl*Qyzl;
	      Qsqxz=Qxxl*Qxzl-Qxzl*Qsqxx+Qxyl*Qyzl;
	      Qsqyz=Qxyl*Qxzl-Qyzl*Qsqxx+Qyyl*Qyzl;
	      Qsqxx=Qxxl*Qxxl+Qxyl*Qxyl+Qxzl*Qxzl;
	      Qsqyy=Qyyl*Qyyl+Qxyl*Qxyl+Qyzl*Qyzl;
	      TrQ2=Qsqxx+Qsqyy+Qsqzz;
	      
	      Hxx=DEH2xx[i][j][k];
	      Hxy=DEH2xy[i][j][k];
	      Hxz=DEH2xz[i][j][k];
	      Hyy=DEH2yy[i][j][k];
	      Hyz=DEH2yz[i][j][k];
	      
	      sigxx=2.0/3.0*xi*((1.0+3.0*Qxxl)*(Hxx*(1.0-2.0*Qxxl-Qyyl)-Hyy*(Qxxl+2.0*Qyyl)-2.0*Hyz*Qyzl)+(Hxy*Qxyl+Hxz*Qxzl)*(1.0-6.0*Qxxl));	/* there should be a minus sign in actual equation */
	      sigxy=xi*(Hxy*(2.0/3.0+Qxxl-4.0*Qxyl*Qxyl+Qyyl)-Hxx*Qxyl*(-1.0+4.0*Qxxl+2.0*Qyyl)-Hyy*Qxyl*(-1.0+4.0*Qyyl+2.0*Qxxl)+Hxz*(-4.0*Qxyl*Qxzl+Qyzl)+Hyz*(Qxzl-4.0*Qxyl*Qyzl));
	      sigyy=2.0/3.0*xi*((1.0+3.0*Qyyl)*(Hyy*(1.0-Qxxl-2.0*Qyyl)-Hxx*(2.0*Qxxl+Qyyl)-2.0*Hxz*Qxzl)+(Hxy*Qxyl+Hyz*Qyzl)*(1.0-6.0*Qyyl));
	      sigxz=xi*(Hxz*(2.0/3.0-4.0*Qxzl*Qxzl-Qyyl)-Hxx*Qxzl*(4.0*Qxxl+2.0*Qyyl)-Hyy*Qxzl*(1.0+4.0*Qyyl+2.0*Qxxl)+Hxy*(Qyzl-4.0*Qxyl*Qxzl)+Hyz*(Qxyl-4.0*Qxzl*Qyzl));
	      sigyz=xi*(Hyz*(2.0/3.0-4.0*Qyzl*Qyzl-Qxxl)-Hyy*Qyzl*(4.0*Qyyl+2.0*Qxxl)-Hxx*Qyzl*(1.0+4.0*Qxxl+2.0*Qyyl)+Hxy*(Qxzl-4.0*Qxyl*Qyzl)+Hxz*(Qxyl-4.0*Qxzl*Qyzl));
	      sigzz=-(sigxx+sigyy); 
	
	      if(phi[i][j][k]>0.0){
	      sigxx += zeta*phi[i][j][k]* Qxxl;
	      sigxy += zeta*phi[i][j][k]* Qxyl;
	      sigxz += zeta*phi[i][j][k]* Qxzl;
	      sigyy += zeta*phi[i][j][k]* Qyyl;
	      sigyz += zeta*phi[i][j][k]* Qyzl;
	      sigzz += zeta*phi[i][j][k]* Qzzl;
	      }	      

	      /* anti-symmetric part of stress tensor */
	      dbdtauxb =  dy_(tauxy,i,j,k) + dz_(tauxz,i,j,k);
	      dbdtauyb = -dx_(tauxy,i,j,k) + dz_(tauyz,i,j,k);
	      dbdtauzb = -dy_(tauyz,i,j,k) - dx_(tauxz,i,j,k);

	      /* thermodynamic stress go- */
	      /*	      Force[0] = dx_(sxx,i,j,k) + dy_(sxy,i,j,k) + dz_(sxz,i,j,k) - friction*ux[i][j][k];
	      Force[1] = dx_(sxy,i,j,k) + dy_(syy,i,j,k) + dz_(syz,i,j,k) - friction*uy[i][j][k];
	      Force[2] = dx_(sxz,i,j,k) + dy_(syz,i,j,k) + dz_(szz,i,j,k) - friction*uz[i][j][k];*/

	      Force[0]=Force[1]=Force[2]=0.0;
	      /*Fh[i][j][k][0]=Fh[i][j][k][1]=Fh[i][j][k][2]=0.0;
	      sigxx=sigxy=sigxz=sigyy=sigyz=sigzz=0.0;*/

	      Force[0] = -phil*dx_(mu,i,j,k)- globalP[0] -friction*ux[i][j][k];
	      Force[1] = -phil*dy_(mu,i,j,k)- globalP[1] -friction*uy[i][j][k];
	      Force[2] = -phil*dz_(mu,i,j,k)- globalP[2] -friction*uz[i][j][k];

	      //Force[1] += bodyforce;

	      A2= (rho*temperature+0.0*phivr*pG[i][j][k])/10.0;
	      A1= A2;
	      A0= rho-14.0*A2;
	      B2= rho/24.0;
	      B1= 8.0*B2;
	      C2= -rho/24.0;
	      C1= 2.0*C2;
	      C0= -2.0*rho/3.0;
	      D2= rho/16.0;
	      D1= 8.0*D2;
	      G2xx= phivr*sigxx/16.0;
	      G2yy= phivr*sigyy/16.0;
	      G2zz= phivr*sigzz/16.0;
	      G2xy= phivr*sigxy/16.0;
	      G2xz= phivr*sigxz/16.0;
	      G2yz= phivr*sigyz/16.0;
	      G1xx= 8.0*G2xx;
	      G1yy= 8.0*G2yy;
	      G1zz= 8.0*G2zz;
	      
	      usq = ux[i][j][k]*ux[i][j][k] + uy[i][j][k]*uy[i][j][k] + uz[i][j][k]*uz[i][j][k];

	      /* i = 0 */
	      //A0 -= densitycorrection;

	      feq[i][j][k][0] = A0 + C0*usq;
	      
	      /* i = 1, 2, 3, 4, 5, 6 */
	      for (l = 1; l <= 6; l++) {
		udote = ux[i][j][k]*e[l][0] + uy[i][j][k]*e[l][1] + uz[i][j][k]*e[l][2];
		omdote  = Force[0]*e[l][0] 	 + Force[1]*e[l][1] 	  + Force[2]*e[l][2];
		omdote += dbdtauxb*e[l][0] 	 + dbdtauyb*e[l][1] 	  + dbdtauzb*e[l][2];
		omdote += Fh[i][j][k][0]*e[l][0] + Fh[i][j][k][1]*e[l][1] + Fh[i][j][k][2]*e[l][2]; 

		feq[i][j][k][l] = A1 + B1*udote + C1*usq + D1*udote*udote
		  + G1xx*e[l][0]*e[l][0] + G1yy*e[l][1]*e[l][1] + 
		  G1zz*e[l][2]*e[l][2]+ tau1*omdote/3.0;
				}

	      /* i = 7, 8, 9, 10, 11, 12, 13, 14 */
	      for (l = 7; l <= 14; l++) {
		udote = ux[i][j][k]*e[l][0] + uy[i][j][k]*e[l][1] + uz[i][j][k]*e[l][2];
		
		omdote  = Force[0]*e[l][0] 	 + Force[1]*e[l][1] 	  + Force[2]*e[l][2];
		omdote += dbdtauxb*e[l][0] 	 + dbdtauyb*e[l][1] 	  + dbdtauzb*e[l][2];
		omdote += Fh[i][j][k][0]*e[l][0] + Fh[i][j][k][1]*e[l][1] + Fh[i][j][k][2]*e[l][2];
		
		feq[i][j][k][l] = A2 + B2*udote + C2*usq + D2*udote*udote
		  + G2xx*e[l][0]*e[l][0] + G2yy*e[l][1]*e[l][1] 
		  + G2zz*e[l][2]*e[l][2] + 2.0*G2xy*e[l][0]*e[l][1] 
		  + 2.0*G2xz*e[l][0]*e[l][2] + 2.0*G2yz*e[l][1]*e[l][2]
		  + tau1*omdote/24.0;
				}
			}
		}
	}
}

/* calculate observables */
void observables(double f[Lx][Ly][Lz][15])
{
	int i, j, k, l;
	int iup,idwn;
	int jup,jdwn;
	int kup,kdwn;
	double phi_total;
	int R;

	/* calculate local density and velocity */
	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				density[i][j][k] = 0.0;

				ux[i][j][k] = 0.0;
				uy[i][j][k] = 0.0;
				uz[i][j][k] = 0.0;

				for (l = 0; l < 15; l++) {
					density[i][j][k] += f[i][j][k][l];

					ux[i][j][k] += f[i][j][k][l]*e[l][0];
					uy[i][j][k] += f[i][j][k][l]*e[l][1];
					uz[i][j][k] += f[i][j][k][l]*e[l][2];
				}
				ux[i][j][k] = ux[i][j][k]/density[i][j][k];
				uy[i][j][k] = uy[i][j][k]/density[i][j][k];
				uz[i][j][k] = uz[i][j][k]/density[i][j][k];
			}
		}
	}
	phi_total = 0.0;

	Vx = 0.0;
	Vy = 0.0;
    Vz = 0.0;

	R = 0;
	R1 = 0;
    	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				phi_total += phi[i][j][k]; 

	 	                Vx += phi[i][j][k]*(ux[i][j][k] + w1*Px[i][j][k]);
	        	        Vy += phi[i][j][k]*(uy[i][j][k] + w1*Py[i][j][k]);
	        	        Vz += phi[i][j][k]*(uz[i][j][k] + w1*Pz[i][j][k]);

				if (phi[i][j][k] > 1.0) {R += 1;}
				if (R > R1) {R1 = R;}
			}
			R = 0;
		}
    	}
	R1 = R1/2;

    	Vx = Vx/phi_total;
    	Vy = Vy/phi_total;
	Vz = Vz/phi_total;
}

/*
 *	Binary Fluid
 */

/* initialization */
void initializephi(void)
{
	int i, j, k;
	int R,seed;
    
    seed=6452;
    srand48(seed);

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0 ; k < Lz; k++) {

#if DROPLET
				R = sqrt((i-Lx/2)*(i-Lx/2) + (j-Ly/2)*(j-Ly/2) + (k-Lz/2)*(k-Lz/2));
				if (R < radius) {
					phi[i][j][k] = 2.0;
				}
				else {
					phi[i][j][k] = 0.0;
				}
#endif
                phi[i][j][k] = phi0;//+0.1*(1.0-2.0*drand48());
			}
		}
	}
}

/* calculate chemical potential */
void phi_function(void)
{
	int i, j, k;
	int iup, jup, kup;
	int iup2, jup2, kup2;
	int idwn, jdwn, kdwn;
	int idwn2, jdwn2, kdwn2;

	double dphidx,dphidy,dphidz;
	double gradphisq,phil,f;
	double vx, vy, vz;
	double vxdphidx, vydphidy, vzdphidz;
	double vxdmdx,   vydmdy,   vzdmdz;

	double Qsqxx,Qsqxy,Qsqxz,Qsqyy,Qsqyz,Qsqzz;
	double Qxxl,Qxyl,Qxzl,Qyyl,Qyzl,Qzzl;
	double TrQ2;

	double Psquare,gradPsq;
	double divu, divP;

	double Jrx[Lx][Ly][Lz], Jry[Lx][Ly][Lz], Jrz[Lx][Ly][Lz];
	double Jax[Lx][Ly][Lz], Jay[Lx][Ly][Lz], Jaz[Lx][Ly][Lz];
	double Jdx[Lx][Ly][Lz], Jdy[Lx][Ly][Lz], Jdz[Lx][Ly][Lz];
	double divJa, divJd;


	/* find chemical potential */
	phi_average = 0.0;
	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
			  

			  Qgradphix[i][j][k] = Qxx[i][j][k]*dx_(phi,i,j,k)+Qxy[i][j][k]*dy_(phi,i,j,k)+Qxz[i][j][k]*dz_(phi,i,j,k);
			  Qgradphiy[i][j][k] = Qxy[i][j][k]*dx_(phi,i,j,k)+Qyy[i][j][k]*dy_(phi,i,j,k)+Qyz[i][j][k]*dz_(phi,i,j,k);
			  Qgradphiz[i][j][k] = Qxz[i][j][k]*dx_(phi,i,j,k)+Qyz[i][j][k]*dy_(phi,i,j,k)+(-Qxx[i][j][k]-Qyy[i][j][k])*dz_(phi,i,j,k);
              laplacianphi[i][j][k] = laplacian_(phi,i,j,k);


	      Jrx[i][j][k]=Jry[i][j][k]=Jrz[i][j][k]=0.0;

	      if(phi[i][j][k]>=0.0){
		Jrx[i][j][k]=sqrt(noisestrength*phi[i][j][k])*sqrt(3.0)*(1.0-2.0*drand48());
		Jry[i][j][k]=sqrt(noisestrength*phi[i][j][k])*sqrt(3.0)*(1.0-2.0*drand48());
		Jrz[i][j][k]=sqrt(noisestrength*phi[i][j][k])*sqrt(3.0)*(1.0-2.0*drand48());
	      }

			}
		}
	}

    //use the bit below for neutral wetting BC
#if BC
    for (i = 0; i < Lx; i++) {
        for (j = 0; j < Ly; j++) {
            laplacianphi[i][j][0] = laplacianphi[i][j][1];
            laplacianphi[i][j][Lz-1]= laplacianphi[i][j][Lz-2];
        }
    }
#endif
    
	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {

				Qxxl=Qxx[i][j][k];
				Qxyl=Qxy[i][j][k];
				Qyyl=Qyy[i][j][k];
				Qxzl=Qxz[i][j][k];
				Qyzl=Qyz[i][j][k];
				Qzzl=-Qxxl-Qyyl;
				
				Qsqxx=Qxxl+Qyyl;
				Qsqzz=Qsqxx*Qsqxx+Qxzl*Qxzl+Qyzl*Qyzl;
				Qsqxy=Qxyl*Qsqxx+Qxzl*Qyzl;
				Qsqxz=Qxxl*Qxzl-Qxzl*Qsqxx+Qxyl*Qyzl;
				Qsqyz=Qxyl*Qxzl-Qyzl*Qsqxx+Qyyl*Qyzl;
				Qsqxx=Qxxl*Qxxl+Qxyl*Qxyl+Qxzl*Qxzl;
				Qsqyy=Qyyl*Qyyl+Qxyl*Qxyl+Qyzl*Qyzl;
				
				TrQ2=Qsqxx+Qsqyy+Qsqzz;
                
                mu[i][j][k] = alpha2phi * phi[i][j][k]+
		              alpha4phi * phi[i][j][k]*phi[i][j][k]*phi[i][j][k]
				  - Kphi*laplacian_(phi,i,j,k)
				  +delta*(- Abulk/6.0*TrQ2
				  - Abulk/3.0*(Qsqxx*Qxxl+2.0*Qsqxy*Qxyl+2.0*Qsqxz*Qxzl+Qsqyy*Qyyl+2.0*Qsqyz*Qyzl+Qsqzz*Qzzl )
					+ Abulk/4.0*TrQ2*TrQ2);
				
				//W contribution
				
				mu[i][j][k] += -2.0*W*(dx_(Qgradphix,i,j,k)+dy_(Qgradphiy,i,j,k)+dz_(Qgradphiz,i,j,k));

				phi_average = phi_average + phi[i][j][k];
                
			}
		}
	}
	phi_average = phi_average/(Lx*Ly*Lz);

    /* boundary condition for mu */
#if BC
    for (i = 0; i < Lx; i++) {
        for (j = 0; j < Ly; j++) {
            for (k = 0; k < Lz; k++) {
                
                if(k == 0) {
                    mu[i][j][k] = mu[i][j][k+1];
                }
                if(k == Lz-1) {
                    mu[i][j][k] = mu[i][j][k-1];
                }
            }
        }
    }
#endif
    
	/* calculate diffusive and advective current go- */
	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				phil = phi[i][j][k];

				Jax[i][j][k] = phi[i][j][k]*(ux[i][j][k] + w1*Px[i][j][k]);
				Jay[i][j][k] = phi[i][j][k]*(uy[i][j][k] + w1*Py[i][j][k]);
				Jaz[i][j][k] = phi[i][j][k]*(uz[i][j][k] + w1*Pz[i][j][k]);

				Jdx[i][j][k] = -M*dx_(mu,i,j,k);
				Jdy[i][j][k] = -M*dy_(mu,i,j,k);
				Jdz[i][j][k] = -M*dz_(mu,i,j,k);
			}
		}
	}
	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) { 				
				iup  = iupa[i];
				iup2 = iupa[iupa[i]];

				jup  = jupa[j];
				jup2 = jupa[jupa[j]];

				kup  = kupa[k];
				kup2 = kupa[kupa[k]];

				idwn  = idwna[i];
				idwn2 = idwna[idwna[i]];

				jdwn  = jdwna[j];
				jdwn2 = jdwna[jdwna[j]];

				kdwn  = kdwna[k];
				kdwn2 = kdwna[kdwna[k]];

				vx = ux[i][j][k] + w1*Px[i][j][k];
				vy = uy[i][j][k] + w1*Py[i][j][k];
				vz = uz[i][j][k] + w1*Pz[i][j][k];

				phil = phi[i][j][k];

				/* upwind scheme */
				vxdphidx = max(vx,0.0)*(2.0*phi[iup][j][k]  + 3.0*phi[i][j][k]   - 6.0*phi[idwn][j][k] +     phi[idwn2][j][k])/6.0 +
					   min(vx,0.0)*(   -phi[iup2][j][k] + 6.0*phi[iup][j][k] - 3.0*phi[i][j][k]    - 2.0*phi[idwn][j][k]) /6.0;
				vydphidy = max(vy,0.0)*(2.0*phi[i][jup][k]  + 3.0*phi[i][j][k]   - 6.0*phi[i][jdwn][k] +     phi[i][jdwn2][k])/6.0 +
					   min(vy,0.0)*(   -phi[i][jup2][k] + 6.0*phi[i][jup][k] - 3.0*phi[i][j][k]    - 2.0*phi[i][jdwn][k]) /6.0;
				vzdphidz = max(vz,0.0)*(2.0*phi[i][j][kup]  + 3.0*phi[i][j][k]   - 6.0*phi[i][j][kdwn] +     phi[i][j][kdwn2])/6.0 +
					   min(vz,0.0)*(   -phi[i][j][kup2] + 6.0*phi[i][j][kup] - 3.0*phi[i][j][k]    - 2.0*phi[i][j][kdwn]) /6.0;

				divu = dx_(ux,i,j,k) + dy_(uy,i,j,k) + dz_(uz,i,j,k);
				divP = dx_(Px,i,j,k) + dy_(Py,i,j,k) + dz_(Pz,i,j,k);

				h_phi[i][j][k] = -vxdphidx - vydphidy - vzdphidz - phil*(divu + w1*divP) + M*laplacian_(mu,i,j,k);
				h_phi[i][j][k] += dx_(Jrx,i,j,k)+dy_(Jry,i,j,k)+dz_(Jrz,i,j,k);


			        h_phi[i][j][k] = 0.0; /*this means no density variation/update */

			}
		}
	}
}

/* predictor-corrector for phi */
void updatephi0(double phipr[Lx][Ly][Lz], double phiold[Lx][Ly][Lz])
{
	int i, j, k;
	double dphidz;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				phipr[i][j][k] = phiold[i][j][k] + dt*h_phi[i][j][k];

				h_phiold[i][j][k] = h_phi[i][j][k];
			}
		}
	}
    
    /* boundary condition for phi - this may need cheching ! */
    for (i = 0; i < Lx; i++) {
        for (j = 0; j < Ly; j++) {
#if BC
            
            phipr[i][j][0] = phipr[i][j][1];
            phipr[i][j][Lz-1] = phipr[i][j][Lz-2];
#endif
        }
    }
    
}
void updatephi(double phinew[Lx][Ly][Lz], double phiold[Lx][Ly][Lz])
{
	int i, j, k;
	double dphidz;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				phinew[i][j][k] = phiold[i][j][k] + 0.5*dt*(h_phi[i][j][k] + h_phiold[i][j][k]);
			}
		}
	}
    
    /* boundary condition for phi - this may need cheching ! */
    for (i = 0; i < Lx; i++) {
        for (j = 0; j < Ly; j++) {
#if BC
            
            phinew[i][j][0] = phinew[i][j][1];
            phinew[i][j][Lz-1] = phinew[i][j][Lz-2];
#endif
        }
    }
    
}

/*
 *	Leslie-Ericksen
 */

/* Initialization */
void initializeQ(void)
{
	int i, j, k;
	int R;
	
	double phase,phase2;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
#if DROPLET
				R = sqrt((i-Lx/2)*(i-Lx/2) + (j-Ly/2)*(j-Ly/2) + (k-Lz/2)*(k-Lz/2));

				if (R < radius) {

					Qxx[i][j][k]=0.3*(1.0-2.0*drand48());
					Qxy[i][j][k]=0.3*(1.0-2.0*drand48());
					Qxz[i][j][k]=0.3*(1.0-2.0*drand48());
					Qyy[i][j][k]=0.3*(1.0-2.0*drand48());
					Qyz[i][j][k]=0.3*(1.0-2.0*drand48());

#if TWIST
	       phase= Pi/180.0*(angxytop-(angxytop-angxybot)/(Lz-1.0)*k);
	       phase2= Pi/180.0*(angztop-(angztop-(angzbot))/(Lz-1.0)*k);
	       
	       amplitude=0.3;

	       Qxx[i][j][k]= amplitude*
		 (3.0/2.0*sin(phase2)*sin(phase2)*cos(phase)*cos(phase)-1.0/2.0);

	       Qxy[i][j][k]= 3.0*amplitude/2.0*
		 (sin(phase2)*sin(phase2)*cos(phase)*sin(phase));

	       Qyy[i][j][k]= amplitude*
		 (3.0/2.0*sin(phase2)*sin(phase2)*sin(phase)*sin(phase)-1.0/2.0);

	       Qxz[i][j][k]= 3.0*amplitude/2.0*(sin(phase2)*cos(phase2)*cos(phase));

	       Qyz[i][j][k]= 3.0*amplitude/2.0*(sin(phase2)*cos(phase2)*sin(phase));


#endif

				}
				else {
					Qxx[i][j][k]=Qxy[i][j][k]=Qxz[i][j][k]=Qyy[i][j][k]=Qyz[i][j][k]=0.0001;
				}
#endif

					Qxx[i][j][k]=0.3*(1.0-2.0*drand48());
			 		Qxy[i][j][k]=0.3*(1.0-2.0*drand48());
					Qxz[i][j][k]=0.3*(1.0-2.0*drand48());
					Qyy[i][j][k]=0.3*(1.0-2.0*drand48());
					Qyz[i][j][k]=0.3*(1.0-2.0*drand48());

			}
		}
	}
}

/* calculate molecular field for Q tensor*/
void Q_function(void)
{
	int i, j, k;
	int iup, jup, kup;
	int iup2, jup2, kup2;
	int idwn, jdwn, kdwn;
	int idwn2, jdwn2, kdwn2;

	double dQxxdx,dQxxdy,dQxxdz,dQxydx,dQxydy,dQxydz,dQyydx,dQyydy,dQyydz;
	double dQxzdx,dQxzdy,dQxzdz,dQyzdx,dQyzdy,dQyzdz;
	double dQyxdx,dQyxdy,dQyxdz;
	double dQzxdx,dQzxdy,dQzxdz,dQzydx,dQzydy,dQzydz;
	double dQzzdx,dQzzdy,dQzzdz;
	double trt,trd2Q;
	double d2Qxxdxdx,d2Qxxdydy,d2Qxxdxdy,d2Qxxdzdz,d2Qxxdxdz,d2Qxxdydz;
	double d2Qyydxdx,d2Qyydydy,d2Qyydxdy,d2Qyydzdz,d2Qyydxdz,d2Qyydydz;
	double d2Qxydxdx,d2Qxydydy,d2Qxydxdy,d2Qxydzdz,d2Qxydxdz,d2Qxydydz;
	double d2Qxzdxdx,d2Qxzdydy,d2Qxzdxdy,d2Qxzdzdz,d2Qxzdxdz,d2Qxzdydz;
	double d2Qyzdxdx,d2Qyzdydy,d2Qyzdxdy,d2Qyzdzdz,d2Qyzdxdz,d2Qyzdydz;
	double DGxx,DGyy,DGxy,DGyx,DGzz,DGxz,DGzx,DGyz,DGzy,TrG,divQx,divQy,divQz;
	
	double Qsqxx,Qsqxy,Qsqxz,Qsqyy,Qsqyz,Qsqzz,Qxxl,Qxyl,Qxzl,Qyyl,Qyzl,Qzzl;
	double Hxx,Hyy,Hxy,Hxz,Hyz,TrQ2,TrDQI;
	double mDQ4xx,mDQ4xy,mDQ4yy,mDQ4xz,mDQ4yz,mDQ4zz,nnxxl,nnyyl;

	double dphidx, dphidy, dphidz;
	double phil, gradphisq, freeE;

	double duxdx, duxdy, duxdz;
	double duydx, duydy, duydz;
	double duzdx, duzdy, duzdz;
	double divu, divP;

	double vx, vy, vz;
	double vxdQxxdx, vydQxxdy, vzdQxxdz;
	double vxdQxydx, vydQxydy, vzdQxydz;
	double vxdQxzdx, vydQxzdy, vzdQxzdz;
	double vxdQyydx, vydQyydy, vzdQyydz;
	double vxdQyzdx, vydQyzdy, vzdQyzdz;

    double Dyy,Dyz,Dzz,omegayz,DD,omega,omega2;

	energy = 0.0;
	for (i = 0; i < Lx; i++) {
	  for (j = 0; j < Ly; j++) {
	    for (k = 0; k < Lz; k++) {
	      iup  = iupa[i];
	      iup2 = iupa[iupa[i]];

	      jup  = jupa[j];
	      jup2 = jupa[jupa[j]];

	      kup  = kupa[k];
	      kup2 = kupa[kupa[k]];

	      idwn  = idwna[i];
	      idwn2 = idwna[idwna[i]];
	      
	      jdwn  = jdwna[j];
	      jdwn2 = jdwna[jdwna[j]];
	      
	      kdwn  = kdwna[k];
	      kdwn2 = kdwna[kdwna[k]];

	      dQxxdx=dx_(Qxx,i,j,k);
	      dQxydx=dx_(Qxy,i,j,k);
	      dQxzdx=dx_(Qxz,i,j,k);
	      dQyxdx=dQxydx;
	      dQyydx=dx_(Qyy,i,j,k);
	      dQyzdx=dx_(Qyz,i,j,k);
	      dQzxdx=dQxzdx;
	      dQzydx=dQyzdx;
	      dQzzdx=-(dQxxdx+dQyydx);

	      dQxxdy=dy_(Qxx,i,j,k);
	      dQxydy=dy_(Qxy,i,j,k);
	      dQxzdy=dy_(Qxz,i,j,k);
	      dQyxdy=dQxydy;
	      dQyydy=dy_(Qyy,i,j,k);
	      dQyzdy=dy_(Qyz,i,j,k);
	      dQzxdy=dQxzdy;
	      dQzydy=dQyzdy;
	      dQzzdy=-(dQxxdy+dQyydy);

	      dQxxdz=dz_(Qxx,i,j,k);
	      dQxydz=dz_(Qxy,i,j,k);
	      dQxzdz=dz_(Qxz,i,j,k);	     
	      dQyxdz=dQxydz;
	      dQyydz=dz_(Qyy,i,j,k);
	      dQyzdz=dz_(Qyz,i,j,k); 
	      dQzxdz=dQxzdz;
	      dQzydz=dQyzdz;
	      dQzzdz=-(dQxxdz+dQyydz); 

	      Qxxl=Qxx[i][j][k];
	      Qxyl=Qxy[i][j][k];
	      Qyyl=Qyy[i][j][k];
	      Qxzl=Qxz[i][j][k];
	      Qyzl=Qyz[i][j][k];
	      Qzzl=-Qxxl-Qyyl;

	      nnxxl= Qxxl+1.0/3.0;
	      nnyyl= Qyyl+1.0/3.0;

	      dphidx=dx_(phi,i,j,k);
	      dphidy=dy_(phi,i,j,k);
	      dphidz=dz_(phi,i,j,k);

	      gradphisq = dphidx*dphidx + dphidy*dphidy + dphidz*dphidz;

	      phil=phi[i][j][k];

	      Qsqxx=Qxxl+Qyyl;
	      Qsqzz=Qsqxx*Qsqxx+Qxzl*Qxzl+Qyzl*Qyzl;
	      Qsqxy=Qxyl*Qsqxx+Qxzl*Qyzl;
	      Qsqxz=Qxxl*Qxzl-Qxzl*Qsqxx+Qxyl*Qyzl;
	      Qsqyz=Qxyl*Qxzl-Qyzl*Qsqxx+Qyyl*Qyzl;
	      Qsqxx=Qxxl*Qxxl+Qxyl*Qxyl+Qxzl*Qxzl;
	      Qsqyy=Qyyl*Qyyl+Qxyl*Qxyl+Qyzl*Qyzl;

	      TrQ2=Qsqxx+Qsqyy+Qsqzz;

	      txx=txy=txz=tyy=tyx=tyz=tzz=tzy=tzx=0;

	      trt=txx+tyy+tzz;

	      /* second derivative contribution to molecular field Hab */

	      DEHxx[i][j][k] = L1 * laplacian_(Qxx,i,j,k) + txx - 1.0/3.0 * trt;
	      DEHxy[i][j][k] = L1 * laplacian_(Qxy,i,j,k) + (txy+tyx)/2.0;
	      DEHxz[i][j][k] = L1 * laplacian_(Qxz,i,j,k) + (txz+tzx)/2.0;
	      DEHyy[i][j][k] = L1 * laplacian_(Qyy,i,j,k) + tyy - 1.0/3.0 * trt;
	      DEHyz[i][j][k] = L1 * laplacian_(Qyz,i,j,k) + (tyz+tzy)/2.0;

	      DEHxx[i][j][k] -= W * (dphidx*dphidx - 1.0/3.0 * gradphisq);
	      DEHxy[i][j][k] -= W * (dphidx*dphidy);
	      DEHxz[i][j][k] -= W * (dphidx*dphidz);
	      DEHyy[i][j][k] -= W * (dphidy*dphidy - 1.0/3.0 * gradphisq);
	      DEHyz[i][j][k] -= W * (dphidy*dphidz);

	      /* bulk contribution to molecular field Hab */

	      gammac=gamma+delta*phi[i][j][k];

	      DEHxx[i][j][k]+= Abulk*(-(1.0-gammac/3.0)*Qxxl+gammac*(Qsqxx-TrQ2/3.0)-gammac*Qxxl*TrQ2);
	      DEHxy[i][j][k]+= Abulk*(-(1.0-gammac/3.0)*Qxyl+gammac*Qsqxy-gammac*Qxyl*TrQ2);
	      DEHyy[i][j][k]+= Abulk*(-(1.0-gammac/3.0)*Qyyl+gammac*(Qsqyy-TrQ2/3.0)-gammac*Qyyl*TrQ2);
	      DEHxz[i][j][k]+= Abulk*(-(1.0-gammac/3.0)*Qxzl+gammac*Qsqxz-gammac*Qxzl*TrQ2);
	      DEHyz[i][j][k]+= Abulk*(-(1.0-gammac/3.0)*Qyzl+gammac*Qsqyz-gammac*Qyzl*TrQ2);

	      //	      if(i==(Lx/2) && (j==Ly/2) && (k==Lz/2)) fprintf(stderr,"%E %E %E %E %E %E\n",Qxx[i][j][k],Qxy[i][j][k],Qxz[i][j][k],Qyy[i][j][k],Qyz[i][j][k],phi[i][j][k]);

	      Hxx=DEHxx[i][j][k];
	      Hxy=DEHxy[i][j][k];
	      Hxz=DEHxz[i][j][k];
	      Hyy=DEHyy[i][j][k];
	      Hyz=DEHyz[i][j][k];
	      
	      DEH2xx[i][j][k]=DEHxx[i][j][k];
	      DEH2xy[i][j][k]=DEHxy[i][j][k];
	      DEH2xz[i][j][k]=DEHxz[i][j][k];
	      DEH2yy[i][j][k]=DEHyy[i][j][k];
	      DEH2yz[i][j][k]=DEHyz[i][j][k];

	     /*
	      *	elastic energy
	      */
	      energy += dQxxdx*dQxxdx + dQxxdy*dQxxdy + dQxxdz*dQxxdz;
	      energy += dQxydx*dQxydx + dQxydy*dQxydy + dQxydz*dQxydz;
	      energy += dQxzdx*dQxzdx + dQxzdy*dQxzdy + dQxzdz*dQxzdz;
	      energy += dQyydx*dQyydx + dQyydy*dQyydy + dQyydz*dQyydz;
	      energy += dQyzdx*dQyzdx + dQyzdy*dQyzdy + dQyzdz*dQyzdz;

	      /* compute (minus) stress tensor */

	      DGxx=(dQxxdx*dQxxdx+2.0*dQxydx*dQxydx+dQyydx*dQyydx+
		    2.0*dQxzdx*dQxzdx+2.0*dQyzdx*dQyzdx+
		    (dQxxdx+dQyydx)*(dQxxdx+dQyydx));
	      DGyy= (dQxxdy*dQxxdy+2.0*dQxydy*dQxydy+dQyydy*dQyydy+
		     2.0*dQxzdy*dQxzdy+2.0*dQyzdy*dQyzdy+
		     (dQxxdy+dQyydy)*(dQxxdy+dQyydy));
	      DGzz= (dQxxdz*dQxxdz+2.0*dQxydz*dQxydz+dQyydz*dQyydz+
		     2.0*dQxzdz*dQxzdz+2.0*dQyzdz*dQyzdz+
		     (dQxxdz+dQyydz)*(dQxxdz+dQyydz));
	      DGxy= (dQxxdx*dQxxdy+2.0*dQxydx*dQxydy+dQyydx*dQyydy+
		     2.0*dQxzdx*dQxzdy+2.0*dQyzdx*dQyzdy+
		     (dQxxdx+dQyydx)*(dQxxdy+dQyydy));
	      DGxz= (dQxxdx*dQxxdz+2.0*dQxydx*dQxydz+dQyydx*dQyydz+
		     2.0*dQxzdx*dQxzdz+2.0*dQyzdx*dQyzdz+
		     (dQxxdx+dQyydx)*(dQxxdz+dQyydz));
	      DGyz= (dQxxdy*dQxxdz+2.0*dQxydy*dQxydz+dQyydy*dQyydz+
		     2.0*dQxzdy*dQxzdz+2.0*dQyzdy*dQyzdz+
		     (dQxxdy+dQyydy)*(dQxxdz+dQyydz));

	      gammac=gamma+delta*phil;

	      freeE = 0.25*alpha2phi*phil*phil*(phil-2.0)*(phil-2.0) + 0.5*Kphi*gradphisq;
	      freeE += Abulk/2.0*(1.0-gammac/3.0)*TrQ2;
	      freeE += Abulk*gammac/4.0*TrQ2*TrQ2;
	      freeE += -Abulk*gammac/3.0*(Qsqxx*Qxxl+2.0*Qsqxy*Qxyl+2.0*Qsqxz*Qxzl+Qsqyy*Qyyl+2.0*Qsqyz*Qyzl+Qsqzz*Qzzl );
	      freeE += L1/2.0*(dQxxdx*dQxxdx + dQxxdy*dQxxdy + dQxxdz*dQxxdz);
	      freeE += L1/2.0*(dQxydx*dQxydx + dQxydy*dQxydy + dQxydz*dQxydz);
	      freeE += L1/2.0*(dQxzdx*dQxzdx + dQxzdy*dQxzdy + dQxzdz*dQxzdz);
	      freeE += L1/2.0*(dQyydx*dQyydx + dQyydy*dQyydy + dQyydz*dQyydz);
	      freeE += L1/2.0*(dQyzdx*dQyzdx + dQyzdy*dQyzdy + dQyzdz*dQyzdz);

	      DG2xx[i][j][k]=-L1*DGxx;
	      DG2yy[i][j][k]=-L1*DGyy;
	      DG2zz[i][j][k]=-L1*DGzz;
	      DG2xy[i][j][k]=-L1*DGxy;
	      DG2xz[i][j][k]=-L1*DGxz;
	      DG2yz[i][j][k]=-L1*DGyz;
	      DGyx=L1*DGxy;
	      DGzx=L1*DGxz;
	      DGzy=L1*DGyz;

	      DG2xx[i][j][k] += -Kphi*dphidx*dphidx - phi[i][j][k]*mu[i][j][k] + freeE;
	      DG2xy[i][j][k] += -Kphi*dphidx*dphidy;
	      DG2xz[i][j][k] += -Kphi*dphidx*dphidz;
	      DG2yy[i][j][k] += -Kphi*dphidy*dphidy - phi[i][j][k]*mu[i][j][k] + freeE;
	      DG2yz[i][j][k] += -Kphi*dphidy*dphidz;
	      DG2zz[i][j][k] += -Kphi*dphidz*dphidz - phi[i][j][k]*mu[i][j][k] + freeE;

	      sxx[i][j][k] = -Kphi*dphidx*dphidx - phi[i][j][k]*mu[i][j][k] + freeE;
	      sxy[i][j][k] = -Kphi*dphidx*dphidy;
	      sxz[i][j][k] = -Kphi*dphidx*dphidz;
	      syy[i][j][k] = -Kphi*dphidy*dphidy - phi[i][j][k]*mu[i][j][k] + freeE;
	      syz[i][j][k] = -Kphi*dphidy*dphidz;
	      szz[i][j][k] = -Kphi*dphidz*dphidz - phi[i][j][k]*mu[i][j][k] + freeE;

	      // needs W term in stress sxx[i][j][k] -= 2.0*W*(Qxx*dphidx*dphidx+Qxyl*dphidx*dphidy+Qxzl*dphidx*dphidz);

	      /* this will go to the bodyforce */
			
	      Fh[i][j][k][0]=-(Hxx*dQxxdx+2.0*Hxy*dQxydx
		+2.0*Hxz*dQxzdx+Hyy*dQyydx+2.0*Hyz*dQyzdx
	        +(-Hyy-Hxx)*(-dQxxdx-dQyydx));
	      Fh[i][j][k][1]=-(Hxx*dQxxdy+2.0*Hxy*dQxydy
		+2.0*Hxz*dQxzdy+Hyy*dQyydy+2.0*Hyz*dQyzdy
		+(-Hyy-Hxx)*(-dQxxdy-dQyydy));
	      Fh[i][j][k][2]=-(Hxx*dQxxdz+2.0*Hxy*dQxydz
		+2.0*Hxz*dQxzdz+Hyy*dQyydz+2.0*Hyz*dQyzdz
		+(-Hyy-Hxx)*(-dQxxdz-dQyydz));

/* work out antisymmetric part of the stress tensor, this will go to the bodyforce */

	      tauxy[i][j][k]= -phivr*
		(Qxy[i][j][k]*(DEHxx[i][j][k]-DEHyy[i][j][k])-
		 DEHxy[i][j][k]*(Qxx[i][j][k]-Qyy[i][j][k])+
		 DEHxz[i][j][k]*Qyz[i][j][k]-DEHyz[i][j][k]*Qxz[i][j][k])
		+0.0*phivr*(DG2xy[i][j][k]-DGyx)/2.0;
	      tauxz[i][j][k]= -phivr*
		(Qxz[i][j][k]*(2.0*DEHxx[i][j][k]+DEHyy[i][j][k])-
		 DEHxz[i][j][k]*(2.0*Qxx[i][j][k]+Qyy[i][j][k])+
		 DEHxy[i][j][k]*Qyz[i][j][k]-DEHyz[i][j][k]*Qxy[i][j][k])
		+0.0*phivr*(DG2xz[i][j][k]-DGzx)/2.0;
	      tauyz[i][j][k]= -phivr*
		(Qyz[i][j][k]*(2.0*DEHyy[i][j][k]+DEHxx[i][j][k])-
		 DEHyz[i][j][k]*(2.0*Qyy[i][j][k]+Qxx[i][j][k])+
		 DEHxy[i][j][k]*Qxz[i][j][k]-DEHxz[i][j][k]*Qxy[i][j][k])
		+0.0*phivr*(DG2yz[i][j][k]-DGzy)/2.0;

	      duxdx = dx_(ux,i,j,k);
	      duxdy = dy_(ux,i,j,k);
	      duxdz = dz_(ux,i,j,k);
				
	      duydx = dx_(uy,i,j,k);
	      duydy = dy_(uy,i,j,k);
	      duydz = dz_(uy,i,j,k);

	      duzdx = dx_(uz,i,j,k);
	      duzdy = dy_(uz,i,j,k);
	      duzdz = dz_(uz,i,j,k);
            
          Dyy=duydy;
          Dyz=0.5*(duydz+duzdy);
          Dzz=duzdz;
          omegayz=0.5*(duzdy-duydz);
          DD=0.5*(Dyy*Dyy+2.0*Dyz*Dyz+Dzz*Dzz);
          DD=sqrt(DD);
          omega2=0.5*(2.0*omegayz*omegayz);
          omega=sqrt(omega2);
            
	      //flowparameter[i][j][k]=0.0;

	      //if((DD+omega2)>0.0001) 

	      flowparameter[i][j][k]=(DD-omega)/(DD+omega);

	      compressibility[i][j][k]=(duydy+duzdz);

	      //	      if(flowparameter[i][j][k]>100) fprintf(stderr,"%lf %lf %lf %lf %lf %lf %lf\n",DD,omega2,duydy,duydy,duydz,duzdz,flowparameter[i][j][k]);

	      /* -(Q+1/3) Tr(D.(Q+1/3)) term*/
	      TrDQI= -(duxdx+duydy+duzdz)/3.0-
		(Qxxl*duxdx+Qyyl*duydy-(Qxxl+Qyyl)*duzdz+Qxyl*(duxdy+duydx)+
		 Qxzl*(duxdz+duzdx)+Qyzl*(duydz+duzdy));
	      mDQ4xy=Qxyl*TrDQI;
	      mDQ4yy=nnyyl*TrDQI;
	      mDQ4xz=Qxzl*TrDQI; 
	      mDQ4yz=Qyzl*TrDQI;
	      mDQ4zz= (-Qxxl-Qyyl+1.0/3.0)*TrDQI;
	      mDQ4xx=nnxxl*TrDQI;

	      /*add molecular field times Gamma and S term */

	      Hxx=G*Hxx+duxdy*Qxyl*(1.0+xi)+duxdx*2.0*nnxxl*xi+
		duydx*Qxyl*(xi-1.0)+duxdz*Qxzl*(1.0+xi)+duzdx*Qxzl*(xi-1.0)+
		2.0*xi*mDQ4xx;

	      Hxy=G*Hxy+0.5*(nnxxl*(xi-1.0)+nnyyl*(1.0+xi))*duxdy+
		Qxyl*xi*(duydy+duxdx)+0.5*(nnxxl*(1.0+xi)+nnyyl*(xi-1.0))*duydx+
		0.5*duxdz*Qyzl*(1.0+xi)+0.5*duzdx*Qyzl*(xi-1.0)+
		0.5*duydz*Qxzl*(1.0+xi)+0.5*duzdy*Qxzl*(xi-1.0)+2.0*xi*mDQ4xy;
	      Hxz=G*Hxz+(-Qxxl-0.5*Qyyl*(1.0+xi)+xi/3.0)*duxdz+
		Qxzl*xi*(duxdx+duzdz)+duzdx*(Qxxl+0.5*Qyyl*(1.0-xi)+xi/3.0)+
		0.5*Qyzl*((1.0+xi)*duxdy+(xi-1.0)*duydx)+
		0.5*Qxyl*((xi-1.0)*duydz+(1.0+xi)*duzdy)+2.0*xi*mDQ4xz;
	      Hyy=G*Hyy+duxdy*Qxyl*(xi-1.0)+duydy*2.0*nnyyl*xi+
		duydx*Qxyl*(xi+1.0)+duydz*Qyzl*(1.0+xi)+duzdy*Qyzl*(xi-1.0)+
		2.0*xi*mDQ4yy;
	      Hyz=G*Hyz+(-Qyyl-0.5*Qxxl*(1.0+xi)+xi/3.0)*duydz+
		Qyzl*xi*(duydy+duzdz)+duzdy*(Qyyl+0.5*Qxxl*(1.0-xi)+xi/3.0)+
		0.5*Qxzl*((1.0+xi)*duydx+(xi-1.0)*duxdy)+
		0.5*Qxyl*((xi-1.0)*duxdz+(1.0+xi)*duzdx)+2.0*xi*mDQ4yz;
	      
	      /* add advection now*/      

	      vx = ux[i][j][k];
	      vy = uy[i][j][k];
	      vz = uz[i][j][k];

	      /* upwind scheme */

	      vxdQxxdx = max(vx,0.0)*(2.0*Qxx[iup][j][k]  + 3.0*Qxx[i][j][k]   - 6.0*Qxx[idwn][j][k] +     Qxx[idwn2][j][k])/6.0 +
		min(vx,0.0)*(   -Qxx[iup2][j][k] + 6.0*Qxx[iup][j][k] - 3.0*Qxx[i][j][k]    - 2.0*Qxx[idwn][j][k]) /6.0;
	      vydQxxdy = max(vy,0.0)*(2.0*Qxx[i][jup][k]  + 3.0*Qxx[i][j][k]   - 6.0*Qxx[i][jdwn][k] +     Qxx[i][jdwn2][k])/6.0 +
		min(vy,0.0)*(   -Qxx[i][jup2][k] + 6.0*Qxx[i][jup][k] - 3.0*Qxx[i][j][k]    - 2.0*Qxx[i][jdwn][k]) /6.0;
	      vzdQxxdz = max(vz,0.0)*(2.0*Qxx[i][j][kup]  + 3.0*Qxx[i][j][k]   - 6.0*Qxx[i][j][kdwn] +     Qxx[i][j][kdwn2])/6.0 +
		min(vz,0.0)*(   -Qxx[i][j][kup2] + 6.0*Qxx[i][j][kup] - 3.0*Qxx[i][j][k]    - 2.0*Qxx[i][j][kdwn]) /6.0;

	      vxdQxydx = max(vx,0.0)*(2.0*Qxy[iup][j][k]  + 3.0*Qxy[i][j][k]   - 6.0*Qxy[idwn][j][k] +     Qxy[idwn2][j][k])/6.0 +
		min(vx,0.0)*(   -Qxy[iup2][j][k] + 6.0*Qxy[iup][j][k] - 3.0*Qxy[i][j][k]    - 2.0*Qxy[idwn][j][k]) /6.0;
	      vydQxydy = max(vy,0.0)*(2.0*Qxy[i][jup][k]  + 3.0*Qxy[i][j][k]   - 6.0*Qxy[i][jdwn][k] +     Qxy[i][jdwn2][k])/6.0 +
		min(vy,0.0)*(   -Qxy[i][jup2][k] + 6.0*Qxy[i][jup][k] - 3.0*Qxy[i][j][k]    - 2.0*Qxy[i][jdwn][k]) /6.0;
	      vzdQxydz = max(vz,0.0)*(2.0*Qxy[i][j][kup]  + 3.0*Qxy[i][j][k]   - 6.0*Qxy[i][j][kdwn] +     Qxy[i][j][kdwn2])/6.0 +
		min(vz,0.0)*(   -Qxy[i][j][kup2] + 6.0*Qxy[i][j][kup] - 3.0*Qxy[i][j][k]    - 2.0*Qxy[i][j][kdwn]) /6.0;

	      vxdQxzdx = max(vx,0.0)*(2.0*Qxz[iup][j][k]  + 3.0*Qxz[i][j][k]   - 6.0*Qxz[idwn][j][k] +     Qxz[idwn2][j][k])/6.0 +
		min(vx,0.0)*(   -Qxz[iup2][j][k] + 6.0*Qxz[iup][j][k] - 3.0*Qxz[i][j][k]    - 2.0*Qxz[idwn][j][k]) /6.0;
	      vydQxzdy = max(vy,0.0)*(2.0*Qxz[i][jup][k]  + 3.0*Qxz[i][j][k]   - 6.0*Qxz[i][jdwn][k] +     Qxz[i][jdwn2][k])/6.0 +
		min(vy,0.0)*(   -Qxz[i][jup2][k] + 6.0*Qxz[i][jup][k] - 3.0*Qxz[i][j][k]    - 2.0*Qxz[i][jdwn][k]) /6.0;
	      vzdQxzdz = max(vz,0.0)*(2.0*Qxz[i][j][kup]  + 3.0*Qxz[i][j][k]   - 6.0*Qxz[i][j][kdwn] +     Qxz[i][j][kdwn2])/6.0 +
		min(vz,0.0)*(   -Qxz[i][j][kup2] + 6.0*Qxz[i][j][kup] - 3.0*Qxz[i][j][k]    - 2.0*Qxz[i][j][kdwn]) /6.0;
	      
	      vxdQyydx = max(vx,0.0)*(2.0*Qyy[iup][j][k]  + 3.0*Qyy[i][j][k]   - 6.0*Qyy[idwn][j][k] +     Qyy[idwn2][j][k])/6.0 +
		min(vx,0.0)*(   -Qyy[iup2][j][k] + 6.0*Qyy[iup][j][k] - 3.0*Qyy[i][j][k]    - 2.0*Qyy[idwn][j][k]) /6.0;
	      vydQyydy = max(vy,0.0)*(2.0*Qyy[i][jup][k]  + 3.0*Qyy[i][j][k]   - 6.0*Qyy[i][jdwn][k] +     Qyy[i][jdwn2][k])/6.0 +
		min(vy,0.0)*(   -Qyy[i][jup2][k] + 6.0*Qyy[i][jup][k] - 3.0*Qyy[i][j][k]    - 2.0*Qyy[i][jdwn][k]) /6.0;
	      vzdQyydz = max(vz,0.0)*(2.0*Qyy[i][j][kup]  + 3.0*Qyy[i][j][k]   - 6.0*Qyy[i][j][kdwn] +     Qyy[i][j][kdwn2])/6.0 +
		min(vz,0.0)*(   -Qyy[i][j][kup2] + 6.0*Qyy[i][j][kup] - 3.0*Qyy[i][j][k]    - 2.0*Qyy[i][j][kdwn]) /6.0;

	      vxdQyzdx = max(vx,0.0)*(2.0*Qyz[iup][j][k]  + 3.0*Qyz[i][j][k]   - 6.0*Qyz[idwn][j][k] +     Qyz[idwn2][j][k])/6.0 +
		min(vx,0.0)*(   -Qyz[iup2][j][k] + 6.0*Qyz[iup][j][k] - 3.0*Qyz[i][j][k]    - 2.0*Qyz[idwn][j][k]) /6.0;
	      vydQyzdy = max(vy,0.0)*(2.0*Qyz[i][jup][k]  + 3.0*Qyz[i][j][k]   - 6.0*Qyz[i][jdwn][k] +     Qyz[i][jdwn2][k])/6.0 +
		min(vy,0.0)*(   -Qyz[i][jup2][k] + 6.0*Qyz[i][jup][k] - 3.0*Qyz[i][j][k]    - 2.0*Qyz[i][jdwn][k]) /6.0;
	      vzdQyzdz = max(vz,0.0)*(2.0*Qyz[i][j][kup]  + 3.0*Qyz[i][j][k]   - 6.0*Qyz[i][j][kdwn] +     Qyz[i][j][kdwn2])/6.0 +
		min(vz,0.0)*(   -Qyz[i][j][kup2] + 6.0*Qyz[i][j][kup] - 3.0*Qyz[i][j][k]    - 2.0*Qyz[i][j][kdwn]) /6.0;


	      Hxx-=vxdQxxdx+vydQxxdy+vzdQxxdz;
	      Hxy-=vxdQxydx+vydQxydy+vzdQxydz;
	      Hxz-=vxdQxzdx+vydQxzdy+vzdQxzdz;
	      Hyy-=vxdQyydx+vydQyydy+vzdQyydz;
	      Hyz-=vxdQyzdx+vydQyzdy+vzdQyzdz;

	      DEHxx[i][j][k]=Hxx;
	      DEHxy[i][j][k]=Hxy;
	      DEHxz[i][j][k]=Hxz;
	      DEHyy[i][j][k]=Hyy;
	      DEHyz[i][j][k]=Hyz;


			}
		}
	}
	energy = energy/(Lx*Ly*Lz);
}

/* predictor-corrector of time evolution of the director field */
void updateP0(double Pxpr[Lx][Ly][Lz],  double Pypr[Lx][Ly][Lz],  double Pzpr[Lx][Ly][Lz],
	      double Pxold[Lx][Ly][Lz], double Pyold[Lx][Ly][Lz], double Pzold[Lx][Ly][Lz])
{
	int i, j, k;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				Pxpr[i][j][k] = Pxold[i][j][k] + dt*h[i][j][k][0];
				Pypr[i][j][k] = Pyold[i][j][k] + dt*h[i][j][k][1];
				Pzpr[i][j][k] = Pzold[i][j][k] + dt*h[i][j][k][2];

				hold[i][j][k][0] = h[i][j][k][0];
				hold[i][j][k][1] = h[i][j][k][1];
				hold[i][j][k][2] = h[i][j][k][2];
			}
		}
	}
}
void updateP(double Pxnew[Lx][Ly][Lz], double Pynew[Lx][Ly][Lz], double Pznew[Lx][Ly][Lz],
	     double Pxold[Lx][Ly][Lz], double Pyold[Lx][Ly][Lz], double Pzold[Lx][Ly][Lz])
{
	int i, j, k;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				Pxnew[i][j][k] = Pxold[i][j][k] + 0.5*dt*(hold[i][j][k][0] + h[i][j][k][0]);
				Pynew[i][j][k] = Pyold[i][j][k] + 0.5*dt*(hold[i][j][k][1] + h[i][j][k][1]);
				Pznew[i][j][k] = Pzold[i][j][k] + 0.5*dt*(hold[i][j][k][2] + h[i][j][k][2]);
			}
		}
	}
}

/* predictor-corrector of time evolution of Q tensor */
void updateQ0(double Qxxpr[Lx][Ly][Lz],  double Qxypr[Lx][Ly][Lz],  double Qxzpr[Lx][Ly][Lz], double Qyypr[Lx][Ly][Lz],  double Qyzpr[Lx][Ly][Lz],
	      double Qxxold[Lx][Ly][Lz], double Qxyold[Lx][Ly][Lz], double Qxzold[Lx][Ly][Lz], double Qyyold[Lx][Ly][Lz], double Qyzold[Lx][Ly][Lz])
{
	int i, j, k;
	
	double Qsqxx,Qsqxy,Qsqxz,Qsqyy,Qsqyz,Qsqzz;
	double Qxxl,Qxyl,Qxzl,Qyyl,Qyzl,Qzzl;
	double TrQ2,qorder;

	double deltaQ;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				Qxxpr[i][j][k] = Qxxold[i][j][k] + dt*DEHxx[i][j][k];
				Qxypr[i][j][k] = Qxyold[i][j][k] + dt*DEHxy[i][j][k];
				Qxzpr[i][j][k] = Qxzold[i][j][k] + dt*DEHxz[i][j][k];
				Qyypr[i][j][k] = Qyyold[i][j][k] + dt*DEHyy[i][j][k];
				Qyzpr[i][j][k] = Qyzold[i][j][k] + dt*DEHyz[i][j][k];				
 
				DEHxxold[i][j][k]=DEHxx[i][j][k];
				DEHxyold[i][j][k]=DEHxy[i][j][k];
				DEHxzold[i][j][k]=DEHxz[i][j][k];
				DEHyyold[i][j][k]=DEHyy[i][j][k];
				DEHyzold[i][j][k]=DEHyz[i][j][k];
			}
		}
	}
	/* boundary condition for Q */
	for (i = 0; i < Lx; i++) {
	  for (j = 0; j < Ly; j++) {
	    for (k = 0; k < Lz; k++) {
#if BC
	      
	      if(k == 0) {

		Qxxpr[i][j][k]=Qxxpr[i][j][k+1];
		Qxypr[i][j][k]=Qxypr[i][j][k+1];
		Qxzpr[i][j][k]=Qxzpr[i][j][k+1];
		Qyypr[i][j][k]=Qyypr[i][j][k+1];	
		Qyzpr[i][j][k]=Qyzpr[i][j][k+1];
		
	      }
	      
	      if(k == Lz-1) {
		
		Qxxpr[i][j][k]=Qxxpr[i][j][k-1];
		Qxypr[i][j][k]=Qxypr[i][j][k-1];
		Qxzpr[i][j][k]=Qxzpr[i][j][k-1];
		Qyypr[i][j][k]=Qyypr[i][j][k-1];	
		Qyzpr[i][j][k]=Qyzpr[i][j][k-1];
	      }
	      
#endif

#if FIXEDBC

	      gammac=gamma+delta*phi[i][j][k];
	      wallamp=0.0;

	      if (gammac > 2.7){
		wallamp=(0.25+0.75*sqrt((1.0-8.0/(3.0*gammac))));
	      }

	      if(k==0){
		Qxxpr[i][j][k]= wallamp*(sin(angzbot/180.0*Pi)*sin(angzbot/180.0*Pi)*cos(angxybot/180.0*Pi)*cos(angxybot/180.0*Pi)-1.0/3.0);
		Qxypr[i][j][k]= wallamp*sin(angzbot/180.0*Pi)*sin(angzbot/180.0*Pi)*cos(angxybot/180.0*Pi)*sin(angxybot/180.0*Pi);
		Qyypr[i][j][k]= wallamp*(sin(angzbot/180.0*Pi)*sin(angzbot/180.0*Pi)*sin(angxybot/180.0*Pi)*sin(angxybot/180.0*Pi)-1.0/3.0);
		Qxzpr[i][j][k]= wallamp*sin(angzbot/180.0*Pi)*cos(angzbot/180.0*Pi)*cos(angxybot/180.0*Pi);
		Qyzpr[i][j][k]= wallamp*sin(angzbot/180.0*Pi)*cos(angzbot/180.0*Pi)*sin(angxybot/180.0*Pi);
	      }

	      if(k==Lz-1){
		Qxxpr[i][j][k]= wallamp*(sin(angztop/180.0*Pi)*sin(angztop/180.0*Pi)*cos(angxytop/180.0*Pi)*cos(angxytop/180.0*Pi)-1.0/3.0);
		Qxypr[i][j][k]= wallamp*sin(angztop/180.0*Pi)*sin(angztop/180.0*Pi)*cos(angxytop/180.0*Pi)*sin(angxytop/180.0*Pi);
		Qyypr[i][j][k]= wallamp*(sin(angztop/180.0*Pi)*sin(angztop/180.0*Pi)*sin(angxytop/180.0*Pi)*sin(angxytop/180.0*Pi)-1.0/3.0);
		Qxzpr[i][j][k]= wallamp*sin(angztop/180.0*Pi)*cos(angztop/180.0*Pi)*cos(angxytop/180.0*Pi);
		Qyzpr[i][j][k]= wallamp*sin(angztop/180.0*Pi)*cos(angztop/180.0*Pi)*sin(angxytop/180.0*Pi);
	      }

#endif

			}
		}
	}

}
void updateQ(double Qxxnew[Lx][Ly][Lz], double Qxynew[Lx][Ly][Lz], double Qxznew[Lx][Ly][Lz], double Qyynew[Lx][Ly][Lz], double Qyznew[Lx][Ly][Lz],
	     double Qxxold[Lx][Ly][Lz], double Qxyold[Lx][Ly][Lz], double Qxzold[Lx][Ly][Lz], double Qyyold[Lx][Ly][Lz], double Qyzold[Lx][Ly][Lz])
{
	int i, j, k;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				Qxxnew[i][j][k] = Qxxold[i][j][k] + 0.5*dt*(DEHxxold[i][j][k] + DEHxx[i][j][k]);
				Qxynew[i][j][k] = Qxyold[i][j][k] + 0.5*dt*(DEHxyold[i][j][k] + DEHxy[i][j][k]);
				Qxznew[i][j][k] = Qxzold[i][j][k] + 0.5*dt*(DEHxzold[i][j][k] + DEHxz[i][j][k]);
				Qyynew[i][j][k] = Qyyold[i][j][k] + 0.5*dt*(DEHyyold[i][j][k] + DEHyy[i][j][k]);
				Qyznew[i][j][k] = Qyzold[i][j][k] + 0.5*dt*(DEHyzold[i][j][k] + DEHyz[i][j][k]);
			}
		}
	}

	/* boundary condition for Q */
	for (i = 0; i < Lx; i++) {
	  for (j = 0; j < Ly; j++) {
	    for (k = 0; k < Lz; k++) {
#if BC

	      if(k == 0) {
  
		Qxxnew[i][j][k]=Qxxnew[i][j][k+1];
		Qxynew[i][j][k]=Qxynew[i][j][k+1];
		Qxznew[i][j][k]=Qxznew[i][j][k+1];
		Qyynew[i][j][k]=Qyynew[i][j][k+1];	
		Qyznew[i][j][k]=Qyznew[i][j][k+1];
		
	      }
	      if(k == Lz-1) {

		Qxxnew[i][j][k]=Qxxnew[i][j][k-1];
		Qxynew[i][j][k]=Qxynew[i][j][k-1];
		Qxznew[i][j][k]=Qxznew[i][j][k-1];
		Qyynew[i][j][k]=Qyynew[i][j][k-1];	
		Qyznew[i][j][k]=Qyznew[i][j][k-1];
	      }
	      
#endif

#if FIXEDBC

	      gammac=gamma+delta*phi[i][j][k];

	      wallamp=0.0;

	      if (gammac > 2.7){
		wallamp=(0.25+0.75*sqrt((1.0-8.0/(3.0*gammac))));
	      }

	      if(k==0){
		Qxxnew[i][j][k]= wallamp*(sin(angzbot/180.0*Pi)*sin(angzbot/180.0*Pi)*cos(angxybot/180.0*Pi)*cos(angxybot/180.0*Pi)-1.0/3.0);
		Qxynew[i][j][k]= wallamp*sin(angzbot/180.0*Pi)*sin(angzbot/180.0*Pi)*cos(angxybot/180.0*Pi)*sin(angxybot/180.0*Pi);
		Qyynew[i][j][k]= wallamp*(sin(angzbot/180.0*Pi)*sin(angzbot/180.0*Pi)*sin(angxybot/180.0*Pi)*sin(angxybot/180.0*Pi)-1.0/3.0);
		Qxznew[i][j][k]= wallamp*sin(angzbot/180.0*Pi)*cos(angzbot/180.0*Pi)*cos(angxybot/180.0*Pi);
		Qyznew[i][j][k]= wallamp*sin(angzbot/180.0*Pi)*cos(angzbot/180.0*Pi)*sin(angxybot/180.0*Pi);
	      }

	      if(k==Lz-1){
		Qxxnew[i][j][k]= wallamp*(sin(angztop/180.0*Pi)*sin(angztop/180.0*Pi)*cos(angxytop/180.0*Pi)*cos(angxytop/180.0*Pi)-1.0/3.0);
		Qxynew[i][j][k]= wallamp*sin(angztop/180.0*Pi)*sin(angztop/180.0*Pi)*cos(angxytop/180.0*Pi)*sin(angxytop/180.0*Pi);
		Qyynew[i][j][k]= wallamp*(sin(angztop/180.0*Pi)*sin(angztop/180.0*Pi)*sin(angxytop/180.0*Pi)*sin(angxytop/180.0*Pi)-1.0/3.0);
		Qxznew[i][j][k]= wallamp*sin(angztop/180.0*Pi)*cos(angztop/180.0*Pi)*cos(angxytop/180.0*Pi);
		Qyznew[i][j][k]= wallamp*sin(angztop/180.0*Pi)*cos(angztop/180.0*Pi)*sin(angxytop/180.0*Pi);
	      }

#endif
	    }
	  }
	}

}

/*
 *	Plot to a file
 */

void plot(void)
{
	int i, j, k;

	sprintf(filename1, "plot%f.txt", friction);

	output1 = fopen(filename1, "w");

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				fprintf(output1, "%d %d %d %E %E %E %E %E %E %E %E %E\n",
					i, j, k,
					Px[i][j][k], Py[i][j][k], Pz[i][j][k],
					ux[i][j][k], uy[i][j][k], uz[i][j][k],
					mu[i][j][k],
					density[i][j][k],
					phi[i][j][k]);
			}
			fprintf(output1, "\n");
		}
		fprintf(output1, "\n");
	}
	fclose(output1);
}
/* Q tensor output */
void plotQ(void)
{
	int i, j, k;

	int nrots,emax,enxt;
	double m[3][3],d[3],v[3][3];

	sprintf(filename1, "Qtensor.txt");

	output1 = fopen(filename1, "w");

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {


			  m[0][0]=Qxx[i][j][k];
			  m[0][1]=Qxy[i][j][k];
			  m[0][2]=Qxz[i][j][k];
			  m[1][0]=Qxy[i][j][k];
			  m[1][1]=Qyy[i][j][k];
			  m[1][2]=Qyz[i][j][k];
			  m[2][0]=Qxz[i][j][k];
			  m[2][1]=Qyz[i][j][k];
			  m[2][2]= -(m[0][0]+m[1][1]);
			  
			  jacobi(m,d,v,&nrots);

			  if (d[0] > d[1]) {
			    emax=0;
			    enxt=1;
			  }
			  else {
			    emax=1;
			    enxt=0;
			  }
			  if (d[2] > d[emax]) {
			    emax=2;
			  }
			  else if (d[2] > d[enxt]) {
			    enxt=2;
			    }


fprintf(output1, "%d %d %d %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E\n",
				  i, j, k,
				  Qxx[i][j][k],Qxy[i][j][k],Qxz[i][j][k],Qyy[i][j][k],Qyz[i][j][k],
				  ux[i][j][k], uy[i][j][k], uz[i][j][k],
	                          v[0][emax],v[1][emax],v[2][emax],d[emax],
				  mu[i][j][k],
				  density[i][j][k],
	                          phi[i][j][k],flowparameter[i][j][k],compressibility[i][j][k]);
			}
			fprintf(output1, "\n");
		}
		fprintf(output1, "\n");
	}
	fclose(output1);
}


void multipleplot(int n)
{
	int i, j, k;

	sprintf(filename1, "plot.txt.%d", n);

	output1 = fopen(filename1, "w");

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				fprintf(output1, "%d %d %d %E %E %E %E %E %E %E %E %E \n",
					i, j, k,
					Px[i][j][k], Py[i][j][k], Pz[i][j][k],
					ux[i][j][k], uy[i][j][k], uz[i][j][k],
					mu[i][j][k],
					density[i][j][k],
					phi[i][j][k]);
			}
			fprintf(output1, "\n");
		}
		fprintf(output1, "\n");
	}

	fclose(output1);
}
void correction(int n)
{
        double densitychange;
	double momentum[3];
	int i,j,k;

	momentum[0] = 0.0;
	momentum[1] = 0.0;
	momentum[2] = 0.0;

	densitychange=0.0;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
	     		        densitychange = densitychange + (density[i][j][k]-densityinit);
				momentum[0] = momentum[0] + density[i][j][k]*ux[i][j][k];
				momentum[1] = momentum[1] + density[i][j][k]*uy[i][j][k];
				momentum[2] = momentum[2] + density[i][j][k]*uz[i][j][k];

				phi[i][j][k] += phi_average0 - phi_average;
				if (phi[i][j][k] < 0.0) {phi[i][j][k] = 0.0;}
			}
		}
	}
	globalP[0] = momentum[0]/(Lx*Ly*Lz);
	globalP[1] = momentum[1]/(Lx*Ly*Lz);
	globalP[2] = momentum[2]/(Lx*Ly*Lz);
	densitycorrection = densitychange/(Lx*Ly*Lz);

	if(n % stepskip == 0) {
		output1 = fopen("correction.txt", "a");
		fprintf(output1, "%d %E %E %E \n", n, momentum[0], momentum[1], momentum[2]);
		fclose(output1);
	}
}
void add_perturbation(void)
{
	int i,j,k;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				//Px[i][j][k] = Px[i][j][k] + 0.01*(0.5 - drand48())*2.0;
				Py[i][j][k] = Py[i][j][k] + 0.01*(0.5 - drand48())*2.0;
				Pz[i][j][k] = Pz[i][j][k] + 0.01*(0.5 - drand48())*2.0;
			}
		}
	}
} 

/* add perturbation for Q tensor */

void add_Qperturbation(void)
{
	int i,j,k;

	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				
				Qyy[i][j][k] += 0.01*(0.5 - drand48())*2.0;
				Qyz[i][j][k] += 0.01*(0.5 - drand48())*2.0;
			}
		}
	}
} 

/*
 *	mathematical functions
 */
double max(double x, double y)
{
	if (x > y) {
		return x;
	}
	else {
		return y;
	}
}
double min(double x, double y)
{
	if (x > y) {
		return y;
	}
	else {
		return x;
	}
}

/* 3D stencil */
double dx_(double phi[Lx][Ly][Lz], int i, int j, int k)
{
    	double dxphi;
    	int iup = iupa[i];
    	int jup = jupa[j];
    	int kup = kupa[k];

    	int idwn = idwna[i];
    	int jdwn = jdwna[j];
    	int kdwn = kdwna[k];

    	dxphi = (                         phi[iup][j][kup]
    	        + phi[iup][jdwn][k] + 2.0*phi[iup][j][k]    + phi[iup][jup][k]
    	                            +     phi[iup][j][kdwn]                   )/12.0
    	      - (                          phi[idwn][j][kup]
    	        + phi[idwn][jdwn][k] + 2.0*phi[idwn][j][k]    + phi[idwn][jup][k]
    	                             +     phi[idwn][j][kdwn]                    )/12.0;

#if BC
        if (k == 0) {
          dxphi = (phi[iup][jup][k] - phi[idwn][jup][k]
                   + 4.0 * phi[iup][j][k] - 4.0 * phi[idwn][j][k]
                   + phi[iup][jdwn][k] - phi[idwn][jdwn][k]) / 12.0;
        } else if (k == Lz - 1) {
          dxphi = (phi[iup][jup][k] - phi[idwn][jup][k]
                   + 4.0 * phi[iup][j][k] - 4.0 * phi[idwn][j][k]
                   + phi[iup][jdwn][k] - phi[idwn][jdwn][k]) / 12.0;
        }
#endif
    	return dxphi;
}
double dy_(double phi[Lx][Ly][Lz], int i, int j, int k)
{
    	double dyphi;
    	int iup = iupa[i];
    	int jup = jupa[j];
    	int kup = kupa[k];

    	int idwn = idwna[i];
    	int jdwn = jdwna[j];
    	int kdwn = kdwna[k];

    	dyphi = (                         phi[i][jup][kup]
    	        + phi[idwn][jup][k] + 2.0*phi[i][jup][k]    + phi[iup][jup][k]
    	                            +     phi[i][jup][kdwn]                   )/12.0
    	      - (                          phi[i][jdwn][kup]
    	        + phi[idwn][jdwn][k] + 2.0*phi[i][jdwn][k]    + phi[iup][jdwn][k]
    	                             +     phi[i][jdwn][kdwn]                    )/12.0;

#if BC
        if (k == 0) {
          dyphi = (phi[iup][jup][k] - phi[iup][jdwn][k]
                   + 4.0 * phi[i][jup][k] - 4.0 * phi[i][jdwn][k]
                   + phi[idwn][jup][k] - phi[idwn][jdwn][k]) / 12.0;
	} else if (k == Lz - 1) {
          dyphi = (phi[iup][jup][k] - phi[iup][jdwn][k]
                   + 4.0 * phi[i][jup][k] - 4.0 * phi[i][jdwn][k]
                   + phi[idwn][jup][k] - phi[idwn][jdwn][k]) / 12.0;
        }
#endif
    	return dyphi;
}
double dz_(double phi[Lx][Ly][Lz], int i, int j, int k)
{
    	double dzphi;
    	int iup = iupa[i];
    	int jup = jupa[j];
    	int kup = kupa[k];
        int kup2 = kupa[kupa[k]];

    	int idwn = idwna[i];
    	int jdwn = jdwna[j];
    	int kdwn = kdwna[k];
	int kdwn2 = kdwna[kdwna[k]];

    	dzphi = (                         phi[i][jup][kup]
    	        + phi[idwn][j][kup] + 2.0*phi[i][j][kup]    + phi[iup][j][kup]
    	                            +     phi[i][jdwn][kup]                   )/12.0
    	      - (                          phi[i][jup][kdwn]
    	        + phi[idwn][j][kdwn] + 2.0*phi[i][j][kdwn]    + phi[iup][j][kdwn]
    	                             +     phi[i][jdwn][kdwn]                    )/12.0;
#if BC
        if (k == 0) {
          dzphi = (-phi[i][j][kup2] + 4.0 * phi[i][j][kup] - 3.0 * phi[i][j][k]) / 2.0;
        } else if (k == Lz - 1) {
          dzphi = (phi[i][j][kdwn2] - 4.0 * phi[i][j][kdwn] + 3.0 * phi[i][j][k]) / 2.0;
        }
#endif
	return dzphi;
}
double laplacian_(double phi[Lx][Ly][Lz], int i , int j , int k)
{
   	double laplacianphi;
   	int iup = iupa[i];
   	int jup = jupa[j];
   	int kup = kupa[k];
	int kup2 = kupa[kupa[k]];
        int kup3 = kupa[kupa[kupa[k]]];
	
	int idwn = idwna[i];
	int jdwn = jdwna[j];
	int kdwn = kdwna[k];
        int kdwn2 = kdwna[kdwna[k]];
        int kdwn3 = kdwna[kdwna[kdwna[k]]];

    	laplacianphi = (                         phi[iup][j][kup]
        	       + phi[iup][jdwn][k] + 2.0*phi[iup][j][k]    + phi[iup][jup][k]
                       	                   +     phi[iup][j][kdwn]                   )/6.0

                     + (     phi[i][jdwn][kup]  +  2.0*phi[i][j][kup]  +     phi[i][jup][kup]
                       + 2.0*phi[i][jdwn][k]    - 24.0*phi[i][j][k]    + 2.0*phi[i][jup][k]
                       +     phi[i][jdwn][kdwn] +  2.0*phi[i][j][kdwn] +     phi[i][jup][kdwn])/6.0

                     + (                          phi[idwn][j][kup]
                       + phi[idwn][jdwn][k] + 2.0*phi[idwn][j][k]    + phi[idwn][jup][k]
                       	                    +     phi[idwn][j][kdwn]                    )/6.0;

#if BC
        if (k == 0) {
          laplacianphi = (phi[idwn][jup][k] + 4.0 * phi[i][jup][k] + phi[iup][jup][k]
                          + 4.0 * phi[idwn][j][k] - 20.0 * phi[i][j][k] + 4.0 * phi[iup][j][k]
                          + phi[idwn][jdwn][k] + 4.0 * phi[i][jdwn][k] + phi[iup][jdwn][k]) / \
6.0

            - phi[i][j][kup3] + 4.0 * phi[i][j][kup2] - 5.0 * phi[i][j][kup] + 2.0 * phi[i][j]\
	    [k];
        } else if (k == Lz - 1) {
          laplacianphi = (phi[idwn][jup][k] + 4.0 * phi[i][jup][k] + phi[iup][jup][k]
                          + 4.0 * phi[idwn][j][k] - 20.0 * phi[i][j][k] + 4.0 * phi[iup][j][k]
                          + phi[idwn][jdwn][k] + 4.0 * phi[i][jdwn][k] + phi[iup][jdwn][k]) / \
6.0

            - phi[i][j][kdwn3] + 4.0 * phi[i][j][kdwn2] - 5.0 * phi[i][j][kdwn] + 2.0 * phi[i]\
	    [j][k];
        }
#endif
	return laplacianphi;
}

void computedensityfluctuations(int boxsize)
{
    int i,j,k,j1,k1,ny,nz;
    int ndata;
    double av,av2,mass,totav,var;
    
    var=0.0;
    av=0.0;
    av2=0.0;
    ndata=0;
    
    for(i=0 ; i<Lx; i++) {
        for(j=0 ; j<Ly; j++) {
            for(k=0 ; k<Lz; k++) {
                
                mass=0.0;
                for(ny=0 ; ny<boxsize; ny++){
                    for(nz=0; nz<boxsize; nz++){
                        j1=(j+ny)%Ly;
                        k1=(k+nz)%Lz;
                        mass += phi[i][j1][k1];
//                        if(boxsize==30) fprintf(stderr,"%d %d %d %d %d %d %lf\n",j,k,ny,nz,j1,k1,phi[i][j1][k1]);
                    }
                }
		av += mass;
		av2 += mass*mass;

		ndata += 1;
	    }
	}
    }
     
    av2=av2/(float)(ndata);
    av=av/(float)(ndata);
    //fprintf(stderr,"%d %d %lf %lf\n",boxsize,ndata,av2,av);
    var = av2-av*av;

    phiav[boxsize] += av;
    phivar[boxsize] += var;

    //fprintf(stderr,"average var=%lf\n",var);
    
}

void computeaveragevelocity(void)
{
    	

}

void plotdensityfluctuations(void)
{
    int n;
    
    output1 = fopen("densityfluctuations.dat","w");
    ndatatot += 1;
    for (n=5; n<=Ly/2; n++){
        computedensityfluctuations(n);
        fprintf(output1,"%lf %lf %d\n",phiav[n]/(float)(ndatatot),phivar[n]/(float)(ndatatot),n);
        fflush(output1);
    } 
    fclose(output1);
}

void plotbifurcation(void)
{
	output1 = fopen("bifurcation.txt","a");

	fprintf(output1, "%E %E %E %E %d \n", friction, Vy, Vz, umax, R1);

	fclose(output1);
}
void plotV(int n)
{
	sprintf(filename1, "V%f.txt", friction);
	output1 = fopen(filename1, "a");

	fprintf(output1, "%d %E %E \n", n, Vy, Vz);

	fclose(output1);
}
void plotumax(int n)
{
	int i,j,k;

	double umag;
	umax = 0.0;
	for (i = 0; i < Lx; i++) {
		for (j = 0; j < Ly; j++) {
			for (k = 0; k < Lz; k++) {
				umag = sqrt(ux[i][j][k]*ux[i][j][k] + uy[i][j][k]*uy[i][j][k] + uz[i][j][k]*uz[i][j][k]);

				if (umag > umax) {umax = umag;}
			}
		}
	}
	sprintf(filename1, "umax%f.txt", friction);
	output1 = fopen(filename1, "a");

	fprintf(output1, "%d %E \n", n, umax);

	fclose(output1);
}
/*
 *	Stencils for 2DQ9 lattice vectors
 *
				dphidy = ( phi[i][jup][kup]  + 4.0*phi[i][jup][k]  + phi[i][jup][kdwn]
					   - phi[i][jdwn][kup] - 4.0*phi[i][jdwn][k] - phi[i][jdwn][kdwn])/12.0;

				dphidz = ( phi[i][jup][kup]  + 4.0*phi[i][j][kup]  + phi[i][jdwn][kup]
					   - phi[i][jup][kdwn] - 4.0*phi[i][j][kdwn] - phi[i][jdwn][kdwn])/12.0;

				laplacianphi = (     phi[i][jdwn][kup]  +  4.0*phi[i][j][kup]  +     phi[i][jup][kup]
					         + 4.0*phi[i][jdwn][k]    - 20.0*phi[i][j][k]    + 4.0*phi[i][jup][k]
					         +     phi[i][jdwn][kdwn] +  4.0*phi[i][j][kdwn] +     phi[i][jup][kdwn])/6.0;

 */


/* subroutine jacobi to find eigenvalues and eigenvectors from numerical recipes */


#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);
#define n 3
void jacobi(double (*a)[n], double d[], double (*v)[n], int *nrot)
{
  int j,iq,ip,i;
  double tresh,theta,tau,t,sm,s,h,g,c;
  double b[n],z[n];

  for (ip=0;ip<n;ip++) {
    for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  for (ip=0;ip<n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=0;ip< n-1;ip++) {
      for (iq=ip+1;iq<n;iq++)
	sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++) {
	g=100.0*fabs(a[ip][iq]);
	if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip])
	    && (fabs(d[iq])+g) == fabs(d[iq]))
	  a[ip][iq]=0.0;
	else if (fabs(a[ip][iq]) > tresh) {
	  h=d[iq]-d[ip];
	  if ((fabs(h)+g) == fabs(h))
	    t=(a[ip][iq])/h;
	  else {
	    theta=0.5*h/(a[ip][iq]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq]=0.0;
	  for (j=0;j<=ip-1;j++) {
	    ROTATE(a,j,ip,j,iq)
	      }
	  for (j=ip+1;j<=iq-1;j++) {
	    ROTATE(a,ip,j,j,iq)
	      }
	  for (j=iq+1;j<n;j++) {
	    ROTATE(a,ip,j,iq,j)
	      }
	  for (j=0;j<n;j++) {
	    ROTATE(v,j,ip,j,iq)
	      }
	  ++(*nrot);
	}
      }
    }
    for (ip=0;ip<n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  printf("Too many iterations in routine jacobi");
  exit(0);
}
#undef n
#undef ROTATE

