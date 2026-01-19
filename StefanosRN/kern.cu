#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <cuda.h>
#include <qd/dd_real.h>
#include <qd/fpu.h>
#include <gqd.cu>
#include "parm.h"

__device__ gdd_real tqrep[4][M];
__device__ gdd_real tqimp[4][M];
__device__ gdd_real tprep[4][M];
__device__ gdd_real tpimp[4][M];
__device__ gdd_real tqrem[4][M];
__device__ gdd_real tqimm[4][M];
__device__ gdd_real tprem[4][M];
__device__ gdd_real tpimm[4][M];

__device__ double fix;
__device__ gdd_real mass;
__device__ gdd_real dtheta;
__device__ gdd_real dx;
__device__ gdd_real dt;
__device__ gdd_real aa;
__device__ gdd_real ss;
__device__ gdd_real mm;
__device__ int my_rank;
__device__ int p;

__device__ gdd_real Power (gdd_real xx, int yy){
gdd_real prod = make_dd(1.0);
 if ( yy < 0) { 
 for (int ii = 0; ii < -yy; ii++)
 prod = prod/xx;
 } else {
 for (int ii = 0; ii < yy; ii++)
 prod = prod*xx;
 }
return  prod;
}

__device__ gdd_real wenoL (gdd_real f, gdd_real fm, gdd_real fp, gdd_real fm2, gdd_real fp2){

       gdd_real beta0, beta1, beta2, w0, w1, w2, ww0, ww1, ww2, eps, res;

         ww0 = make_dd(0.1);
         ww1 = make_dd(0.6);
         ww2 = make_dd(0.3);

         res = (ww0*(2.0*fm2-7.0*fm+11.0*f) + 
               ww1*(-1.0*fm+5.0*f+2.0*fp) + 
               ww2*(2.0*f+5.0*fp-1.0*fp2)) / 6.0;

         return (res);

}

__device__ gdd_real wenoR (gdd_real f, gdd_real fm, gdd_real fp, gdd_real fm2, gdd_real fp2){

       gdd_real beta0, beta1, beta2, w0, w1, w2, ww0, ww1, ww2, eps, res;

         ww0 = make_dd(0.1);         
         ww1 = make_dd(0.6);
         ww2 = make_dd(0.3);

         res = (ww0*(2.0*fp2-7.0*fp+11.0*f) + 
               ww1*(-1.0*fp+5.0*f+2.0*fm) + 
               ww2*(2.0*f+5.0*fm-1.0*fm2)) / 6.0;

         return (res);

}

/* -------------------------------------------------- */

__global__ void
kernel_init (gdd_real *qre, gdd_real *qim, gdd_real *pre, gdd_real *pim, gdd_real dtheta_in, gdd_real mass_in, gdd_real aa_in, gdd_real ss_in, gdd_real mm_in, gdd_real dx_in, gdd_real dt_in, gdd_real *theta, gdd_real *r_c, gdd_real *bb_c, gdd_real *cc_c, gdd_real *ee_c, gdd_real* a_re_c, gdd_real* a_im_c, gdd_real* d_re_c, gdd_real* d_im_c, gdd_real* v_re_c, gdd_real* v_im_c, int my_rank_in, int p_in)
{

  int k;
  int j;

  gdd_real rpt, ypt, stheta, charge, S, em, Cos, Sin, Csc, Cot;
  gdd_real x, bir,  Ataurho, Arhorho, Ayy, BtauReal,  ff, df; 
  gdd_real BtauImag, BrhoReal, By, BrhoImag, CReal, CImag;  

  k = blockIdx.y * blockDim.y + threadIdx.y;
  j = blockIdx.x * blockDim.x + threadIdx.x;

  mass = mass_in;
  aa = aa_in;
  ss = ss_in;
  mm = mm_in;
  dtheta = dtheta_in;
  dx = dx_in;
  dt = dt_in;
  S = make_dd(Xmax);
  my_rank = my_rank_in;
  p = p_in;

  charge = aa;
  em = mm;

  fix = 4.0;
  if ( (int)(to_double(mm))%2==0 ) { fix = 2.0; }

  qre[idx (k, j)] = make_dd(0.0);
  qim[idx (k, j)] = make_dd(0.0);
  pre[idx (k, j)] = make_dd(0.0);
  pim[idx (k, j)] = make_dd(0.0);

  //if (k >= 1 && k < M - 1)
      {

         //Cos = 1.0 + make_dd(k - 0.0) * dtheta;
         Cos = (theta[k]);
         Sin = sqrt(1.0 - Cos*Cos);
         Csc = 1.0 / Sin;
         Cot = Cos * Csc; 
 
         rpt = Xmin + ((Xmax-Xmin)/p)*my_rank + make_dd(j - 0.0) *dx;
         x = rpt / S;

         ff = 1.0 - 2.0*mass*(1.0/rpt - 1.0/S) + charge*charge*Power(1.0/rpt - 1.0/S, 2);
         df = 2.0*mass/(rpt*rpt) - 2.0*charge*charge*(1.0/rpt - 1.0/S)/(rpt*rpt);

	 bir = x*( 4.0*mass/S*(1.0-x) - x*x + 1.0 + 2.0*x )*( 2.0 + ff*(-2.0*x*(1.0-2.0*mass/S) - (4.0*mass/S+1.0) + x*x) );
     
         Ataurho = 2.0*x*Power(1.0-x,2)*( 1.0 + ff*( -2.0*x*(1.0-2.0*mass/S)-(4.0*mass/S+1.0) + x*x ) ) / bir;
         
         Arhorho = x*Power(1.0-x,4)*ff / bir;
         
         Ayy = make_dd(0.0);
         
         BtauReal = 0.0-x*Power(1.0-x,2)*( df*(4.0*mass/S*(1.0-x) - x*x + 1.0 + 2.0*x) - 2.0*ff*(2.0*mass/S + x - 1.0)/S ) / bir;
         
         BrhoReal = 0.0-Power(1.0-x,3)*( x*(x-1.0)*df + 2.0*x*ff/S ) / bir;
         
         By = make_dd(0.0);
         CReal = 0.0-Power(1.0-x,2)*( x*(1.0-x)*df + ss*(ss+1.0)/S ) / bir / (x*S);
         BtauImag = make_dd(0.0);
         BrhoImag = make_dd(0.0);
         CImag = make_dd(0.0);

         bb_c[idx (k, j)] = 0.0-(Ataurho - sqrt(Ataurho*Ataurho+4.0*Arhorho))/2.0;
         cc_c[idx (k, j)] = Ayy;
         ee_c[idx (k, j)] = Ataurho + bb_c[idx (k, j)];

         a_re_c[idx (k, j)] = 0.0-BtauReal;
         a_im_c[idx (k, j)] = 0.0-BtauImag;

         d_re_c[idx (k, j)] = BrhoReal -bb_c[idx (k, j)]*BtauReal; 
         d_im_c[idx (k, j)] = BrhoImag -bb_c[idx (k, j)]*BtauImag;

         v_re_c[idx (k, j)] = 0.0-CReal;
         v_im_c[idx (k, j)] = 0.0-CImag;

        }

         //Cos = 1.0 + make_dd(k - 0.0) * dtheta;
         Cos = (theta[k]);
         Sin = sqrt(1.0 - Cos*Cos);
         x = Xmin + ((Xmax-Xmin)/p)*my_rank + make_dd(j - 0.0) *dx;

         qre[idx (k, j)] = make_dd(0.0);
         qim[idx (k, j)] = make_dd(0.0);
         pim[idx (k, j)] = make_dd(0.0);

	 //this is L=5
	 //pre[idx (k, j)] = exp(0.0-(x-1.0)*(x-1.0)/0.1)*(63.0*Power(Cos,4)-70.0*Cos*Cos+15.0)*Cos;

	 //this is L=4
	 //pre[idx (k, j)] = exp(0.0-(x-1.15)*(x-1.15)/0.1)*(35.0*Power(Cos,4)-30.0*Cos*Cos+3.0);

	 //this is L=3
	 //pre[idx (k, j)] = exp(0.0-(x-1.0)*(x-1.0)/0.1)*(5.0*Power(Cos,2)-3.0)*Cos;

	 //this is L=2
	 //pre[idx (k, j)] = exp(0.0-(x-1.15)*(x-1.15)/0.1)*(3.0*Cos*Cos-1.0);

	 //this is L=1
	 //pre[idx (k, j)] = exp(0.0-(x-1.0)*(x-1.0)/0.1)*Cos;

	 //this is L=0
	 qre[idx (k, j)] = exp(0.0-(x-1.0)*(x-1.0)/0.1);
         pre[idx (k, j)] = bb_c[idx (k, j)] * qre[idx (k, j)] * (-2.0*(x-1.0)/0.1);
         if (x > 8.0) { pre[idx (k, j)] = make_dd(0.0); qre[idx (k, j)] = make_dd(0.0);}

}


__global__ void
kernel_rhs (gdd_real *rhs_qre, gdd_real *rhs_qim, gdd_real *rhs_pre, gdd_real *rhs_pim, gdd_real *qre, gdd_real *qim, gdd_real *pre, gdd_real *pim, gdd_real *theta, gdd_real *r_c, double timer, gdd_real *bb_c, gdd_real *cc_c, gdd_real *ee_c, gdd_real* a_re_c, gdd_real* a_im_c, gdd_real* d_re_c, gdd_real* d_im_c, gdd_real* v_re_c, gdd_real* v_im_c, gdd_real *d_buff, gdd_real *d_buff2)
{
  int b;
  int a;

  b = blockIdx.y * blockDim.y + threadIdx.y;
  a = blockIdx.x * blockDim.x + threadIdx.x;

  gdd_real ctan, stheta, ctheta, dbdx, dre;
  gdd_real qrey = make_dd(0.0);
  gdd_real qreyy = make_dd(0.0);
  gdd_real qimy = make_dd(0.0);
  gdd_real qimyy = make_dd(0.0);
  gdd_real ll_re, ll_im, f_mid_p, f_mid_m;
  gdd_real qrex, qimx, prex, pimx;

        gdd_real dc1[12];

        dc1[0] = make_dd(-1.0) / 1260.0;
        dc1[1] = make_dd(5.0) / 504.0;
        dc1[2] = make_dd(-5.0) / 84.0;
        dc1[3] = make_dd(5.0) / 21.0;
        dc1[4] = make_dd(-5.0) / 6.0;
        dc1[5] = make_dd(0.0);
        dc1[6] = 0.0-dc1[4];
        dc1[7] = 0.0-dc1[3];
        dc1[8] = 0.0-dc1[2];
        dc1[9] = 0.0-dc1[1];
        dc1[10] = 0.0-dc1[0];

        gdd_real dc2[12];

        dc2[0] = make_dd(1.0) / 6300.0;
        dc2[1] = make_dd(-5.0) / 2016.0; 
        dc2[2] = make_dd(5.0) / 252.0;
        dc2[3] = make_dd(-5.0) / 42.0;
        dc2[4] = make_dd(5.0) / 6.0;
        dc2[5] = make_dd(-5269.0) / 3600.0;
        dc2[6] = dc2[4];
        dc2[7] = dc2[3];
        dc2[8] = dc2[2];
        dc2[9] = dc2[1];
        dc2[10] = dc2[0];


  if (my_rank < p-1) {

  tqrep[1][b] = d_buff[b];
  tqrep[2][b] = d_buff[M+b];
  tqrep[3][b] = d_buff[2*M+b];

  tqimp[1][b] = d_buff[3*M+b];
  tqimp[2][b] = d_buff[3*M+M+b];
  tqimp[3][b] = d_buff[3*M+2*M+b];

  tprep[1][b] = d_buff[6*M+b];
  tprep[2][b] = d_buff[6*M+M+b];
  tprep[3][b] = d_buff[6*M+2*M+b];

  tpimp[1][b] = d_buff[9*M+b];
  tpimp[2][b] = d_buff[9*M+M+b];
  tpimp[3][b] = d_buff[9*M+2*M+b];

  }

  if (my_rank > 0) {

  tqrem[1][b] = d_buff2[b];
  tqrem[2][b] = d_buff2[M+b];
  tqrem[3][b] = d_buff2[2*M+b];

  tqimm[1][b] = d_buff2[3*M+b];
  tqimm[2][b] = d_buff2[3*M+M+b];
  tqimm[3][b] = d_buff2[3*M+2*M+b];

  tprem[1][b] = d_buff2[6*M+b];
  tprem[2][b] = d_buff2[6*M+M+b];
  tprem[3][b] = d_buff2[6*M+2*M+b];

  tpimm[1][b] = d_buff2[9*M+b];
  tpimm[2][b] = d_buff2[9*M+M+b];
  tpimm[3][b] = d_buff2[9*M+2*M+b];

  }


         dbdx = make_dd(0.0);
         if (a >= 3 && a < N-3){
         dbdx = (45.0*bb_c[idx (b, a+1)]-45.0*bb_c[idx (b, a-1)]
                  -9.0*bb_c[idx (b, a+2)]+9.0*bb_c[idx (b, a-2)] 
                  +bb_c[idx (b, a+3)]-bb_c[idx (b, a-3)])/(60.0*dx);
         }

         if (a == 2) {
	 dbdx = (80.0*bb_c[idx (b, a+1)]-24.0*bb_c[idx (b, a-1)]
                       -30.0*bb_c[idx (b, a+2)]+2.0*bb_c[idx (b, a-2)]
                       +8.0*bb_c[idx (b, a+3)]-bb_c[idx (b, a+4)]
                       -35.0*bb_c[idx (b, a)])/(60.0*dx);
         }


         if (a == N-3) {
	 dbdx = (24.0*bb_c[idx (b, a+1)]-80.0*bb_c[idx (b, a-1)]
                       -2.0*bb_c[idx (b, a+2)]+30.0*bb_c[idx (b, a-2)]
                       -8.0*bb_c[idx (b, a-3)]+bb_c[idx (b, a-4)]
                       +35.0*bb_c[idx (b, a)])/(60.0*dx);
         }

         if (a == 1) {
	 dbdx = (150.0*bb_c[idx (b, a+1)]-10.0*bb_c[idx (b, a-1)]
                       -100.0*bb_c[idx (b, a+2)]+2.0*bb_c[idx (b, a+5)]
                       +50.0*bb_c[idx (b, a+3)]-15.0*bb_c[idx (b, a+4)]
                       -77.0*bb_c[idx (b, a)])/(60.0*dx);
         }

        if (a == N-2) {
	dbdx = (10.0*bb_c[idx (b, a+1)]-150.0*bb_c[idx (b, a-1)]
                       -2.0*bb_c[idx (b, a-5)]+100.0*bb_c[idx (b, a-2)]
                       -50.0*bb_c[idx (b, a-3)]+15.0*bb_c[idx (b, a-4)]
                       +77.0*bb_c[idx (b, a)])/(60.0*dx);
        }


        if (a == 0) {
	dbdx = (360.0*bb_c[idx (b, a+1)]-10.0*bb_c[idx (b, a+6)]
                       -450.0*bb_c[idx (b, a+2)]+72.0*bb_c[idx (b, a+5)]
                       +400.0*bb_c[idx (b, a+3)]-225.0*bb_c[idx (b, a+4)]
                       -147.0*bb_c[idx (b, a)])/(60.0*dx);

        }


        if (a == N-1) {
	dbdx = (10.0*bb_c[idx (b, a-6)]-360.0*bb_c[idx (b, a-1)]
                       -72.0*bb_c[idx (b, a-5)]+450.0*bb_c[idx (b, a-2)]
                       -400.0*bb_c[idx (b, a-3)]+225.0*bb_c[idx (b, a-4)]
                       +147.0*bb_c[idx (b, a)])/(60.0*dx);
        }

        dre = d_re_c[idx (b, a)] - dbdx*ee_c[idx (b, a)];
 

  //if (b >= 5 && b < M - 5)
    if (a >= 0 && a < N - 0)
      {

        ctheta = (theta[b]);

        qrey = make_dd(0.0);
        qreyy = make_dd(0.0);
        qimy = make_dd(0.0);
        qimyy = make_dd(0.0);

       if (a == 2) {        
        if (my_rank == 0){
	qrex = (80.0*qre[idx (b, a+1)]-24.0*qre[idx (b, a-1)]
                       -30.0*qre[idx (b, a+2)]+2.0*qre[idx (b, a-2)]
                       +8.0*qre[idx (b, a+3)]-qre[idx (b, a+4)]
                       -35.0*qre[idx (b, a)])/(60.0*dx);
	qimx = (80.0*qim[idx (b, a+1)]-24.0*qim[idx (b, a-1)]
                       -30.0*qim[idx (b, a+2)]+2.0*qim[idx (b, a-2)]
                       +8.0*qim[idx (b, a+3)]-qim[idx (b, a+4)]
                       -35.0*qim[idx (b, a)])/(60.0*dx);
	prex = (80.0*pre[idx (b, a+1)]-24.0*pre[idx (b, a-1)]
                       -30.0*pre[idx (b, a+2)]+2.0*pre[idx (b, a-2)]
                       +8.0*pre[idx (b, a+3)]-pre[idx (b, a+4)]
                       -35.0*pre[idx (b, a)])/(60.0*dx);
	pimx = (80.0*pim[idx (b, a+1)]-24.0*pim[idx (b, a-1)]
                       -30.0*pim[idx (b, a+2)]+2.0*pim[idx (b, a-2)]
                       +8.0*pim[idx (b, a+3)]-pim[idx (b, a+4)]
                       -35.0*pim[idx (b, a)])/(60.0*dx);
        } else {
	qrex = (45.0*qre[idx (b, a+1)]-45.0*qre[idx (b, a-1)]
                       -9.0*qre[idx (b, a+2)]+9.0*qre[idx (b, a-2)]
                       +qre[idx (b, a+3)]-tqrem[1][b])/(60.0*dx);
	qimx = (45.0*qim[idx (b, a+1)]-45.0*qim[idx (b, a-1)]
                       -9.0*qim[idx (b, a+2)]+9.0*qim[idx (b, a-2)]
                       +qim[idx (b, a+3)]-tqimm[1][b])/(60.0*dx);
	prex = (45.0*pre[idx (b, a+1)]-45.0*pre[idx (b, a-1)]
                       -9.0*pre[idx (b, a+2)]+9.0*pre[idx (b, a-2)]
                       +pre[idx (b, a+3)]-tprem[1][b])/(60.0*dx);
	pimx = (45.0*pim[idx (b, a+1)]-45.0*pim[idx (b, a-1)]
                       -9.0*pim[idx (b, a+2)]+9.0*pim[idx (b, a-2)]
                       +pim[idx (b, a+3)]-tpimm[1][b])/(60.0*dx);
        }
        }


        if (a == N-3) {
        if (my_rank == p-1){
	qrex = (24.0*qre[idx (b, a+1)]-80.0*qre[idx (b, a-1)]
                       -2.0*qre[idx (b, a+2)]+30.0*qre[idx (b, a-2)]
                       -8.0*qre[idx (b, a-3)]+qre[idx (b, a-4)]
                       +35.0*qre[idx (b, a)])/(60.0*dx);
	qimx = (24.0*qim[idx (b, a+1)]-80.0*qim[idx (b, a-1)]
                       -2.0*qim[idx (b, a+2)]+30.0*qim[idx (b, a-2)]
                       -8.0*qim[idx (b, a-3)]+qim[idx (b, a-4)]
                       +35.0*qim[idx (b, a)])/(60.0*dx);
	prex = (24.0*pre[idx (b, a+1)]-80.0*pre[idx (b, a-1)]
                       -2.0*pre[idx (b, a+2)]+30.0*pre[idx (b, a-2)]
                       -8.0*pre[idx (b, a-3)]+pre[idx (b, a-4)]
                       +35.0*pre[idx (b, a)])/(60.0*dx);
	pimx = (24.0*pim[idx (b, a+1)]-80.0*pim[idx (b, a-1)]
                       -2.0*pim[idx (b, a+2)]+30.0*pim[idx (b, a-2)]
                       -8.0*pim[idx (b, a-3)]+pim[idx (b, a-4)]
                       +35.0*pim[idx (b, a)])/(60.0*dx);
        } else {
	qrex = (45.0*qre[idx (b, a+1)]-45.0*qre[idx (b, a-1)]
                       -9.0*qre[idx (b, a+2)]+9.0*qre[idx (b, a-2)]
                       +tqrep[1][b]-qre[idx (b, a-3)])/(60.0*dx);
	qimx = (45.0*qim[idx (b, a+1)]-45.0*qim[idx (b, a-1)]
                       -9.0*qim[idx (b, a+2)]+9.0*qim[idx (b, a-2)]
                       +tqimp[1][b]-qim[idx (b, a-3)])/(60.0*dx);
	prex = (45.0*pre[idx (b, a+1)]-45.0*pre[idx (b, a-1)]
                       -9.0*pre[idx (b, a+2)]+9.0*pre[idx (b, a-2)]
                       +tprep[1][b]-pre[idx (b, a-3)])/(60.0*dx);
	pimx = (45.0*pim[idx (b, a+1)]-45.0*pim[idx (b, a-1)]
                       -9.0*pim[idx (b, a+2)]+9.0*pim[idx (b, a-2)]
                       +tpimp[1][b]-pim[idx (b, a-3)])/(60.0*dx);
        }
        }

        if (a == 1) {
        if (my_rank == 0){
	qrex = (150.0*qre[idx (b, a+1)]-10.0*qre[idx (b, a-1)]
                       -100.0*qre[idx (b, a+2)]+2.0*qre[idx (b, a+5)]
                       +50.0*qre[idx (b, a+3)]-15.0*qre[idx (b, a+4)]
                       -77.0*qre[idx (b, a)])/(60.0*dx);
	qimx = (150.0*qim[idx (b, a+1)]-10.0*qim[idx (b, a-1)]
                       -100.0*qim[idx (b, a+2)]+2.0*qim[idx (b, a+5)]
                       +50.0*qim[idx (b, a+3)]-15.0*qim[idx (b, a+4)]
                       -77.0*qim[idx (b, a)])/(60.0*dx);
	prex = (150.0*pre[idx (b, a+1)]-10.0*pre[idx (b, a-1)]
                       -100.0*pre[idx (b, a+2)]+2.0*pre[idx (b, a+5)]
                       +50.0*pre[idx (b, a+3)]-15.0*pre[idx (b, a+4)]
                       -77.0*pre[idx (b, a)])/(60.0*dx);
	pimx = (150.0*pim[idx (b, a+1)]-10.0*pim[idx (b, a-1)]
                       -100.0*pim[idx (b, a+2)]+2.0*pim[idx (b, a+5)]
                       +50.0*pim[idx (b, a+3)]-15.0*pim[idx (b, a+4)]
                       -77.0*pim[idx (b, a)])/(60.0*dx);
        } else {
        qrex = (45.0*qre[idx (b, a+1)]-45.0*qre[idx (b, a-1)]
                       -9.0*qre[idx (b, a+2)]+9.0*tqrem[1][b]
                       +qre[idx (b, a+3)]-tqrem[2][b])/(60.0*dx);
        qimx = (45.0*qim[idx (b, a+1)]-45.0*qim[idx (b, a-1)]
                       -9.0*qim[idx (b, a+2)]+9.0*tqimm[1][b]
                       +qim[idx (b, a+3)]-tqimm[2][b])/(60.0*dx);
        prex = (45.0*pre[idx (b, a+1)]-45.0*pre[idx (b, a-1)]
                       -9.0*pre[idx (b, a+2)]+9.0*tprem[1][b]
                       +pre[idx (b, a+3)]-tprem[2][b])/(60.0*dx);
        pimx = (45.0*pim[idx (b, a+1)]-45.0*pim[idx (b, a-1)]
                       -9.0*pim[idx (b, a+2)]+9.0*tpimm[1][b]
                       +pim[idx (b, a+3)]-tpimm[2][b])/(60.0*dx);
        }
        }

        if (a == N-2) {
        if (my_rank == p-1){
	qrex = (10.0*qre[idx (b, a+1)]-150.0*qre[idx (b, a-1)]
                       -2.0*qre[idx (b, a-5)]+100.0*qre[idx (b, a-2)]
                       -50.0*qre[idx (b, a-3)]+15.0*qre[idx (b, a-4)]
                       +77.0*qre[idx (b, a)])/(60.0*dx);
	qimx = (10.0*qim[idx (b, a+1)]-150.0*qim[idx (b, a-1)]
                       -2.0*qim[idx (b, a-5)]+100.0*qim[idx (b, a-2)]
                       -50.0*qim[idx (b, a-3)]+15.0*qim[idx (b, a-4)]
                       +77.0*qim[idx (b, a)])/(60.0*dx);
	prex = (10.0*pre[idx (b, a+1)]-150.0*pre[idx (b, a-1)]
                       -2.0*pre[idx (b, a-5)]+100.0*pre[idx (b, a-2)]
                       -50.0*pre[idx (b, a-3)]+15.0*pre[idx (b, a-4)]
                       +77.0*pre[idx (b, a)])/(60.0*dx);
	pimx = (10.0*pim[idx (b, a+1)]-150.0*pim[idx (b, a-1)]
                       -2.0*pim[idx (b, a-5)]+100.0*pim[idx (b, a-2)]
                       -50.0*pim[idx (b, a-3)]+15.0*pim[idx (b, a-4)]
                       +77.0*pim[idx (b, a)])/(60.0*dx);
        } else {
	qrex = (45.0*qre[idx (b, a+1)]-45.0*qre[idx (b, a-1)]
                       -9.0*tqrep[1][b]+9.0*qre[idx (b, a-2)]
                       +tqrep[2][b]-qre[idx (b, a-3)])/(60.0*dx);
	qimx = (45.0*qim[idx (b, a+1)]-45.0*qim[idx (b, a-1)]
                       -9.0*tqimp[1][b]+9.0*qim[idx (b, a-2)]
                       +tqimp[2][b]-qim[idx (b, a-3)])/(60.0*dx);
	prex = (45.0*pre[idx (b, a+1)]-45.0*pre[idx (b, a-1)]
                       -9.0*tprep[1][b]+9.0*pre[idx (b, a-2)]
                       +tprep[2][b]-pre[idx (b, a-3)])/(60.0*dx);
	pimx = (45.0*pim[idx (b, a+1)]-45.0*pim[idx (b, a-1)]
                       -9.0*tpimp[1][b]+9.0*pim[idx (b, a-2)]
                       +tpimp[2][b]-pim[idx (b, a-3)])/(60.0*dx);
        }
        }

       if (a == 0) {
        if (my_rank == 0) {
	qrex = (360.0*qre[idx (b, a+1)]-10.0*qre[idx (b, a+6)]
                       -450.0*qre[idx (b, a+2)]+72.0*qre[idx (b, a+5)]
                       +400.0*qre[idx (b, a+3)]-225.0*qre[idx (b, a+4)]
                       -147.0*qre[idx (b, a)])/(60.0*dx);
	qimx = (360.0*qim[idx (b, a+1)]-10.0*qim[idx (b, a+6)]
                       -450.0*qim[idx (b, a+2)]+72.0*qim[idx (b, a+5)]
                       +400.0*qim[idx (b, a+3)]-225.0*qim[idx (b, a+4)]
                       -147.0*qim[idx (b, a)])/(60.0*dx);
	prex = (360.0*pre[idx (b, a+1)]-10.0*pre[idx (b, a+6)]
                       -450.0*pre[idx (b, a+2)]+72.0*pre[idx (b, a+5)]
                       +400.0*pre[idx (b, a+3)]-225.0*pre[idx (b, a+4)]
                       -147.0*pre[idx (b, a)])/(60.0*dx);
	pimx = (360.0*pim[idx (b, a+1)]-10.0*pim[idx (b, a+6)]
                       -450.0*pim[idx (b, a+2)]+72.0*pim[idx (b, a+5)]
                       +400.0*pim[idx (b, a+3)]-225.0*pim[idx (b, a+4)]
                       -147.0*pim[idx (b, a)])/(60.0*dx);
        } else {
        qrex = (45.0*qre[idx (b, a+1)]-45.0*tqrem[1][b]
                       -9.0*qre[idx (b, a+2)]+9.0*tqrem[2][b]
                       +qre[idx (b, a+3)]-tqrem[3][b])/(60.0*dx);
        qimx = (45.0*qim[idx (b, a+1)]-45.0*tqimm[1][b]
                       -9.0*qim[idx (b, a+2)]+9.0*tqimm[2][b]
                       +qim[idx (b, a+3)]-tqimm[3][b])/(60.0*dx);
        prex = (45.0*pre[idx (b, a+1)]-45.0*tprem[1][b]
                       -9.0*pre[idx (b, a+2)]+9.0*tprem[2][b]
                       +pre[idx (b, a+3)]-tprem[3][b])/(60.0*dx);
        pimx = (45.0*pim[idx (b, a+1)]-45.0*tpimm[1][b]
                       -9.0*pim[idx (b, a+2)]+9.0*tpimm[2][b]
                       +pim[idx (b, a+3)]-tpimm[3][b])/(60.0*dx);
        }        
        }


        if (a == N-1) {
        if (my_rank == p-1){
	qrex = (10.0*qre[idx (b, a-6)]-360.0*qre[idx (b, a-1)]
                       -72.0*qre[idx (b, a-5)]+450.0*qre[idx (b, a-2)]
                       -400.0*qre[idx (b, a-3)]+225.0*qre[idx (b, a-4)]
                       +147.0*qre[idx (b, a)])/(60.0*dx);
	qimx = (10.0*qim[idx (b, a-6)]-360.0*qim[idx (b, a-1)]
                       -72.0*qim[idx (b, a-5)]+450.0*qim[idx (b, a-2)]
                       -400.0*qim[idx (b, a-3)]+225.0*qim[idx (b, a-4)]
                       +147.0*qim[idx (b, a)])/(60.0*dx);
	prex = (10.0*pre[idx (b, a-6)]-360.0*pre[idx (b, a-1)]
                       -72.0*pre[idx (b, a-5)]+450.0*pre[idx (b, a-2)]
                       -400.0*pre[idx (b, a-3)]+225.0*pre[idx (b, a-4)]
                       +147.0*pre[idx (b, a)])/(60.0*dx);
	pimx = (10.0*pim[idx (b, a-6)]-360.0*pim[idx (b, a-1)]
                       -72.0*pim[idx (b, a-5)]+450.0*pim[idx (b, a-2)]
                       -400.0*pim[idx (b, a-3)]+225.0*pim[idx (b, a-4)]
                       +147.0*pim[idx (b, a)])/(60.0*dx);
        } else {
	qrex = (45.0*tqrep[1][b]-45.0*qre[idx (b, a-1)]
                       -9.0*tqrep[2][b]+9.0*qre[idx (b, a-2)]
                       +tqrep[3][b]-qre[idx (b, a-3)])/(60.0*dx);
	qimx = (45.0*tqimp[1][b]-45.0*qim[idx (b, a-1)]
                       -9.0*tqimp[2][b]+9.0*qim[idx (b, a-2)]
                       +tqimp[3][b]-qim[idx (b, a-3)])/(60.0*dx);
	prex = (45.0*tprep[1][b]-45.0*pre[idx (b, a-1)]
                       -9.0*tprep[2][b]+9.0*pre[idx (b, a-2)]
                       +tprep[3][b]-pre[idx (b, a-3)])/(60.0*dx);
	pimx = (45.0*tpimp[1][b]-45.0*pim[idx (b, a-1)]
                       -9.0*tpimp[2][b]+9.0*pim[idx (b, a-2)]
                       +tpimp[3][b]-pim[idx (b, a-3)])/(60.0*dx);
        }
        }

        gdd_real  qredss = make_dd(0.0);
        gdd_real  qimdss = make_dd(0.0);
        gdd_real  predss = make_dd(0.0);
        gdd_real  pimdss = make_dd(0.0);

/*
        if ((a >= 3)&&(a < N-3)) {
	qrex = (45.0*qre[idx (b, a+1)]-45.0*qre[idx (b, a-1)]
                       -9.0*qre[idx (b, a+2)]+9.0*qre[idx (b, a-2)]
                       +qre[idx (b, a+3)]-qre[idx (b, a-3)])/(60.0*dx);
	qimx = (45.0*qim[idx (b, a+1)]-45.0*qim[idx (b, a-1)]
                       -9.0*qim[idx (b, a+2)]+9.0*qim[idx (b, a-2)]
                       +qim[idx (b, a+3)]-qim[idx (b, a-3)])/(60.0*dx);
	prex = (45.0*pre[idx (b, a+1)]-45.0*pre[idx (b, a-1)]
                       -9.0*pre[idx (b, a+2)]+9.0*pre[idx (b, a-2)]
                       +pre[idx (b, a+3)]-pre[idx (b, a-3)])/(60.0*dx);
	pimx = (45.0*pim[idx (b, a+1)]-45.0*pim[idx (b, a-1)]
                       -9.0*pim[idx (b, a+2)]+9.0*pim[idx (b, a-2)]
                       +pim[idx (b, a+3)]-pim[idx (b, a-3)])/(60.0*dx);
        }

         if ((a >= 4)&&(a < N-4)) {
         qredss = qre[idx (b,a+4)]+qre[idx (b,a-4)] - 8.0*(qre[idx (b,a+3)]+qre[idx (b,a-3)]) + 70.0*qre[idx (b,a)]
           - 56.0*(qre[idx (b,a+1)]+qre[idx (b,a-1)]) + 28.0*(qre[idx (b,a+2)]+qre[idx (b,a-2)]);   
         qimdss = qim[idx (b,a+4)]+qim[idx (b,a-4)] - 8.0*(qim[idx (b,a+3)]+qim[idx (b,a-3)]) + 70.0*qim[idx (b,a)]
           - 56.0*(qim[idx (b,a+1)]+qim[idx (b,a-1)]) + 28.0*(qim[idx (b,a+2)]+qim[idx (b,a-2)]);   
         predss = pre[idx (b,a+4)]+pre[idx (b,a-4)] - 8.0*(pre[idx (b,a+3)]+pre[idx (b,a-3)]) + 70.0*pre[idx (b,a)]
           - 56.0*(pre[idx (b,a+1)]+pre[idx (b,a-1)]) + 28.0*(pre[idx (b,a+2)]+pre[idx (b,a-2)]);   
         pimdss = pim[idx (b,a+4)]+pim[idx (b,a-4)] - 8.0*(pim[idx (b,a+3)]+pim[idx (b,a-3)]) + 70.0*pim[idx (b,a)]
           - 56.0*(pim[idx (b,a+1)]+pim[idx (b,a-1)]) + 28.0*(pim[idx (b,a+2)]+pim[idx (b,a-2)]);   
         }

         qredss = 0.1 * qredss / 16.0 / dt;
         qimdss = 0.1 * qimdss / 16.0 / dt;
         predss = 0.1 * predss / 16.0 / dt;
         pimdss = 0.1 * pimdss / 16.0 / dt;
*/

        if ((a >= 3)&&(a < N-3)) {
        f_mid_p = wenoL(qre[idx (b,a)], qre[idx (b,a-1)], qre[idx (b,a+1)], qre[idx (b,a-2)], qre[idx (b,a+2)]); 
        f_mid_m = wenoL(qre[idx (b,a-1)], qre[idx (b,a-2)], qre[idx (b,a)], qre[idx (b,a-3)], qre[idx (b,a+1)]); 
        qrex = (f_mid_p - f_mid_m) / dx;

        f_mid_p = wenoL(qim[idx (b,a)], qim[idx (b,a-1)], qim[idx (b,a+1)], qim[idx (b,a-2)], qim[idx (b,a+2)]); 
        f_mid_m = wenoL(qim[idx (b,a-1)], qim[idx (b,a-2)], qim[idx (b,a)], qim[idx (b,a-3)], qim[idx (b,a+1)]); 
        qimx = (f_mid_p - f_mid_m) / dx;

        f_mid_m = wenoR(pre[idx (b,a)], pre[idx (b,a-1)], pre[idx (b,a+1)], pre[idx (b,a-2)], pre[idx (b,a+2)]); 
        f_mid_p = wenoR(pre[idx (b,a+1)], pre[idx (b,a)], pre[idx (b,a+2)], pre[idx (b,a-1)], pre[idx (b,a+3)]); 
        prex = (f_mid_p - f_mid_m) / dx;

        f_mid_m = wenoR(pim[idx (b,a)], pim[idx (b,a-1)], pim[idx (b,a+1)], pim[idx (b,a-2)], pim[idx (b,a+2)]); 
        f_mid_p = wenoR(pim[idx (b,a+1)], pim[idx (b,a)], pim[idx (b,a+2)], pim[idx (b,a-1)], pim[idx (b,a+3)]); 
        pimx = (f_mid_p - f_mid_m) / dx;
        }

	rhs_qre[idx (b, a)] = 0.0-bb_c[idx (b, a)] * qrex + pre[idx (b, a)] -qredss;
	rhs_qim[idx (b, a)] = 0.0-bb_c[idx (b, a)] * qimx + pim[idx (b, a)] -qimdss;

        //printf ("%d %g %g\n", a, rhs_qre[idx (b, a)], rhs_qim[idx (b, a)]);

	rhs_pre[idx (b, a)] = 
	  ee_c[idx (b, a)] * prex                 -predss
	  + dre * qrex - d_im_c[idx (b, a)] * qimx
	  - a_re_c[idx (b, a)] * pre[idx (b, a)] + a_im_c[idx (b, a)] * pim[idx (b, a)]
	  - v_re_c[idx (b, a)] * qre[idx (b, a)] + v_im_c[idx (b, a)] * qim[idx (b, a)];

	rhs_pim[idx (b, a)] = 
	  ee_c[idx (b, a)] * pimx                 -pimdss
	  + dre * qimx + d_im_c[idx (b, a)] * qrex
	  - a_re_c[idx (b, a)] * pim[idx (b, a)] - a_im_c[idx (b, a)] * pre[idx (b, a)]
	  - v_re_c[idx (b, a)] * qim[idx (b, a)] - v_im_c[idx (b, a)] * qre[idx (b, a)];

      }
}

__global__ void
kernel_update1 (gdd_real *rhs_qre, gdd_real *rhs_qim, gdd_real *rhs_pre, gdd_real *rhs_pim, gdd_real *qre, gdd_real *qim, gdd_real *pre, gdd_real *pim, gdd_real *qre_c, gdd_real *qim_c, gdd_real *pre_c, gdd_real *pim_c)
{
  int b;
  int a;

  b = blockIdx.y * blockDim.y + threadIdx.y;
  a = blockIdx.x * blockDim.x + threadIdx.x;

  //if (b >= 5 && b < M - 5)
    if (a >= 0 && a < N - 0)
      {
	qre_c[idx (b, a)] = qre[idx (b, a)] + 0.5 * dt * rhs_qre[idx (b, a)];
	qim_c[idx (b, a)] = qim[idx (b, a)] + 0.5 * dt * rhs_qim[idx (b, a)];
	pre_c[idx (b, a)] = pre[idx (b, a)] + 0.5 * dt * rhs_pre[idx (b, a)];
	pim_c[idx (b, a)] = pim[idx (b, a)] + 0.5 * dt * rhs_pim[idx (b, a)];
      }
}

__global__ void
kernel_boundary_a (gdd_real *qre_c, gdd_real *qim_c, gdd_real *pre_c, gdd_real *pim_c)
{
  int a;

  a = blockIdx.x * blockDim.x + threadIdx.x;

}


__global__ void
kernel_boundary_b (gdd_real *qre_c, gdd_real *qim_c, gdd_real *pre_c, gdd_real *pim_c, gdd_real *d_buff, gdd_real *d_buff2)
{

  int b;
  b = blockIdx.x * blockDim.x + threadIdx.x;

  if (my_rank > 0) {

  d_buff[b] = qre_c[idx (b, 0)];
  d_buff[M+b] = qre_c[idx (b, 1)];
  d_buff[2*M+b] = qre_c[idx (b, 2)];

  d_buff[3*M+b] = qim_c[idx (b, 0)];
  d_buff[3*M+M+b] = qim_c[idx (b, 1)];
  d_buff[3*M+2*M+b] = qim_c[idx (b, 2)];

  d_buff[6*M+b] = pre_c[idx (b, 0)];
  d_buff[6*M+M+b] = pre_c[idx (b, 1)];
  d_buff[6*M+2*M+b] = pre_c[idx (b, 2)];

  d_buff[9*M+b] = pim_c[idx (b, 0)];
  d_buff[9*M+M+b] = pim_c[idx (b, 1)];
  d_buff[9*M+2*M+b] = pim_c[idx (b, 2)];

  }

  if (my_rank < p-1) {

  d_buff2[b] = qre_c[idx (b, N-1)];
  d_buff2[M+b] = qre_c[idx (b, N-2)];
  d_buff2[2*M+b] = qre_c[idx (b, N-3)];

  d_buff2[3*M+b] = qim_c[idx (b, N-1)];
  d_buff2[3*M+M+b] = qim_c[idx (b, N-2)];
  d_buff2[3*M+2*M+b] = qim_c[idx (b, N-3)];

  d_buff2[6*M+b] = pre_c[idx (b, N-1)];
  d_buff2[6*M+M+b] = pre_c[idx (b, N-2)];
  d_buff2[6*M+2*M+b] = pre_c[idx (b, N-3)];

  d_buff2[9*M+b] = pim_c[idx (b, N-1)];
  d_buff2[9*M+M+b] = pim_c[idx (b, N-2)];
  d_buff2[9*M+2*M+b] = pim_c[idx (b, N-3)];

  }

}


__global__ void
kernel_update2 (gdd_real *rhs_qre, gdd_real *rhs_qim, gdd_real *rhs_pre, gdd_real *rhs_pim, gdd_real *qre, gdd_real *qim, gdd_real *pre, gdd_real *pim)
{
  int b;
  int a;

  b = blockIdx.y * blockDim.y + threadIdx.y;
  a = blockIdx.x * blockDim.x + threadIdx.x;

  //if (b >= 5 && b < M - 5)
    if (a >= 0 && a < N - 0)
      {
	qre[idx (b, a)] = qre[idx (b, a)] + dt * rhs_qre[idx (b, a)];
	qim[idx (b, a)] = qim[idx (b, a)] + dt * rhs_qim[idx (b, a)];
	pre[idx (b, a)] = pre[idx (b, a)] + dt * rhs_pre[idx (b, a)];
	pim[idx (b, a)] = pim[idx (b, a)] + dt * rhs_pim[idx (b, a)];
      }


}


__global__ void
pull_data (gdd_real *qre, gdd_real *qim, double *d_qre_buff, double *d_qim_buff)
{
  int b;
  int a;

  b = blockIdx.y * blockDim.y + threadIdx.y;
  a = blockIdx.x * blockDim.x + threadIdx.x;

  d_qre_buff[idx (b, a)] = to_double(qre[idx (b, a)]);
  d_qim_buff[idx (b, a)] = to_double(qim[idx (b, a)]);

}



/* -------------------------------------------------- */
