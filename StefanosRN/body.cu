#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <cuda.h>
#include <mpi.h>
#include <qd/dd_real.h>
#include <qd/fpu.h>
#include <gqd.cu>
#include "parm.h"

__global__ void kernel_init (gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real, gdd_real, gdd_real, gdd_real, gdd_real, gdd_real, gdd_real, gdd_real *, double *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, int, int);
__global__ void kernel_rhs (gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *,  gdd_real *, gdd_real *);
__global__ void kernel_update1 (gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *);
__global__ void kernel_boundary_a (gdd_real *, gdd_real *, gdd_real *, gdd_real *);
__global__ void kernel_boundary_b (gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *);
__global__ void kernel_update2 (gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *, gdd_real *);
__global__ void pull_data (gdd_real *, gdd_real *, double *, double *);

double project (double *, double, double);
#include "kern.cu"

/* -------------------------------------------------- */

void
body (dd_real *r_c, dd_real *theta, dd_real dtheta, 
      dd_real mass, dd_real aa, dd_real ss, dd_real mm, dd_real dx,
      dd_real dt, int p, int my_rank)
{

  cudaSetDevice( my_rank%4 );

  int l; 
  gdd_real gmass, gaa, gss, gmm, gdx, gdtheta, gdt;
  double tr_c[N], ttheta[M], sbuff[M*12*2], rbuff[M*12*2];
  double qre_buff[M*N], qim_buff[M*N], qq, ell;
  double psi_re0[M], psi_im0[M], psi_re1[M], psi_im1[M], psi_re2[M], psi_im2[M];
  double psi_re3[M], psi_im3[M], psi_re4[M], psi_im4[M], psi_re5[M], psi_im5[M];

  ell = to_double(mm);
  if ( abs(ss) > abs(mm) ) 
  ell = to_double(ss);
  ell = fabs(ell);

  for (l = 0; l < N; l++){
    tr_c[l] = to_double(r_c[l]);
  }

  printf (" tr_c[N / 2] = %g \n",  tr_c[N / 2]);

  for (l = 0; l < M; l++)
    ttheta[l] = to_double(theta[l]);

  printf (" ttheta[M / 2] = %g \n",  ttheta[M / 2]);

  GDDStart( my_rank%4 );

        gmass.x = mass.x[0];
        gmass.y = mass.x[1];

        gaa.x = aa.x[0];
        gaa.y = aa.x[1];

        gss.x = ss.x[0];
        gss.y = ss.x[1];

        gmm.x = mm.x[0];
        gmm.y = mm.x[1];

        gdx.x = dx.x[0];
        gdx.y = dx.x[1];

        gdt.x = dt.x[0];
        gdt.y = dt.x[1];

        gdtheta.x = dtheta.x[0];
        gdtheta.y = dtheta.x[1];

  gdd_real *d_r_c, *d_theta, *d_buff, *d_buff2;
  double *d_qre_buff, *d_qim_buff, *d_tr_c, *d_ttheta;
  cudaMalloc ((void **) &d_r_c, N * sizeof (gdd_real));
  cudaMalloc ((void **) &d_theta, M * sizeof (gdd_real));
  cudaMalloc ((void **) &d_tr_c, N * sizeof (double));
  cudaMalloc ((void **) &d_ttheta, M * sizeof (double));
  cudaMalloc ((void **) &d_qre_buff, M * N * sizeof (double));
  cudaMalloc ((void **) &d_qim_buff, M * N * sizeof (double));
  cudaMalloc ((void **) &d_buff, M * 12 * sizeof (gdd_real));
  cudaMalloc ((void **) &d_buff2, M * 12 * sizeof (gdd_real));

//Dynamic memory experiment ... 

  gdd_real *d_qre_c, *d_qim_c, *d_pre_c, *d_pim_c;
  cudaMalloc ((void **) &d_qre_c, M * N * sizeof (gdd_real));
  cudaMalloc ((void **) &d_qim_c, M * N * sizeof (gdd_real));
  cudaMalloc ((void **) &d_pre_c, M * N * sizeof (gdd_real));
  cudaMalloc ((void **) &d_pim_c, M * N * sizeof (gdd_real));

  gdd_real *d_qre, *d_qim, *d_pre, *d_pim;
  cudaMalloc ((void **) &d_qre, M * N * sizeof (gdd_real));
  cudaMalloc ((void **) &d_qim, M * N * sizeof (gdd_real));
  cudaMalloc ((void **) &d_pre, M * N * sizeof (gdd_real));
  cudaMalloc ((void **) &d_pim, M * N * sizeof (gdd_real));

  gdd_real *d_rhs_qre, *d_rhs_qim, *d_rhs_pre, *d_rhs_pim;
  cudaMalloc ((void **) &d_rhs_qre, M * N * sizeof (gdd_real));
  cudaMalloc ((void **) &d_rhs_qim, M * N * sizeof (gdd_real));
  cudaMalloc ((void **) &d_rhs_pre, M * N * sizeof (gdd_real));
  cudaMalloc ((void **) &d_rhs_pim, M * N * sizeof (gdd_real));

  gdd_real *d_bb_c, *d_cc_c, *d_ee_c, *d_a_re_c, *d_a_im_c, *d_d_re_c, *d_d_im_c, *d_v_re_c, *d_v_im_c;
  cudaMalloc ((void **) &d_bb_c, M * N * sizeof (gdd_real));
  cudaMalloc ((void **) &d_cc_c, M * N * sizeof (gdd_real));
  cudaMalloc ((void **) &d_ee_c, M * N * sizeof (gdd_real));
  cudaMalloc ((void **) &d_a_re_c, M * N * sizeof (gdd_real));
  cudaMalloc ((void **) &d_a_im_c, M * N * sizeof (gdd_real));
  cudaMalloc ((void **) &d_d_re_c, M * N * sizeof (gdd_real));
  cudaMalloc ((void **) &d_d_im_c, M * N * sizeof (gdd_real));
  cudaMalloc ((void **) &d_v_re_c, M * N * sizeof (gdd_real));
  cudaMalloc ((void **) &d_v_im_c, M * N * sizeof (gdd_real));

  //End dynamic memory portion

  cudaMemcpy (d_tr_c, tr_c, N * sizeof (double), cudaMemcpyHostToDevice);
  cudaMemcpy (d_ttheta, ttheta, M * sizeof (double), cudaMemcpyHostToDevice);
  cudaMemcpy (d_r_c, r_c, N * sizeof (gdd_real), cudaMemcpyHostToDevice);
  cudaMemcpy (d_theta, theta, M * sizeof (gdd_real), cudaMemcpyHostToDevice);

  dim3 threads (256, 1);
  dim3 grid (N / threads.x, M / threads.y);

  dim3 threads_boundary (256);
  dim3 grid_boundary (N / threads_boundary.x);

  dim3 threads_boundary_m (1);
  dim3 grid_boundary_m (M / threads_boundary_m.x);

  kernel_init <<< grid, threads >>> (d_qre, d_qim, d_pre, d_pim, gdtheta, gmass, gaa, gss, gmm, gdx, gdt, d_theta, d_r_c, d_bb_c, d_cc_c, d_ee_c, d_a_re_c, d_a_im_c, d_d_re_c, d_d_im_c, d_v_re_c, d_v_im_c, my_rank, p);
  cudaThreadSynchronize ();

  double timer = 0.0;
  double start = timer;
  double endtime = T;
  int nt = to_int ((endtime - start) / dt);
  int k, j;
  char string[32];

  FILE* fpout[Nf];
  FILE* fprad;

  if (my_rank == 0){
  for (k = 0; k <= 4; k++){
     sprintf (string, "TimeDAT/%d", k);
     strcat (string, ".dat");
     fpout[k] = fopen(string,"w");
    }
  }

  if (my_rank == p-1){
  for (k = 5; k <= 5; k++){
     sprintf (string, "TimeDAT/%d", k);
     strcat (string, ".dat");
     fpout[k] = fopen(string,"w");
    }
  }

  for (k = 0; k < nt; k++)
    {

     // ........ MPI transfers begin ......

     if (my_rank > 0) {
      cudaMemcpy (sbuff, d_buff, M * 12 * sizeof (gdd_real),
		      cudaMemcpyDeviceToHost);
      MPI_Send(
      /* data         = */ sbuff, 
      /* count        = */ M * 12 * 2, 
      /* datatype     = */ MPI_DOUBLE, 
      /* destination  = */ my_rank-1, 
      /* tag          = */ 0, 
      /* communicator = */ MPI_COMM_WORLD);
     } 

     if (my_rank < p-1) {
      MPI_Recv(
      /* data         = */ rbuff, 
      /* count        = */ M * 12 * 2, 
      /* datatype     = */ MPI_DOUBLE, 
      /* source       = */ my_rank+1, 
      /* tag          = */ 0, 
      /* communicator = */ MPI_COMM_WORLD,
      /* status       = */ MPI_STATUS_IGNORE);
      cudaMemcpy (d_buff, rbuff, M * 12 * sizeof (gdd_real),
		      cudaMemcpyHostToDevice);
     } 

     if (my_rank < p-1) {
      cudaMemcpy (sbuff, d_buff2, M * 12 * sizeof (gdd_real),
		      cudaMemcpyDeviceToHost);
      MPI_Send(
      /* data         = */ sbuff, 
      /* count        = */ M * 12 * 2, 
      /* datatype     = */ MPI_DOUBLE, 
      /* destination  = */ my_rank+1, 
      /* tag          = */ 0, 
      /* communicator = */ MPI_COMM_WORLD);
     } 

     if (my_rank > 0) {
      MPI_Recv(
      /* data         = */ rbuff, 
      /* count        = */ M * 12 * 2, 
      /* datatype     = */ MPI_DOUBLE, 
      /* destination  = */ my_rank-1, 
      /* tag          = */ 0, 
      /* communicator = */ MPI_COMM_WORLD,
      /* status       = */ MPI_STATUS_IGNORE);
      cudaMemcpy (d_buff2, rbuff, M * 12 * sizeof (gdd_real),
		      cudaMemcpyHostToDevice);
     } 

     // ........ MPI transfers end ......

      kernel_rhs <<< grid, threads >>>(d_rhs_qre, d_rhs_qim, d_rhs_pre, d_rhs_pim, d_qre, d_qim, d_pre, d_pim, d_theta, d_r_c, timer, d_bb_c, d_cc_c, d_ee_c, d_a_re_c, d_a_im_c, d_d_re_c, d_d_im_c, d_v_re_c, d_v_im_c, d_buff, d_buff2);
      cudaThreadSynchronize ();

      kernel_update1 <<< grid, threads >>> (d_rhs_qre, d_rhs_qim, d_rhs_pre, d_rhs_pim, d_qre, d_qim, d_pre, d_pim, d_qre_c, d_qim_c, d_pre_c, d_pim_c);
      cudaThreadSynchronize ();

      kernel_boundary_a <<< grid_boundary, threads_boundary >>> (d_qre_c, d_qim_c, d_pre_c, d_pim_c);
      cudaThreadSynchronize ();

      kernel_boundary_b <<< grid_boundary_m, threads_boundary_m >>> (d_qre_c, d_qim_c, d_pre_c, d_pim_c, d_buff, d_buff2);
      cudaThreadSynchronize ();

     // ........ MPI transfers begin ......

     if (my_rank > 0) {
      cudaMemcpy (sbuff, d_buff, M * 12 * sizeof (gdd_real),
		      cudaMemcpyDeviceToHost);
      MPI_Send(
      /* data         = */ sbuff, 
      /* count        = */ M * 12 * 2, 
      /* datatype     = */ MPI_DOUBLE, 
      /* destination  = */ my_rank-1, 
      /* tag          = */ 0, 
      /* communicator = */ MPI_COMM_WORLD);
     } 

     if (my_rank < p-1) {
      MPI_Recv(
      /* data         = */ rbuff, 
      /* count        = */ M * 12 * 2, 
      /* datatype     = */ MPI_DOUBLE, 
      /* source       = */ my_rank+1, 
      /* tag          = */ 0, 
      /* communicator = */ MPI_COMM_WORLD,
      /* status       = */ MPI_STATUS_IGNORE);
      cudaMemcpy (d_buff, rbuff, M * 12 * sizeof (gdd_real),
		      cudaMemcpyHostToDevice);
     } 

     if (my_rank < p-1) {
      cudaMemcpy (sbuff, d_buff2, M * 12 * sizeof (gdd_real),
		      cudaMemcpyDeviceToHost);
      MPI_Send(
      /* data         = */ sbuff, 
      /* count        = */ M * 12 * 2, 
      /* datatype     = */ MPI_DOUBLE, 
      /* destination  = */ my_rank+1, 
      /* tag          = */ 0, 
      /* communicator = */ MPI_COMM_WORLD);
     } 

     if (my_rank > 0) {
      MPI_Recv(
      /* data         = */ rbuff, 
      /* count        = */ M * 12 * 2, 
      /* datatype     = */ MPI_DOUBLE, 
      /* destination  = */ my_rank-1, 
      /* tag          = */ 0, 
      /* communicator = */ MPI_COMM_WORLD,
      /* status       = */ MPI_STATUS_IGNORE);
      cudaMemcpy (d_buff2, rbuff, M * 12 * sizeof (gdd_real),
		      cudaMemcpyHostToDevice);
     } 

     // ........ MPI transfers end ......

      kernel_rhs <<< grid, threads >>>(d_rhs_qre, d_rhs_qim, d_rhs_pre, d_rhs_pim, d_qre_c, d_qim_c, d_pre_c, d_pim_c, d_theta, d_r_c, timer, d_bb_c, d_cc_c, d_ee_c, d_a_re_c, d_a_im_c, d_d_re_c, d_d_im_c, d_v_re_c, d_v_im_c, d_buff, d_buff2);
      cudaThreadSynchronize ();

      kernel_update2 <<< grid, threads >>> (d_rhs_qre, d_rhs_qim, d_rhs_pre, d_rhs_pim, d_qre, d_qim, d_pre, d_pim);
      cudaThreadSynchronize ();

      kernel_boundary_a <<< grid_boundary, threads_boundary >>>(d_qre, d_qim, d_pre, d_pim);
      cudaThreadSynchronize ();

      kernel_boundary_b <<< grid_boundary_m, threads_boundary_m >>>(d_qre, d_qim, d_pre, d_pim, d_buff, d_buff2);
      cudaThreadSynchronize ();

      if (!(k % (int) Ft)){

      pull_data <<< grid, threads >>> (d_qre, d_qim, d_qre_buff, d_qim_buff);
      cudaThreadSynchronize ();

      cudaMemcpy (qre_buff, d_qre_buff, M * N * sizeof (double),
		      cudaMemcpyDeviceToHost);
      cudaMemcpy (qim_buff, d_qim_buff, M * N * sizeof (double),
		      cudaMemcpyDeviceToHost);

      dd_real tmp;
      tmp = mass + sqrt(mass*mass-aa*aa);
      tmp = tmp * Xmax / (tmp + Xmax);
      qq = 1.0;

      for (int r = 0; r < M; r++){

      j = 0 + to_int((tmp - Xmin)/dx + 0.5);
    
      psi_re0[r] = qre_buff[idx (r, j)];
      psi_im0[r] = qim_buff[idx (r, j)];

      psi_re1[r] = (-147.0*qre_buff[idx (r, j)] + 360.0*qre_buff[idx (r, j+1)] - 450.0*qre_buff[idx (r, j+2)] + 400.0*qre_buff[idx (r, j+3)] - 225.0*qre_buff[idx (r, j+4)] + 72.0*qre_buff[idx (r, j+5)] - 10.0*qre_buff[idx (r, j+6)] ) / (60.0*to_double(dx)) / qq;
      psi_im1[r] = (-147.0*qim_buff[idx (r, j)] + 360.0*qim_buff[idx (r, j+1)] - 450.0*qim_buff[idx (r, j+2)] + 400.0*qim_buff[idx (r, j+3)] - 225.0*qim_buff[idx (r, j+4)] + 72.0*qim_buff[idx (r, j+5)] - 10.0*qim_buff[idx (r, j+6)] ) / (60.0*to_double(dx)) / qq;

      psi_re2[r] = (469.0*qre_buff[idx (r, j)] - 2007.0*qre_buff[idx (r, j+1)] + 3955.5*qre_buff[idx (r, j+2)] - 4745.0*qre_buff[idx (r, j+3)] + 3690.0*qre_buff[idx (r, j+4)] - 1809.0*qre_buff[idx (r, j+5)] + 509.5*qre_buff[idx (r, j+6)] - 63.0*qre_buff[idx (r, j+7)] ) / (90.0*to_double(dx*dx)) / qq;
      psi_im2[r] = (469.0*qim_buff[idx (r, j)] - 2007.0*qim_buff[idx (r, j+1)] + 3955.5*qim_buff[idx (r, j+2)] - 4745.0*qim_buff[idx (r, j+3)] + 3690.0*qim_buff[idx (r, j+4)] - 1809.0*qim_buff[idx (r, j+5)] + 509.5*qim_buff[idx (r, j+6)] - 63.0*qim_buff[idx (r, j+7)] ) / (90.0*to_double(dx*dx)) / qq;

      psi_re3[r] = (-2403.0*qre_buff[idx (r, j)] + 13960.0*qre_buff[idx (r, j+1)] - 36706.0*qre_buff[idx (r, j+2)] + 57384.0*qre_buff[idx (r, j+3)] - 58280.0*qre_buff[idx (r, j+4)] + 39128.0*qre_buff[idx (r, j+5)] - 16830.0*qre_buff[idx (r, j+6)] + 4216.0*qre_buff[idx (r, j+7)] - 469.0*qre_buff[idx (r, j+8)]  ) / (240.0*to_double(dx*dx*dx)) / qq;
      psi_im3[r] = (-2403.0*qim_buff[idx (r, j)] + 13960.0*qim_buff[idx (r, j+1)] - 36706.0*qim_buff[idx (r, j+2)] + 57384.0*qim_buff[idx (r, j+3)] - 58280.0*qim_buff[idx (r, j+4)] + 39128.0*qim_buff[idx (r, j+5)] - 16830.0*qim_buff[idx (r, j+6)] + 4216.0*qim_buff[idx (r, j+7)] - 469.0*qim_buff[idx (r, j+8)]  ) / (240.0*to_double(dx*dx*dx)) / qq;

      //psi_re4[r] = (4275.0*qre_buff[idx (r, j)] - 30668.0*qre_buff[idx (r, j+1)] + 99604.0*qre_buff[idx (r, j+2)] - 192624.0*qre_buff[idx (r, j+3)] + 244498.0*qre_buff[idx (r, j+4)] - 210920.0*qre_buff[idx (r, j+5)] + 123348.0*qre_buff[idx (r, j+6)] - 47024.0*qre_buff[idx (r, j+7)] + 10579.0*qre_buff[idx (r, j+8)] - 1068.0*qre_buff[idx (r, j+9)]  ) / (240.0*to_double(dx*dx*dx*dx)) / qq;
      //psi_im4[r] = (4275.0*qim_buff[idx (r, j)] - 30668.0*qim_buff[idx (r, j+1)] + 99604.0*qim_buff[idx (r, j+2)] - 192624.0*qim_buff[idx (r, j+3)] + 244498.0*qim_buff[idx (r, j+4)] - 210920.0*qim_buff[idx (r, j+5)] + 123348.0*qim_buff[idx (r, j+6)] - 47024.0*qim_buff[idx (r, j+7)] + 10579.0*qim_buff[idx (r, j+8)] - 1068.0*qim_buff[idx (r, j+9)] ) / (240.0*to_double(dx*dx*dx*dx)) / qq;
     
      j = 0 + to_int((2.0 - Xmin)/dx + 0.5);
      psi_re4[r] = qre_buff[idx (r, j)];
      psi_im4[r] = qim_buff[idx (r, j)];

      psi_re5[r] = qre_buff[idx (r, N-1)];
      psi_im5[r] = qim_buff[idx (r, N-1)];

      }

         if (my_rank == 0){
         fprintf(fpout[0],"%f %g %g \n",timer,psi_re0[0],psi_im0[0]);
         fprintf(fpout[1],"%f %g %g \n",timer,psi_re1[0],psi_im1[0]);
         fprintf(fpout[2],"%f %g %g \n",timer,psi_re2[0],psi_im2[0]);
         fprintf(fpout[3],"%f %g %g \n",timer,psi_re3[0],psi_im3[0]);
         fprintf(fpout[4],"%f %g %g \n",timer,psi_re4[0],psi_im4[0]);
         printf("%f %g \n", timer,psi_re0[0]);
         }

         sprintf (string, "SpaceDAT/%d", (int)(timer + 0.5) + (my_rank+1) * 100000);
         strcat (string, ".dat");
         fprad = fopen(string,"w");

         for (int q = 0; q < N; q = q + 100)
         fprintf(fprad,"%f %g %g \n",tr_c[q],qre_buff[idx (0, q)],qim_buff[idx (0, q)]);

         fclose(fprad);

         if (my_rank == p-1)
         fprintf(fpout[5],"%f %g %g \n",timer,psi_re5[0],psi_im5[0]);
                  
    }

      timer += to_double(dt);

    }

    GDDEnd();

    if (my_rank == 0){
    for (int q = 0; q <= 5; q++)
    fclose (fpout[q]);
    }

    if (my_rank == p-1)
    fclose (fpout[5]);

    cudaFree (d_r_c);
    cudaFree (d_theta);
    cudaFree (d_tr_c);
    cudaFree (d_ttheta);
    cudaFree (d_qre_buff);
    cudaFree (d_qim_buff);
    cudaFree (d_buff);
    cudaFree (d_buff2);

}

