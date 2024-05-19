#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qd/dd_real.h>
#include <qd/fpu.h>
#include "parm.h"


/* Physical Parameters */

dd_real mass, aa, ss, mm;

/* grid */

dd_real x_c[N], r_c[N], theta[M];

/* miscellaneous */

dd_real dx, dtheta, dt, delta;

void grid (int, int);		/* set up the grid */
void body (dd_real*, dd_real* , dd_real, dd_real, 
           dd_real, dd_real, dd_real, dd_real, dd_real, int, int);

int
main (int argc, char **argv)
{				/* Main Code */

  int p, my_rank;

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &p);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  //unsigned int old_cw;
  //fpu_fix_start(&old_cw);

  if (my_rank==0) {
  printf ("\n");
  printf ("WELCOME TO TEUKOLSKY CODE!\n");
  printf ("\n");
  }

  mass = "1.0";
  aa = "1.0";
  ss = "1e-60";
  mm = "1e-60";

  dx = (dd_real)((Xmax - Xmin) / N / p); /* Radial grid increment */
  dtheta = (dd_real)(-2.0 / M);	 /* Angular grid increment */
  dt = 0.25 * dx;

  grid (my_rank, p);			 /* Setting up the grid */

  body (r_c, theta, dtheta, mass, aa, ss, mm, dx,
	dt, p, my_rank);

  if (my_rank==0) {
  printf ("All done.\n");
  printf ("Bye bye!\n");
  }

  MPI_Finalize ();

}

void
grid (int my_rank, int p)
{

  int l, m;
  dd_real tol, diff, drdx, rp, rm, r_old, r_new, x_o;
  dd_real rpt, rin, Rrpt;
 
  for (l = 0; l < N; l++)
    {				/* Setting up the grid */

      x_c[l] = Xmin + ((Xmax-Xmin)/p)*my_rank + (dd_real) (l-0) * dx;
      r_c[l] = Xmin + ((Xmax-Xmin)/p)*my_rank + (dd_real) (l-0) * dx;

    }


  for (m = 0; m < M; m++)
    {

      theta[m] = 1.0 + (dd_real) m *dtheta;

    }

}
