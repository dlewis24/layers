/** 
  \file fit-layer/model.c

  Function for solving the forward problem with the 3-layer model.

  \author David Lewis, CABI, NKI
  \copyright GNU Public License
  \date 2012-2013

 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "header.h"


/**
  \brief Calculates the concentration as a function of space and 
  time and returns the probe concentration as a function of time.

  This function solves the diffusion equation in each layer \f$ k \f$ ,

\f[
\frac{\partial c_k(\vec r, t)}{\partial t} 
  = D_{free} \theta_k \nabla^2{c_k(\vec r, t)} 
  + s/\alpha_k - \kappa_k c_k(\vec r, t)
  \quad ,
\f]

subject to the continuity conditions at the interfaces between layers

\f[
c(\vec r_k, t) = c(\vec r_{k+1}, t)
\quad 
\f]

\f[
\alpha_k \theta_k \nabla c(\vec r_k, t) 
= \alpha_{k+1} \theta_{k+1} \nabla c(\vec r_{k+1}, t)
\quad 
\f]

and the boundary condition \f$ c(\vec r, t) = 0 \f$ (total 
absorption) at the top, the bottom, and the side of the cylinder, 
where

	\f$ c_k(\vec r, t) \f$ = concentration in layer \f$ k \f$

	\f$ D_{free} \f$ = free diffusion coefficient

	\f$ \theta_k \f$ = permeability in layer \f$ k \f$

	\f$ s_k(\vec r, t) \f$ = source in layer \f$ k \f$

	\f$ \alpha_k \f$ = extracellular volume fraction in layer \f$ k \f$

	\f$ \kappa_k \f$ = nonspecific clearance factor in layer \f$ k \f$

This function calls convolve3() to compute the Laplacian in cylindrical 
coordinates.

  \param[in] nt Number of support points in time
  \param[in] nz Number of rows of concentration matrix
  \param[in] nr Number of columns of concentration matrix
  \param[in] iprobe z-index of probe location
  \param[in] jprobe r-index of probe location 
  \param[in] iz1 z-index of SR-SP boundary
  \param[in] iz2 z-index of SP-SO boundary
  \param[in] nolayer Flag for whether to model a single homogeneous environment (true) or to model a 3-layer environment
  \param[in] dt Spacing in time (\f$ t_{i+1} = t_i + \Delta t \f$)
  \param[in] dr Spacing in r (in this program, \f$ \Delta z = \Delta r \f$)
  \param[in] sd Source delay (time before source starts)
  \param[in] st Duration of source
  \param[in] alpha_so Extracellular volume fraction in SO layer
  \param[in] theta_so Permeability in SO layer
  \param[in] kappa_so Nonspecific clearance factor in SO layer
  \param[in] alpha_sp Extracellular volume fraction in SP layer
  \param[in] theta_sp Permeability in SP layer
  \param[in] kappa_sp Nonspecific clearance factor in SP layer
  \param[in] alpha_sr Extracellular volume fraction in SR layer
  \param[in] theta_sr Permeability in SR layer
  \param[in] kappa_sr Nonspecific clearance factor in SR layer
  \param[in] dfree Free diffusion coefficient
  \param[in] t Time array 
  \param[in] s Source array
  \param[in] invr Array for \f$ 1/r \f$ values
  \param[out] p Probe array (concentration as a function of time)
 */

void calc_diffusion_curve_layer_fit_layer(int nt, int nz, int nr, int iprobe, int jprobe, int iz1, int iz2, int nolayer, double dt, double dr, double sd, double st, double alpha_so, double theta_so, double kappa_so, double alpha_sp, double theta_sp, double kappa_sp, double alpha_sr, double theta_sr, double kappa_sr, double dfree, double *t, double *s, double *invr, double *p)
{
	int i, j, k;
	double dstar_so = theta_so * dfree;
	double dstar_sp = theta_sp * dfree;
	double dstar_sr = theta_sr * dfree;
	double const_so1 = dstar_so * dt / SQR(dr);
	double const_so2 = dstar_so * dt / (2.0 * dr);
	double const_sp1 = dstar_sp * dt / SQR(dr);
	double const_sp2 = dstar_sp * dt / (2.0 * dr);
	double const_sr1 = dstar_sr * dt / SQR(dr);
	double const_sr2 = dstar_sr * dt / (2.0 * dr);

	/* Arrays internal to this function */
	double *c;   	/* concentration */
	double *c_sr;	/* concentration (SR layer) */
	double *c_sp;  	/* concentration (SP layer)*/
	double *c_so;	/* concentration (SO layer)*/
	double *cb_sr;	/* Average of rows at SR layer boundary */
	double *cb_so;	/* Average of rows at SO layer boundary */
	double *dc;   	/* Change in concentration */
	double *dc_sr; 
	double *dc_sp; 
	double *dc_so;

	/* Arrays for concentrations
	   r=0 is at c[1,*] 
	   c[0,*] is for extended symmetrical values */
	c = create_array(nz*(nr+1), "concentration");

	/* Arrays for concentration changes (d*) and for layers */
	dc = create_array(nz*(nr+1), "dc");

	cb_sr = create_array(nr+1, "cb_sr");
	cb_so = create_array(nr+1, "cb_so");

	c_sr = create_array((iz1+2)*(nr+1), "c_sr");
	dc_sr = create_array((iz1+2)*(nr+1), "dc_sr");

	c_sp = create_array((iz2-iz1+2)*(nr+1), "c_sp");
	dc_sp = create_array((iz2-iz1+2)*(nr+1), "dc_sp");

	c_so = create_array((nz-iz2)*(nr+1), "c_so");
	dc_so = create_array((nz-iz2)*(nr+1), "dc_so");


	/* Initialize the concentration at t=0 */
	for (i=0; i<nz*(nr+1); i++)
		c[i] = s[i];

	/* Source delay: Concentration = 0 for sd seconds */
	int nds = lround(sd/dt);
	if (nds >= nt)
		error("nds=%d, nt=%d. Delay start should be < total expt time", 
			nds, nt);
	for (k=0; k<nds; k++) 
		p[k] = 0.0;


	/* Loop over time */
	for (k=nds; k<nt; k++) {
		p[k] = c[INDEX(iprobe,jprobe)];    	/* record c at time t[k] */

		/* Calculate c at time t[k+1] = t[k] + dt */
		if (nolayer == FALSE) {	/* 3-layer model */
			/*
			 * The SR/SP/SO arrays extend in z one more row 
			 * to include extrapolated boundary values (the true 
			 * boundary is right in the middle of node positions)
			 */
			for (j=0; j<nr+1; j++) {
				cb_sr[j] = (   dstar_sr * alpha_sr * c[INDEX(iz1,j)] 
				             + dstar_sp * alpha_sp * c[INDEX(iz1+1,j)]   )
				           / ( dstar_sr * alpha_sr + dstar_sp * alpha_sp );

				cb_so[j] = (   dstar_sp * alpha_sp * c[INDEX(iz2,j)] 
				             + dstar_so * alpha_so * c[INDEX(iz2+1,j)]   )
				           / ( dstar_sp * alpha_sp + dstar_so * alpha_so );
			}

			for (j=0; j<nr+1; j++) {
				for (i=0; i<iz1+1; i++) 
					c_sr[INDEX(i,j)] = c[INDEX(i,j)];
				c_sr[INDEX(iz1+1,j)] = 2.0 * cb_sr[j] - c[INDEX(iz1,j)];

				c_sp[INDEX(0,j)] = 2.0 * cb_sr[j] - c[INDEX(iz1+1,j)];
				for (i=iz1+1; i<iz2+1; i++) 
					c_sp[INDEX(i-iz1,j)] = c[INDEX(i,j)];
				c_sp[INDEX(iz2-iz1+1,j)] = 2.0 * cb_so[j] - c[INDEX(iz2,j)];

				c_so[INDEX(0,j)] = 2.0 * cb_so[j] - c[INDEX(iz2+1,j)];
				for (i=iz2+1; i<nz; i++) 
					c_so[INDEX(i-iz2,j)] = c[INDEX(i,j)];
			}

			/* Calculate the delta-c matrices */
		    convolve3(iz1+2, nr+1, c_sr, const_sr1, const_sr2, invr, dc_sr);

		    convolve3(iz2-iz1+2, nr+1, c_sp, const_sp1, const_sp2, invr, 
			                                                     dc_sp);
		    convolve3(nz-iz2, nr+1, c_so, const_so1, const_so2, invr, dc_so);

			/* Update the concentration matrix */
			for (j=0; j<nr+1; j++) {
				for (i=0; i<iz1+1; i++) 
					c[INDEX(i,j)] = c_sr[INDEX(i,j)] 
					             + dc_sr[INDEX(i,j)];

				for (i=iz1+1; i<iz2+1; i++) 
					c[INDEX(i,j)] = c_sp[INDEX(i-iz1,j)] 
					             + dc_sp[INDEX(i-iz1,j)];

				for (i=iz2+1; i<nz; i++) 
					c[INDEX(i,j)] = c_so[INDEX(i-iz2,j)]
					             + dc_so[INDEX(i-iz2,j)];
			}

		} else {	/* 1-layer model */
			/* Print a warning to the user */
			if (k==0) printf("\nNOTE: nolayer = %d, "
							"so using the 1 layer model\n\n", nolayer);

			/* Calculate the delta-c matrix */
		    convolve3(nz, nr+1, c, const_sr1, const_sr2, invr, dc);

			/* Update the concentration matrix */
			for (i=0; i<nz*(nr+1); i++)
				c[i] += dc[i];

		}


		/* If t < st, the source gets added to the concentration matrix 
	       for the next time-step */
		if (t[k] + dt/2.0 < sd + st) {
			for (i=0; i<nz*(nr+1); i++)
				c[i] += s[i];
		} 

		/* Model the non-specific clearance */
		for (j=0; j<nr+1; j++) {
			for (i=0; i<iz1+1; i++) 
				c[INDEX(i,j)] *= (1. - kappa_sr * dt);
			for (i=iz1+1; i<iz2+1; i++) 
				c[INDEX(i,j)] *= (1. - kappa_sp * dt);
			for (i=iz2+1; i<nz; i++) 
				c[INDEX(i,j)] *= (1. - kappa_so * dt);
		}

		/* Set the i=0 row to be the same as the i=2 row 
		   (symmetry about r=0 (i=1)) */ 
		for (i=0; i<nz; i++)
			c[INDEX(i,0)] = c[INDEX(i,2)];

	} /* End of k for loop */


	/* Deallocate arrays */
	free(c);
	free(dc);
	free(cb_sr);
	free(cb_so);
	free(c_sr);
	free(c_sp);
	free(c_so);
	free(dc_sr);
	free(dc_sp);
	free(dc_so);

	return;
}
