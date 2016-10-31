/**
  \file 3layer/rti-theory.c

  Functions for fitting the homogeneous model to the data 
  (traditional fit) to obtain the apparent parameters and 
  the characteristic curve.

  \author David Lewis, CABI, NKI
  \copyright GNU Public License
  \date 2012-2013

 */

#include "header.h"

/**
  \brief Mean squared error function for simplex fitting 

  The 3layer program uses a minimization function in GSL to 
  minimize the mean squared error between the diffusion curve 
  calculated from the multilayer model and a diffusion curve 
  calculated from a formula that assumes an isotropic 
  homogeneous environment. It does so to calculate the 
  apparent parameters and characteristic diffusion curve 
  (for comparison purposes).

  This function calculates the mean squared error between 
  the two curves. The GSL minimization routine needs the 
  first parameter to be a vector of doubles which represent 
  the parameters to be fit and the second parameter to be a 
  pointer to parameters of the function to fit (rti_theory). 

  The second parameter actually points to a struct of parameters. 
  The characteristic curve calculated from rti_theory is stored 
  in the p_theory array of this struct.

  \param[in] x Vector of doubles representing parameters to fit (alpha and theta)
  \param[in,out] params Vector of parameters needed by rti_theory

  \return Mean squared error between multilayer curve and characteristic curve
 */
double calc_mse_rti(const gsl_vector *x, void *params)
{
	int i = -1;
	mse_rti_params_struct_type *p = (mse_rti_params_struct_type *) params;
	int nt = p->nt;


	p->alpha = gsl_vector_get(x, 0);
	p->theta = gsl_vector_get(x, 1);
	if (p->alpha <= 0.001) p->alpha = 0.001;
	if (p->theta <= 0.001) p->theta = 0.001;
   
	/* Call the rti_theory function */
	rti_theory(nt, p->spdist, p->samplitude, p->sdelay, p->sduration, p->kappa, p->dfree, 
	           p->alpha, p->theta, p->t, p->p_theory);

	double mse = 0.;

	for (i=1; i<nt; i++) {
		mse += SQR(p->p_model[i] - p->p_theory[i]);
	}
	mse /= nt;

	return mse;
}


/**
  \brief Calculates RTI data for diffusion in an isotropic, 
         homogeneous environment (direct calculation from an equation)

  \todo Change formula to take kappa into account

  \param[in] nt Number of time points of calculation
  \param[in] spdist Distance between source and probe
  \param[in] samplitude Amplitude of source
  \param[in] sdelay Source delay (time before source starts)
  \param[in] sduration Duration of source
  \param[in] kappa Nonspecific clearance factor 
  \param[in] dfree Free diffusion coefficient
  \param[in] alpha Extracellular volume fraction 
  \param[in] theta Permeability 
  \param[in] t Time array
  \param[out] p_theory Probe array (concentration as a function of time)

 */
void rti_theory(int nt, double spdist, double samplitude, double sdelay, double sduration, double kappa, double dfree, double alpha, double theta, double *t, double *p_theory)
{
	int i;

	double dstar = theta * dfree;
	double ampl = samplitude / (4.0 * PI * alpha * dstar * spdist);

	for (i=0; i<nt; i++) {
		if (t[i] <= sdelay) {
			p_theory[i] = 0.;
		} else {
			p_theory[i] = ampl * erfc( spdist /
									(2.0 * sqrt(dstar * (t[i] - sdelay))) );
			if (t[i] > sdelay + sduration)
				p_theory[i] = p_theory[i] - ampl * erfc( spdist /
							(2.0 * sqrt(dstar * (t[i] - (sdelay + sduration)))) );
		}
	}
}

