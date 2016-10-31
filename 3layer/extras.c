/**
  \file 3layer/extras.c

  Extra functions that are useful for 3layer.

  \author David Lewis, CABI, NKI
  \copyright GNU Public License
  \date 2012-2013

 */

// includes
#include "header.h"

/** 
  \brief Print an error message to stderr and exit with status EXIT_FAILURE.

  \param [in] errorstring String containing error message (variable-length argument list, like printf)
 */
void error(char *errorstring, ...)
{
	va_list ap;

	va_start(ap, errorstring);
	fprintf(stderr, "Error: ");
	vfprintf(stderr, errorstring, ap);
	fprintf(stderr, "\n");
	va_end(ap);

	exit(EXIT_FAILURE);
}


/**
  \brief Print usage message and exit with status EXIT_FAILURE.

  \param [in] program String containing the name of the program
 */

void print_usage(char *program)
{
    fprintf(stderr, "Usage: %s [options] <input_file> \n\n", program);
    fprintf(stderr, 
		"Reads an input parameter file <input_file> with parameters \n"
		"from a layered environment (3 layers: SR, SP, and SO of the CA1 \n"
		"region of the hippocampus) and calculates a diffusion curve.\n"
		"The input filename extension can be anything, but don't use \".dat\".\n"
		"The output file has the same name as the input file but with the \n"
		"extension \".dat\".\n\n"
		);
    fprintf(stderr, "z positions are relative to the source at z=0.\n\n");

    fprintf(stderr, "Options:\n"
        "\t-h, --help              print usage message\n"
        "\t-v, --verbose           be verbose\n"
		"\t-g, --global_kappa      use the same kappa in all layers (= kappa_sp)\n"
        "\t--nr <nr>               specify nr\n"
        "\t--nz <nr>               specify nz\n"
		"\t--nt <nt>               specify nt\n"
		"\t--nt_scale <factor>     specify scale factor for nt\n"
		"\t--probe_z <probe_z>     specify probe_z (in microns)\n"
		"\t--probe_r <probe_r>     specify probe_r (in microns)\n"
		"\t--ez1 <ez1>             specify z-position of bottom of cylinder (<0)\n"
		"\t--ez2 <ez2>             specify z-position of top of cylinder (>0)\n"
        "\t--alpha_so <alpha_so>   specify alpha_so\n"
        "\t--alpha_sp <alpha_sp>   specify alpha_sp\n"
        "\t--alpha_sr <alpha_sr>   specify alpha_sr\n"
        "\t--theta_so <theta_so>   specify theta_so\n"
        "\t--theta_sp <theta_sp>   specify theta_sp\n"
        "\t--theta_sr <theta_sr>   specify theta_sr\n"
		"\t--kappa_so <kappa_so>   specify kappa_so\n"
		"\t--kappa_sp <kappa_sp>   specify kappa_sp\n"
		"\t--kappa_sr <kappa_sr>   specify kappa_sr\n"
		"\t--kappa_outside <k_out> specify kappa_outside (mutually excl. with -g)\n"
        "\t--alpha_start <a_start> specify initial guess for apparent alpha\n"
        "\t--theta_start <t_start> specify initial guess for apparent theta\n"
        "\t--alpha_step <a_step>   specify initial step for apparent alpha\n"
        "\t--theta_step <t_step>   specify initial step for apparent theta\n"
        "\t--tmax <tmax>           specify total duration of experiment\n"
        "\t--fit_tol <fit_tol>     specify stopping criterion (simplex size)\n"
        "\t--itermax <itermax>     specify stopping criterion (max iterations)\n"
        "\t--outfile <outfile>     specify output file (parameters and curves)\n"
        "\t--pathfile <pathfile>   specify simplex path output file (just \n"
        "\t                        one vertex of the simplex per iteration)\n"
        "\t--images <basename>     specify basename of output conc images\n"
        "\t--image_spacing <delta> specify time spacing between output images\n"
        "\t--additional_sources \"<string>\" specify additional sources\n"
		"\t    <string> = <num_additional_sources> <source_params>\n"
		"\t    <source_params> = <sz1> <sr1> <crnt1> [<sz2> <sr2> <crnt2> ...]\n"
        );
    exit(EXIT_FAILURE);
}


/** 
  \brief Create an array of doubles of the specified size.

   Also initialize the array to all 0.

  \param [in] N Number of elements in array
  \param [in] string String containing a description of the array (for debugging purposes)

  \return Pointer to the new array
 */
double *create_array(int N, char *string)
{
	double *a;
	int i;

	a = (double *)malloc(sizeof(double) * N);
	if (a == NULL)
		error("Cannot allocate memory for %s array", string);

	for (i=0; i<N; i++)
		a[i] = 0.0;

	return a;
}

