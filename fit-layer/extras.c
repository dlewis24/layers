/**
  \file fit-layer/extras.c

  Extra functions that are useful for fit-layer.

  \author David Lewis, CABI, NKI
  \copyright GNU Public License
  \date 2012-2013

 */

// includes
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
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
void print_usage_fit_layer(char *program)
{
    fprintf(stderr, "Usage: %s [options] <input_file> \n\n", program);
    fprintf(stderr, 
		"Reads an input parameter/data file <input_file> with parameters and \n"
		"RTI data from a layered environment (3 layers: SR, SP, and SO of the \n"
		"CA1 region of the hippocampus) and fits alpha, theta, and kappa of SP.\n"
		"The input filename extension can be anything, but don't use \".dat\".\n"
		"The output file has the same name as the input file but with the \n"
		"extension \".dat\".\n\n"
		);
    fprintf(stderr, "Options:\n"
        "\t-h, --help              print usage message\n"
        "\t-v, --verbose           be verbose\n"
        "\t-g, --global_kappa      use the same kappa in all layers (= kappa_sp)\n"
        "\t--nr <nr>               specify nr\n"
        "\t--nz <nz>               specify nz\n"
        "\t--nt <nt>               specify nt\n"
        "\t--nt_scale <factor>     specify scale factor for nt\n"
		"\t--ez1 <ez1>             specify z-position of bottom of cylinder (<0)\n"
		"\t--ez2 <ez2>             specify z-position of top of cylinder (>0)\n"
        "\t--alpha_so <alpha_so>   specify alpha_so\n"
        "\t--alpha_sp <alpha_sp>   specify initial alpha_sp\n"
        "\t--alpha_sr <alpha_sr>   specify alpha_sr\n"
        "\t--theta_so <theta_so>   specify theta_so\n"
        "\t--theta_sp <theta_sp>   specify initial theta_sp\n"
        "\t--theta_sr <theta_sr>   specify theta_sr\n"
        "\t--kappa_so <kappa_so>   specify kappa_so\n"
        "\t--kappa_sp <kappa_sp>   specify initial kappa_sp\n"
        "\t--kappa_sr <kappa_sr>   specify kappa_sr\n"
        "\t--kappa_outside <k_out> specify kappa_outside (mutually excl. with -g)\n"
        "\t--alpha_step <a_step>   specify initial step in alpha_sp direction\n"
        "\t--theta_step <t_step>   specify initial step in theta_sp direction\n"
        "\t--kappa_step <k_step>   specify initial step in kappa_sp direction\n"
        "\t--minalpha <minalpha>   specify minimum value of alpha_sp\n"
        "\t--maxalpha <maxalpha>   specify maximum value of alpha_sp\n"
        "\t--mintheta <mintheta>   specify minimum value of theta_sp\n"
        "\t--maxtheta <maxtheta>   specify maximum value of theta_sp\n"
        "\t--minkappa <minkappa>   specify minimum value of kappa_sp\n"
        "\t--maxkappa <maxkappa>   specify maximum value of kappa_sp\n"
        "\t--tmax <tmax>           specify total duration of experiment\n"
        "\t--fit_tol <fit_tol>     specify stopping criterion (simplex size)\n"
        "\t--itermax <itermax>     specify stopping criterion (max iterations)\n"
        "\t--outfile <outfile>     specify output file (parameters and curves)\n"
        "\t--pathfile <pathfile>   specify simplex path output file (just \n"
        "\t                        one vertex of the simplex per iteration)\n"
        );
    exit(EXIT_FAILURE);
}

/**
  \brief Make sure that the filename is not too long

   If the input string is not to big, copy it to the output string. 
   Otherwise, exit with an error message.  This function is used 
   to check that the user input filename is not too long. 

  \param [in] in Input string = the user input filename (or base filename)
  \param [out] out Output string = the input string, if not too long
 */
void check_filename(char *in, char *out)
{
	// The input and output filenames are 4 characters longer than 
	// the base filename, hence the "- 4"
	if (strlen(in) < (FILENAME_MAX - 4)) 
		strcpy(out, in);
	else 
		error("Filename length is too long");
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


/**
  \brief Create a string containing the command used to run the program. 

  It can be very useful for the program's output file to have a 
  comment line that contains the command used to run the program.
  This function assembles that command into a string from the 
  argc and argv parameters so it can be added to the output file.

  \param [in] argc Argument count (parameter passed to main)
  \param [in] argv Argument vector (parameter passed to main)
  \param [out] command String containing the command line

  \return Number of elements on the command line
***********************************************************************/
int assemble_command(int argc, char *argv[], char *command)
{   
    int i, length;
 
	strcpy(command,argv[0]);
	strcat(command, " ");
	length=strlen(argv[0]) + 1;
	for (i=1; i<argc; i++) {
		length += strlen(argv[i]) + 1;
		if (length > MAX_COMMAND_LENGTH) {
			printf("Warning: length of command too long to put in struct \n");
			strcat(command, "...");
			break;
		}
		strcat(command, argv[i]);
		strcat(command, " ");
	}

	return (i);
}

