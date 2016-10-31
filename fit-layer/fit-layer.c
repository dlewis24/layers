/**
  \file fit-layer.c

  Fits the 3-layer model to RTI data to determine alpha,  
  theta, and kappa of the SP layer.

  Usage:

  	fit-layer [options] \<input_file\>

  where \<input_file\> is the name of the input data file.
  The output file will have the same basename as the input
  file but will have the extension '.dat'. The input file
  should not have the extension '.dat'.

  To get a list of options, run fit-layer with no input 
  filename.

  The input data file consists of:
  - Header
    - Comment lines begin with '#'
    - Lines assigning parameter values look like 
      "parameter = value [optional trailing text]"
  - Two blank lines
  - Data = 2 (or more) columns with a heading:
      - First column = time (s)
      - Second column = concentration (mM)

  Notes:
  - Units are always m^2/s for dfree, nA for current, s for duration,
    and microns for distance -- no matter what appears in the optional 
    trailing text
  - If there are more than 2 columns, only the first two columns are 
    read; additional columns are skipped
  - Parameter values specified on the command line override parameter 
    values specified in the input file
  - Based on the IDL program layer.pro by Jan Hrabe, CABI, NKI
 
  \author David Lewis, CABI, NKI
  \copyright GNU Public License
  \date 2012-2013
*/

// includes
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <gsl/gsl_multimin.h>
#include "header.h"

/// Typedef for struct for passing parameters and arrays to mse function
typedef struct {
    int nt;                ///< Number of support points in time.
	int nd;                ///< Number of data points.
	int nz;                ///< Number of support points in z (rows of concentration matrix).
	int nr;                ///< Number of support points in r (columns of concentration matrix).
	int iprobe;            ///< z-index of probe location.
	int jprobe;            ///< r-index of probe location.
	int iz1;               ///< z-index of SR-SP boundary.
	int iz2;               ///< z-index of SP-SO boundary.
	int nolayer;           ///< Flag for no layer (homogenous environment).
	int opt_global_kappa;  ///< True if user specifies that kappa be the same in all layers.
    double dt;             ///< Spacing in time.
    double dr;             ///< Spacing in r (in this program, same as spacing in z).
    double sd;             ///< Source delay (time before source starts).
    double st;             ///< Duration of source.
	double alpha_so;       ///< Extracellular volume fraction in SO layer.
	double theta_so;       ///< Permeability in SO layer.
	double kappa_so;       ///< Nonspecific clearance factor in SO layer.
	double alpha_sp;       ///< Extracellular volume fraction in SP layer.
	double theta_sp;       ///< Permeability in SP layer.
	double kappa_sp;       ///< Nonspecific clearance factor in SP layer.
	double alpha_sr;       ///< Extracellular volume fraction in SR layer.
	double theta_sr;       ///< Permeability in SR layer.
	double kappa_sr;       ///< Nonspecific clearance factor in SR layer.
	double minalpha;       ///< Lower boundary for alpha_sp (add penalty if alpha_sp is out of bounds).
	double maxalpha;       ///< Upper boundary for alpha_sp (add penalty if alpha_sp is out of bounds).
	double mintheta;       ///< Lower boundary for theta_sp (add penalty if theta_sp is out of bounds).
	double maxtheta;       ///< Upper boundary for theta_sp (add penalty if theta_sp is out of bounds).
	double minkappa;       ///< Lower boundary for kappa_sp (add penalty if kappa_sp is out of bounds).
	double maxkappa;       ///< Upper boundary for kappa_sp (add penalty if kappa_sp is out of bounds).
    double dfree;          ///< Free diffusion coefficient.
    double *t;             ///< Time array for model.
    double *s;             ///< Source array.
    double *invr;          ///< Array of 1/r values.
    double *t_data;        ///< Time array for data.
    double *p_data;        ///< Probe concentration data.
    double *p;             ///< Probe concentration calculated from model.
} param_struct_type;


/**
  \brief Mean squared error function for simplex fitting. 

  Calls calc_diffusion_curve_layer_fit_layer(), 
  which calculates the diffusion curve from the layer model. 
  Calculates and returns the MSE between the model curve 
  and the data from experiment.

  Note that the data typically has about 1000 samples, but 
  the model curve generally has several thousand points. 
  In order to calculate the mean squared error, the model 
  curve is downsampled to the same number of sample points 
  as in the data. 

  \author Dave Lewis, CABI, NKI

  \param [in,out] x Vector of parameters to fit (alpha, theta, kappa of SP)
  \param [in,out] params Struct of parameters and arrays (e.g., geometry of environment, nt, time array, data array, model curve array)

  \return Mean squared error between model and data; the model data is returned in an array in the params struct
*/
double calc_mse_fit_layer(const gsl_vector *x, void *params)
{
	int i = -1;
	param_struct_type *p = (param_struct_type *) params;
	int nt = p->nt;
	int nd = p->nd;
	int p_index = -1;
	double index_scale = -1.;


	p->alpha_sp = gsl_vector_get(x, 0);
	p->theta_sp = gsl_vector_get(x, 1);
	p->kappa_sp = gsl_vector_get(x, 2);
	if (p->alpha_sp <= 0.001) p->alpha_sp = 0.001;
	if (p->theta_sp <= 0.001) p->theta_sp = 0.001;

	if (p->opt_global_kappa) {
		p->kappa_sr = p->kappa_sp;
		p->kappa_so = p->kappa_sp;
	}
   
	calc_diffusion_curve_layer_fit_layer(p->nt, p->nz, p->nr, 
		p->iprobe, p->jprobe, p->iz1, p->iz2, 
		p->nolayer, p->dt, p->dr, p->sd, p->st, 
		p->alpha_so, p->theta_so, p->kappa_so, 
		p->alpha_sp, p->theta_sp, p->kappa_sp, 
		p->alpha_sr, p->theta_sr, p->kappa_sr, 
		p->dfree, p->t, p->s, p->invr, p->p);

	double mse = 0.;

	if (nt > nd) {
		index_scale = (double) nt / (double) nd;
		for (i=1; i<nd; i++) {
			p_index = (lround) (i * index_scale);
			mse += SQR(p->p[p_index] - p->p_data[i]);
		}
		mse /= nd;
	} else {
		index_scale = (double) nd / (double) nt;
		for (i=1; i<nt; i++) {
			p_index = (lround) (i * index_scale);
			mse += SQR(p->p[i] - p->p_data[p_index]);
		}
		mse /= nt;
	}

	double penalty_factor = 10.;  // Hard-coding this 

	if (p->alpha_sp < p->minalpha) 
		mse += (p->minalpha - p->alpha_sp) * penalty_factor;
	if (p->alpha_sp > p->maxalpha) 
		mse += (p->alpha_sp - p->maxalpha) * penalty_factor;
	if (p->theta_sp < p->mintheta) 
		mse += (p->mintheta - p->theta_sp) * penalty_factor;
	if (p->theta_sp > p->maxtheta) 
		mse += (p->theta_sp - p->maxtheta) * penalty_factor;
	if (p->kappa_sp < p->minkappa) 
		mse += (p->minkappa - p->kappa_sp) * penalty_factor;
	if (p->kappa_sp > p->maxkappa) 
		mse += (p->kappa_sp - p->maxkappa) * penalty_factor;

	return mse;
}


/// Main program
int main(int argc, char *argv[])
{
	// Input parameters 
	int opt = -1;
	int opt_index = -1;
	int opt_help = FALSE;
	int opt_verbose = FALSE;
	int opt_pathfile = FALSE;
	int num_args_left = -1;
	FILE *file_ptr = NULL;
	FILE *pathfile_ptr = NULL;
	char basename[FILENAME_MAX];
	memset(basename, '\0', FILENAME_MAX);
	char infilename[FILENAME_MAX];
	memset(infilename, '\0', FILENAME_MAX);
	char outfilename[FILENAME_MAX];
	memset(outfilename, '\0', FILENAME_MAX);
	char pathfilename[FILENAME_MAX];
	memset(pathfilename, '\0', FILENAME_MAX);
	char string[MAX_LINELENGTH];
	memset(string, '\0', MAX_LINELENGTH);
	char parameter[MAX_LINELENGTH];
	memset(parameter, '\0', MAX_LINELENGTH);
	char value[MAX_LINELENGTH];
	memset(value, '\0', MAX_LINELENGTH);
	char junk_char = '\0';
	int i, j, k, l; 
	i = j = k = -1;
	int linelength = -1;
	int found_header_end = FALSE;
	int found_data_end = FALSE;

	time_t start_time = -1;
	time_t end_time = -1;
	double total_time = -1.;
	double program_version = 0.2;

	// struct for holding comments from the input parameter file 
	struct {
		int n;
		char line[MAXNUM_COMMENTLINES][MAX_LINELENGTH];
		char command[MAX_COMMAND_LENGTH];
	} comments;
	comments.n = 0;
	memset(comments.command, '\0', MAX_COMMAND_LENGTH);

	// Geometry 
	double rmax = 1000.0e-6;  // Maximum extent in r (m) 
	double zmax = 2000.0e-6;  // Maximum extent in z (m) 
	int specified_zmax = FALSE;
	double lz1 = -1.;  // Lower boundary of SP in z direction (m)
	int specified_lz1 = FALSE;  // True if user specifies lz1
	double lz2 = -1.;  // Upper boundary of SP in z direction (m)
	int specified_lz2 = FALSE;  // True if user specifies lz1
	int iz1 = -1;  // z-index of lower boundary of SP 
	int iz2 = -1;  // z-index of upper boundary of SP 
	int nolayer = FALSE; // Flag for no layer (homogenous environment)

	double ez1 = -1.;   // Location of lower edge of cylinder 
	int specified_ez1 = FALSE;
	double ez2 = -1.;   // Location of upper edge of cylinder 
	int specified_ez2 = FALSE;

	// Discretization steps in r, z, and t 
	int nr = 500;  // Number of discrete samples in r
	int nz = 1000;  // Number of discrete samples in z
	int nt = -1;  // Number of discrete samples in t
	// Default determined from von Neumann stability criterion 

	int specified_nt = FALSE;
	double nt_scale = -1.;
	int specified_nt_scale = FALSE;
	int ns = -1;
	int nds = -1;	// Number of time points in delay period before source 
	double dr = -1.;
	double dz = -1.;
	double dt = -1.;

	// Source 
	double trn = 0.35;  // Transport number of the source electrode
	double crnt = 80.0e-9;  // Iontophoretic current (A)
	double sa = -1.;  // Source amplitude
	double tmax = 150.0;  // Maximum extent in time (s) 
	double sd = 10.0;  // Delay before source begins (s)
	double st = 50.0;  // Source duration (s)
	double sr = 0.0;  // Source r-position (m) (always 0) 
	double sz = -1.;  // Source z-position (m) (should be input as 0)
	/* In previous incarnations of this program, the user could 
	   specify the source position, but in this version it is 0 
	   (source position defines the origin).  Since we're transitioning 
	   from one convention to another, the user should input sz 
	   explicitly and the value given here won't matter. */

	int isource, jsource;  // z- and r-indices of source position
	isource = jsource = -1;

	// Probe 
	double pr = 0.0;  // Probe r coordinate
	double pz = -1.;  // Probe z coordinate
	int specified_pz = FALSE;  // True if user specifies pz
	int iprobe, jprobe;  // z- and r-indices of probe position
	iprobe = jprobe = -1;

	// ECS parameters -- got defaults from paper -- p. 12 and Table 1 
	// of manuscript submitted in summer 2011 
	int opt_global_kappa = FALSE;  // True if user specifies that kappa be the same in all layers
	double alpha_so = 0.218;
	double theta_so = 0.447;
	double kappa_so = 0.007;
	double alpha_sp = 0.2;
	double theta_sp = 0.4;
	double kappa_sp = 0.01;
	double alpha_sr = 0.218;
	double theta_sr = 0.447;
	double kappa_sr = 0.007;
	double kappa_outside = 0.007;  // Used when user specifies kappa outside the SP layer (= kappa_so = kappa_sr)
	int specified_kappa_outside = FALSE;  // True when user specifies 
	                                      // kappa outside the SP layer
	// dstar = dfree * theta 
	// dstar_max is used in calculating dt from the von Neumann criterion
	double dstar_so, dstar_sp, dstar_sr, dstar_max;  
	dstar_so = dstar_sp = dstar_sr = dstar_max = -1.;

	double dfree = 1.24e-09;  // Free diffusion coefficient 

	// Arrays 
	double *t = NULL;  // Time
//	double *s = NULL;  // Source(z,r)
	double *p = NULL;  // Probe (calculated from 3-layer model)
	double *alphas = NULL;  // alpha(z,r)
//	double *invr = NULL;  // 1/r

	// Data arrays 
	double *tdata = NULL;  // Time (for data values, from input file) 
	double *pdata = NULL;  // Data (from input file)
	int	nd;  // Number of data points 

	// Parameters to send to calc_mse_fit_layer 
	param_struct_type param_struct;

	param_struct.nt = -1;
	param_struct.nd = -1;
	param_struct.nz = -1;
	param_struct.nr = -1;
	param_struct.iprobe = -1;
	param_struct.jprobe = -1;
	param_struct.iz1 = -1;
	param_struct.iz2 = -1;
	param_struct.nolayer = -1;
	param_struct.opt_global_kappa = -1;

	param_struct.dt = -1.;
	param_struct.dr = -1.;
	param_struct.sd = -1.;
	param_struct.st = -1.;
	param_struct.alpha_so = -1.;
	param_struct.theta_so = -1.;
	param_struct.kappa_so = -1.;
	param_struct.alpha_sp = -1.;
	param_struct.theta_sp = -1.;
	param_struct.kappa_sp = -1.;
	param_struct.alpha_sr = -1.;
	param_struct.theta_sr = -1.;
	param_struct.kappa_sr = -1.;
	param_struct.minalpha = -1.;
	param_struct.maxalpha = -1.;
	param_struct.mintheta = -1.;
	param_struct.maxtheta = -1.;
	param_struct.minkappa = -1.;
	param_struct.maxkappa = -1.;
	param_struct.dfree = -1.;

	param_struct.t = NULL;
	param_struct.s = NULL;
	param_struct.invr = NULL;
	param_struct.p = NULL;
	param_struct.t_data = NULL;
	param_struct.p_data = NULL;


	// Parameters for curve fitting 
	double minalpha = 0.001;  // Add penalty if alpha_sp outside range
	double maxalpha = 0.25;
	double mintheta = 0.001;  // Add penalty if theta_sp outside range
	double maxtheta = 0.75;
	double minkappa = 0.0;    // Add penalty if kappa_sp outside range
	double maxkappa = 0.1;
	double alpha_fit = -1.;  // Value of alpha_sp from fit
	double theta_fit = -1.;  // Value of theta_sp from fit
	double kappa_fit = -1.;  // Value of kappa_sp from fit
	double mse = -1.;        // Mean squared error from fit
	gsl_vector *steps = NULL;  // Step sizes for simplex
	gsl_vector *simplex = NULL;  // Simplex for minimization
	double alpha_step = 0.1;  // Initial step size for alpha_sp
	double theta_step = 0.2;  // Initial step size for theta_sp
	double kappa_step = 0.002;  // Initial step size for kappa_sp
	// The following are used by the GSL minimization algorithm
	const gsl_multimin_fminimizer_type *fit_algorithm = 
		gsl_multimin_fminimizer_nmsimplex;
	gsl_multimin_fminimizer *fit_state = NULL;
	gsl_multimin_function fit_func;
	size_t fit_iter = 0;  // Counter for iterations of minimization algo.
	int itermax = 100;  // Maximum number of iterations of min. algorithm
	int fit_status = -1;  
	double fit_size = -1.;
	double fit_tol = 1.e-4;  // Fit tolerance:  stopping criterion; 
	                         // a lower fit_tol gives a more precise 
	                         // (but not necessarily more accurate) fit


	// Get start time of program 
	start_time = time(NULL);

	/* If no arguments are given, or if there is on argument 
	   beginning with - (like -h), print the usage statement 
	    ((argc == 2) &&(argv[argc-1][0] == '-'))) 
Note: Previously it parsed the command line and then the input file. 
It used the last argument as the name of the input file. 
Now it parses the input file first, because that makes it easy for 
the command-line parameters to override the parameters in the input 
file (the desired behavior). But now it doesn't handle options that 
come after the input filename.
Currently it just checks if the first character of the last arg is -
*/
	if ((argc < 2) || (argv[argc-1][0] == '-')) 
		print_usage_fit_layer(argv[0]);


/****************************************************
 Read the parameters and the data from the input file
 ****************************************************/

	/* Get the name of the input file or the basename. Possibilities:
		User specifies	Name of input file	Name of output file
		--------------	------------------	-------------------
		<basename>     	<basename>.txt   	<basename>.dat  
		<basename>.txt 	<basename>.txt   	<basename>.dat
	   Search for the '.' in the input name. If none was found, 
	   the basename was specified; append .txt to get the 
	   input filename and append .dat to get the output filename. 
	   If the '.' was found, the input filename was specified; 
	   for the output filename, replace the suffix with ".out".  */
	check_filename(argv[argc-1], infilename);
	strcpy(outfilename, infilename);

	char *strng;
	strng = index(outfilename, '.');  // Find the '.' in the filename 
	if (strng == NULL) { 	// If no '.', then add the extensions 
		strcat(infilename, ".txt");
		strcat(outfilename, ".dat");
	} else {
		strcpy(strng, ".dat");  // Replace the suffix with ".dat" 
	}


	// Open file 
	if ((file_ptr = fopen(infilename,"r")) == NULL) {
		fprintf(stderr, "Error opening input file %s\n", infilename);
		exit(EXIT_FAILURE);
	}

	// Read input parameter file and parse the parameters 
	comments.n = 0;	// counter for the comment lines in the parameter file 
	for (i=0; i<MAXNUM_LINES; i++) {
		if (fgets(string,MAX_LINELENGTH,file_ptr) == NULL) 
			break;  // Assume that EOF was found, rather than an error 

		// If the line that was read into 'string' was not EOF 
		// and there was no error, continue; get length of 'string' 
		linelength = strlen(string);

		// If the line starts with '#', it's a comment; 
		// copy the comment, but don't parse it 
		if (string[0] == '#') {
			if (comments.n > MAXNUM_COMMENTLINES)
				fprintf(stderr, "Warning: Maximum # of comment lines "
					"exceeded.\nWill not copy more comment lines to "
					"the output file.\n");
			else
				strcpy(comments.line[comments.n], string);
			comments.n++;
		// If the line is 1 character long, stop reading header 
		} else if (linelength < 3) {
			found_header_end = TRUE;
			break;
		// If the line is too long, just skip it 
		} else if (linelength >= MAX_LINELENGTH-1) {
			fprintf(stderr, "Warning: Line %d seems to be too long\n", i+1);
			// Put a null character at end of string and feel good 
			string[MAX_LINELENGTH-1] = '\0'; 
			// Ignore the rest of the line 
			if (fscanf(file_ptr, "%*[^\n] %[\n]", &junk_char) == EOF)
				error("fscanf returned EOF");
		// Read parameters 
		} else {
			if (sscanf(string, "%s = %s", parameter, value) == EOF)
				error("scanf returned EOF");
			if (STREQ(parameter, "dfree")) {
				dfree = atof(value);
				/* Normally dfree = 1.24e-9 or so. However, in some 
				   parameter files dfree might be read as 1.24 instead 
				   (e.g. if the parameter file has 1.24e_09, or if 
				   the user just inputs 1.24). So if the value for 
				   dfree is unrealistically high, assume that it 
				   needs to be multiplied by 1e-9. */
				if (dfree > 0.01) dfree *= 1e-9;
			}
			if (STREQ(parameter, "trn")) trn = atof(value);
			if (STREQ(parameter, "current")) {
				crnt = atof(value);
				crnt *= 1e-9;	// Input in nA; convert to A 
			}
			if (STREQ(parameter, "delay")) sd = atof(value);
			if (STREQ(parameter, "duration")) st = atof(value);
			if (STREQ(parameter, "source_z")) {
				sz = atof(value);
				// sz *= 1e-6; 	// Input in microns; convert to m 
				if (! IS_ZERO(sz)) // The source location should be 0. 
					error("source_z = %f microns but should be 0 "
						"(or not specified in the output file)", sz);
			}
			if (STREQ(parameter, "probe_z")) {
				pz = atof(value);
				specified_pz = TRUE;
				pz *= 1e-6;	// Input in microns; convert to m 
			}
			if (STREQ(parameter, "probe_r")) {
				pr = atof(value);
				pr *= 1e-6;	// Input in microns; convert to m 
			}
			if (STREQ(parameter, "nolayer")) nolayer = atoi(value);
			if (STREQ(parameter, "lz1")) {
				lz1 = atof(value);
				specified_lz1 = TRUE;
				lz1 *= 1e-6;	// Input in microns; convert to m 
			}
			if (STREQ(parameter, "lz2")) {
				lz2 = atof(value);
				specified_lz2 = TRUE;
				lz2 *= 1e-6;	// Input in microns; convert to m 
			}
			if (STREQ(parameter, "ez1")) {
				ez1 = atof(value);
				specified_ez1 = TRUE;
				ez1 *= 1e-6;    // Input in microns; convert to m 
			}
			if (STREQ(parameter, "ez2")) {
				ez2 = atof(value);
				specified_ez2 = TRUE;
				ez2 *= 1e-6;    // Input in microns; convert to m 
			}
			if (STREQ(parameter, "alpha_so")) alpha_so = atof(value);
			if (STREQ(parameter, "alpha_sr")) alpha_sr = atof(value);
			if (STREQ(parameter, "theta_so")) theta_so = atof(value);
			if (STREQ(parameter, "theta_sr")) theta_sr = atof(value);
			if (STREQ(parameter, "kappa_so")) kappa_so = atof(value);
			if (STREQ(parameter, "kappa_sr")) kappa_sr = atof(value);
/* Removing kappa_outside option from input file
   - could lead to confusion
   - not very useful
   - I probably have not used it in the input file
			if (STREQ(parameter, "kappa_outside")) {
				kappa_outside = atof(value);
				specified_kappa_outside = TRUE;
			}
*/
			if (STREQ(parameter, "nt")) {
				nt = atoi(value);
				specified_nt = TRUE;
			}
			if (STREQ(parameter, "nt_scale")) {
				nt_scale = atof(value);
				specified_nt_scale = TRUE;
			}
			if (STREQ(parameter, "nr")) nr = atoi(value);
			if (STREQ(parameter, "nz")) nz = atoi(value);
			if (STREQ(parameter, "rmax")) {
				rmax = atof(value);
				rmax *= 1e-6;	// Input in microns; convert to m 
			}
			if (STREQ(parameter, "zmax")) {
				zmax = atof(value);
				zmax *= 1e-6;	// Input in microns; convert to m 
			}
			if (STREQ(parameter, "tmax")) tmax = atof(value);
		}
	}

	if (!found_header_end)
		error("Did not find blank line after header");
else
if (opt_verbose)
printf("Found end of header\n");

	// The next line should also be blank 
	if (fgets(string,MAX_LINELENGTH,file_ptr) == NULL) 
		error("EOF (or error) reached before reading data");

	linelength = strlen(string);
	if (linelength > 2)
		error("The line after the header has %d characters "
			"(should be 1 or 2)", linelength);

	// The next line should be the header for the data 
	if (fgets(string,MAX_LINELENGTH,file_ptr) == NULL) 
		error("EOF (or error) reached before reading data");


	// Read data 
	// Make the data arrays long enough to hold MAXNUM_LINES values 
	tdata = create_array(MAXNUM_LINES, "tdata");
	pdata = create_array(MAXNUM_LINES, "pdata");
	if (opt_verbose)
		printf("Reading data from file\n");
	for (j=0; j<MAXNUM_LINES; j++) 
		if (fscanf(file_ptr, "%lf%lf%*[^\n] %[\n]",
			&tdata[j], &pdata[j], &junk_char) == EOF) 
		{
			found_data_end = TRUE;
			break; // Stop reading at EOF 
		}

	if (!found_data_end)
		error("Read maximum number of lines in output file (%d)\n"
			"but did not reach end of file", j);
else
if (opt_verbose)
printf("Found end of data\n");
 
	nd = j; // nd = Number of data points 

if (opt_verbose)
printf("The number of data points is %d\n", nd);

	// Close the file 
	fclose(file_ptr);


/******************************
 Parse the command line options 
 ******************************/
	static struct option long_opts[] = {
		{"help", no_argument, NULL, 'h'},
		{"verbose", no_argument, NULL, 'v'},
		{"global_kappa", no_argument, NULL, 'g'},
		{"nr", required_argument, NULL, 0},
		{"nz", required_argument, NULL, 0},
		{"nt", required_argument, NULL, 0},
		{"nt_scale", required_argument, NULL, 0},
		{"ez1", required_argument, NULL, 0},
		{"ez2", required_argument, NULL, 0},
		{"alpha_so", required_argument, NULL, 0},
		{"alpha_sp", required_argument, NULL, 0},
		{"alpha_sr", required_argument, NULL, 0},
		{"theta_so", required_argument, NULL, 0},
		{"theta_sp", required_argument, NULL, 0},
		{"theta_sr", required_argument, NULL, 0},
		{"kappa_so", required_argument, NULL, 0},
		{"kappa_sp", required_argument, NULL, 0},
		{"kappa_sr", required_argument, NULL, 0},
		{"kappa_outside", required_argument, NULL, 0},
		{"alpha_step", required_argument, NULL, 0},
		{"theta_step", required_argument, NULL, 0},
		{"kappa_step", required_argument, NULL, 0},
		{"minalpha", required_argument, NULL, 0},
		{"maxalpha", required_argument, NULL, 0},
		{"mintheta", required_argument, NULL, 0},
		{"maxtheta", required_argument, NULL, 0},
		{"minkappa", required_argument, NULL, 0},
		{"maxkappa", required_argument, NULL, 0},
		{"tmax", required_argument, NULL, 0},
		{"fit_tol", required_argument, NULL, 0},
		{"itermax", required_argument, NULL, 0},
		{"outfile", required_argument, NULL, 0},
		{"pathfile", required_argument, NULL, 0},
		{NULL, no_argument, NULL, 0}
	};

	while ((opt = getopt_long(argc, argv, "hvg", 
	                          long_opts, &opt_index)) != -1) {
		switch (opt) {

		case 'h':   // user needs help -- print usage 
			opt_help = TRUE;
			break;

		case 'v':   // user wants verbose output 
			opt_verbose = TRUE;
			break;

		case 'g':   // user wants kappa to be the same in all regions 
			opt_global_kappa = TRUE;
			break;

		case 0:   // long option has no corresponding short option 
			if (STREQ("nr", long_opts[opt_index].name)) {
				nr = atoi(optarg);
			} else if (STREQ("nz", long_opts[opt_index].name)) {
				nz = atoi(optarg);
			} else if (STREQ("nt", long_opts[opt_index].name)) {
				nt = atoi(optarg);
				specified_nt = TRUE;
			} else if (STREQ("nt_scale", long_opts[opt_index].name)) {
				nt_scale = atof(optarg);
				specified_nt_scale = TRUE;
			} else if (STREQ("ez1", long_opts[opt_index].name)) {
				ez1 = atof(optarg);
				specified_ez1 = TRUE;
				ez1 *= 1e-6;    // Input in microns; convert to m 
			} else if (STREQ("ez2", long_opts[opt_index].name)) {
				ez2 = atof(optarg);
				specified_ez2 = TRUE;
				ez2 *= 1e-6;    // Input in microns; convert to m 
			} else if (STREQ("alpha_so", long_opts[opt_index].name)) {
				alpha_so = atof(optarg);
			} else if (STREQ("alpha_sp", long_opts[opt_index].name)) {
				alpha_sp = atof(optarg);
			} else if (STREQ("alpha_sr", long_opts[opt_index].name)) {
				alpha_sr = atof(optarg);
			} else if (STREQ("theta_so", long_opts[opt_index].name)) {
				theta_so = atof(optarg);
			} else if (STREQ("theta_sp", long_opts[opt_index].name)) {
				theta_sp = atof(optarg);
			} else if (STREQ("theta_sr", long_opts[opt_index].name)) {
				theta_sr = atof(optarg);
			} else if (STREQ("kappa_so", long_opts[opt_index].name)) {
				kappa_so = atof(optarg);
			} else if (STREQ("kappa_sp", long_opts[opt_index].name)) {
				kappa_sp = atof(optarg);
			} else if (STREQ("kappa_sr", long_opts[opt_index].name)) {
				kappa_sr = atof(optarg);
			} else if (STREQ("kappa_outside", long_opts[opt_index].name)) {
				kappa_outside = atof(optarg);
				specified_kappa_outside = TRUE;
			} else if (STREQ("alpha_step", long_opts[opt_index].name)) {
				alpha_step = atof(optarg);
			} else if (STREQ("theta_step", long_opts[opt_index].name)) {
				theta_step = atof(optarg);
			} else if (STREQ("kappa_step", long_opts[opt_index].name)) {
				kappa_step = atof(optarg);
			} else if (STREQ("minalpha", long_opts[opt_index].name)) {
				minalpha = atof(optarg);
			} else if (STREQ("maxalpha", long_opts[opt_index].name)) {
				maxalpha = atof(optarg);
			} else if (STREQ("mintheta", long_opts[opt_index].name)) {
				mintheta = atof(optarg);
			} else if (STREQ("maxtheta", long_opts[opt_index].name)) {
				maxtheta = atof(optarg);
			} else if (STREQ("minkappa", long_opts[opt_index].name)) {
				minkappa = atof(optarg);
			} else if (STREQ("maxkappa", long_opts[opt_index].name)) {
				maxkappa = atof(optarg);
			} else if (STREQ("tmax", long_opts[opt_index].name)) {
				tmax = atof(optarg);
			} else if (STREQ("fit_tol", long_opts[opt_index].name)) {
				fit_tol = atof(optarg);
			} else if (STREQ("itermax", long_opts[opt_index].name)) {
				itermax = atoi(optarg);
			} else if (STREQ("outfile", long_opts[opt_index].name)) {
				check_filename(optarg, outfilename);
			} else if (STREQ("pathfile", long_opts[opt_index].name)) {
				check_filename(optarg, pathfilename);
				opt_pathfile = TRUE;
			}
			break;

		case ':':   // missing argument for an option 
			error("Missing argument for a command-line option '-%c'", 
				optopt);
			break;

		case '?':   // unknown option 
			error("Unknown command line option '-%c'", optopt);
			break;

		default:
			error("Problem in parsing command line options");
			break;
		}
	}

	// See if any options are left 
	num_args_left = argc - optind;
	if ((num_args_left != 1) || opt_help)
		print_usage_fit_layer(argv[0]);


	if (opt_verbose) {
		printf("The name of the input file is %s\n", infilename);
		printf("The name of the output file will be %s\n", outfilename);
		if (opt_pathfile) 
			printf("The name of the simplex path file will be %s\n", 
				pathfilename);
	}

	// Check for conflicts
	if (STREQ(infilename, outfilename)) 
		error("The input and output filenames cannot be the same.");
	if (STREQ(infilename, pathfilename)) 
		error("The input and simplex path filenames cannot be the same.");
	if (STREQ(outfilename, pathfilename)) 
		error("The output and simplex path filenames cannot be the same.");

    if (specified_ez1 && !specified_ez2)
        error("You specified ez1 but did not specify ez2");
    if (specified_ez2 && !specified_ez1)
        error("You specified ez2 but did not specify ez1");
    if (specified_ez1 && specified_zmax)
        error("You specified ez1 and ez2, so you should not specify zmax");

	if (specified_ez1) {
		if (ez1 > 0) error("Bottom of cylinder ez1 = %f > 0\n", ez1);
		if (ez2 < 0) error("Top of cylinder ez2 = %f < 0\n", ez2);
		if (ez1 > lz1) error("Bottom of cylinder ez1 = %f > lz1 = %f\n",
		                      ez1, lz1);
		if (ez2 < lz2) error("Top of cylinder ez2 = %f < lz2 = %f\n",
		                      ez2, lz2);
	}


	// Default values of parameters that depend on zmax 
	if (specified_pz == FALSE) {
		pz = 120.0e-6; // m 
if (opt_verbose)
		printf("Warning: probe location set to default value of %g m = "
			"%f microns\n (relative to source)", 
			pz, 1e6*pz);
	}

	if (specified_lz1 == FALSE) {
		lz1 = - 50.0e-6 / 2.0; 
if (opt_verbose)
		printf("Warning: lz1 set to default value of %g m = "
			"%f microns\n (relative to source)", 
			lz1, 1e6*lz1);
	}

	if (specified_lz2 == FALSE) {
		lz2 = lz1 + 50.0e-6; 
if (opt_verbose)
		printf("Warning: lz2 set to default value of %g m = "
			"%f microns\n (relative to source)", 
			lz2, 1e6*lz2);
	}

	// Set kappa in SR and SO to kappa_outside 
	if (specified_kappa_outside) {
		kappa_sr = kappa_outside;
		kappa_so = kappa_outside;
		if (opt_global_kappa) {
			error("You're fitting for global kappa but specified "
				"kappa_outside.\nWhen you fit for global kappa, "
				"kappa_sr and kappa_so are set = kappa_sp.");
		}
	}

	/* If nolayer option, the SR parameters are used for the whole
	   environment; the SO and SP parameters are not used. However,
	   their values are still listed in the output file. Set their
	   values equal to the SR values here so that the correct values
	   are printed to the output file. */
    if (nolayer) {
        alpha_so = alpha_sr;
        alpha_sp = alpha_sr;
        theta_so = theta_sr;
        theta_sp = theta_sr;
        kappa_so = kappa_sr;
        kappa_sp = kappa_sr;
        if (opt_verbose)
            printf("\nNOTE: nolayer option given; the diffusion "
                    "parameters of \nthe homogeneous environment "
                    "are set to the SR values\n");
    }

	// If using a global value for kappa, set kappa_sr=kappa_so=kappa_sp 
	if (opt_global_kappa) {
		kappa_sr = kappa_sp;
		kappa_so = kappa_sp;
		if (opt_verbose)
			printf("NOTE: kappa will be the same in all layers (-g)\n"
					"kappa_sr and kappa_so set to kappa_sp\n");
	}

/*
 * Change coordinates. In the coordinate system used in the
 * input file, source_z is always 0, so lz1, lz2, and pz are
 * measured relative to the source. This was done because in
 * the experiments distances are measured relative to the
 * source.
 *
 * In the coordinate system used in the model, the bottom of
 * the cylinder is at z=0 and the top of the cylinder is at
 * z=zmax. The length of the cylinder is therefore zmax.
 * So we shift the z-coordinates here.
 *
 * By default the z-coordinates are shifted such that the
 * SP layer is centered in the volume, ie the midplane of SP is
 * at z=zmax/2. However, the user might instead want to specify
 * the location of the ends of the cylinder relative to the
 * source (ez1 and ez2), which requires a different z-shift.
 */
	double coord_shift;

	if (specified_ez1) {
		zmax = ez2 + (-ez1);    // Calculate the cylinder length 
		coord_shift = (-ez1);
	} else {
		coord_shift = (zmax - (lz1+lz2))/2.;
	}
	sz = coord_shift;
	pz += coord_shift;
	lz1 += coord_shift;
	lz2 += coord_shift;


	// Discretization intervals in r, z 
	dr = rmax / nr;
	dz = zmax / nz;
	if (fabs(dr - dz) > 1.0e-15) {
		dr = dz;
		rmax = dr * nr;
	}

	sz = round(sz / dz) * dz;
	pz = round(pz / dz) * dz;
	pr = round(pr / dr) * dr;

	// Layer geometry 
	iz1 = (long) (lz1 / dz);
	lz1 = iz1 * dz + dz / 2.0;

	iz2 = (long) (lz2 / dz);
	lz2 = iz2 * dz + dz / 2.0;


	// D* 
	dstar_so = theta_so * dfree;
	dstar_sp = theta_sp * dfree;
	dstar_sr = theta_sr * dfree;
	dstar_max = MAX(dstar_so, dstar_sp);
	dstar_max = MAX(dstar_max, dstar_sr);


	// Check if layer thickness is numerically reasonable 
	if ( ((iz2 - iz1) < 2) && (nolayer == 0) ) 
		error("Layer has too few discrete steps to continue.");

	// Calculate time step from nt or from von Neumann criterion 
	if (specified_nt == TRUE)
		dt = tmax / nt;
	else
		dt = 0.9 * dr*dr / (6.0 * dstar_max);

	// Scale nt if the user specified to do so -- but do it via dt 
	if (specified_nt_scale == TRUE) {
		if (IS_ZERO(nt_scale)) error("nt_scale = 0");
		if (nt_scale < 0.) error("nt_scale < 0");
		dt /= nt_scale;
	}

	// Calculate tmax and st as multiples of dt 
	nt = lround(tmax / dt);
	tmax = dt * nt;
	ns = lround(st / dt);
	st = dt * ns;
	nds = lround(sd / dt);
	sd = dt * nds;

	if (sd >= tmax) 
		error("Source delay (%f) should be < tmax (%f)", sd, tmax);
	if (st >= tmax) 
		error("Source duration (%f) should be < tmax (%f)", st, tmax);
	if (sd+st >= tmax) 
		error("Source delay (%f) + duration (%f) should be < tmax (%f)", 
			sd, st, tmax);
		

	// Calculate source amplitude in mol/s from current and transport number
	sa = crnt * trn / FARADAY; // source strength in mol/s 
	                           // (not a concentration) 


	// Assemble string with command that user input 
	i = assemble_command(argc, argv, comments.command);
	if (opt_verbose)
		printf("\nIn main(): The command used was\n\t%s\n(%d words)\n\n", 
			comments.command, i);


	// Put start time in string 
	strncpy(string, ctime(&start_time), MAX_LINELENGTH-1);
 
	// Feedback to user to check adjusted values 
	if (opt_verbose) {
		printf("Output from fit-layer.c, version %.1f:\n", program_version);
		printf("# Note that the z-values (sz, pz, lz1, and lz2) ");
		printf("have been shifted \n# by %f microns ", 1.0e6*coord_shift);
		if (specified_ez1) {
			printf("to have the volume go from z=0 to z=zmax.\n");
		} else {
			printf("to center the SP layer in the volume.\n");
		}
		printf("nr x nz = %d x %d\n", nr, nz);
		printf("rmax x zmax = %f x %f microns\n", 1.0e6 * rmax, 1.0e6 * zmax);
		printf("dr x dz = %f x %f microns\n", 1.0e6 * dr, 1.0e6 * dz);
		printf("(sr, sz) = (%f, %f) microns\n", 1.0e6 * sr, 1.0e6 * sz);
		printf("(pr, pz) = (%f, %f) microns\n", 1.0e6 * pr, 1.0e6 * pz);
		printf("Electrode distance = %f microns\n", 
			1.0e6 * sqrt(SQR(pr-sr) + SQR(pz-sz)));
		printf("(iz1, iz2) = (%d, %d)\n", iz1, iz2);
		printf("(lz1, lz2) = (%f, %f) microns\n", 1.0e6 * lz1, 1.0e6 * lz2);
		printf("Layer thickness = %f microns\n", 1.0e6 * (lz2 - lz1));
		printf("Layer discrete steps = %d\n", iz2 - iz1);
		printf("Nolayer flag = %d\n" , nolayer);
		printf("dfree = %g m^2/s\n", dfree);
		printf("alpha_so = %.4f, theta_so = %.4f, "
			"lambda_so = %.4f, kappa_so = %.6f\n", 
			alpha_so, theta_so, 1.0/sqrt(theta_so), kappa_so);
		printf("Starting alpha_sp = %.4f, theta_sp = %.4f, "
			"lambda_sp = %.4f, kappa_sp = %.6f\n", 
			alpha_sp, theta_sp, 1.0/sqrt(theta_sp), kappa_sp);
		printf("Starting alpha_step = %.4f, theta_step = %.4f\n", 
			alpha_step, theta_step);
		printf("Constraints: minalpha = %.8f, maxalpha = %.8f\n", 
			minalpha, maxalpha);
		printf("Constraints: mintheta = %.8f, maxtheta = %.8f\n", 
			mintheta, maxtheta);
		printf("Constraints: minkappa = %.8f, maxkappa = %.8f\n", 
			minkappa, maxkappa);
		printf("Stopping criteria: simplex size < %g or # iterations = %d\n", 
			fit_tol, itermax);
		printf("alpha_sr = %.4f, theta_sr = %.4f, "
			"lambda_sr = %.4f, kappa_sr = %.6f\n", 
			alpha_sr, theta_sr, 1.0/sqrt(theta_sr), kappa_sr);
		if (opt_global_kappa) 
			printf("NOTE: kappa_sr and kappa_so set to kappa_sp (-g)\n");
		printf("nt = %d\n", nt);
		printf("tmax = %f s\n", tmax);
		printf("dt = %f ms\n", 1.0e3 * dt);
		printf("von Neumann dt / (dr^2/(6*dstar)) = %f\n", 
			dt * 6.0 * dstar_max / (dr*dr));
		printf("ns = %d\n", ns);
		printf("source delay sd = %f s\n", sd);
		printf("source duration st = %f s\n", st);
		printf("Current = %g nA\n", 1.0e9 * crnt);
		printf("Transport number = %f\n", trn);
		printf("Start time = %s", string);
	}


	// Open output data file 
	if ((file_ptr = fopen(outfilename,"w")) == NULL) {
		fprintf(stderr, "Error opening output file %s\n", outfilename);
		exit(EXIT_FAILURE);
	}

	// Write out command line to output file 
	fprintf(file_ptr, "# Fit-layer Output File\n");
	fprintf(file_ptr, "# ~~~~~~~~~~~~~~~~~~~~~\n");
	fprintf(file_ptr, "# Command used to run program:\n");
	fprintf(file_ptr, "# %s\n", comments.command);

	// Write comments from input file to output file 
	if (comments.n > 0) {
		fprintf(file_ptr, "# --------------------------------------\n");
		fprintf(file_ptr, "# Comments from input parameter file:\n");
		for (j=0; j<comments.n; j++)
			fprintf(file_ptr, "%s", comments.line[j]);
		fprintf(file_ptr, "# --------------------------------------\n");
	}


	// Write parameters to output file 
	fprintf(file_ptr, "# Output from fit-layer.c, version %.1f:\n", program_version);
	fprintf(file_ptr, "# Note that the z-values (sz, pz, lz1, and lz2) ");
	fprintf(file_ptr, "have been shifted \n# by %f microns ", 1.0e6*coord_shift);
	if (specified_ez1) {
		fprintf(file_ptr, "to have the volume go from z=0 to z=zmax.\n");
	} else {
		fprintf(file_ptr, "to center the SP layer in the volume.\n");
	}
	fprintf(file_ptr, "# nr x nz = %d x %d\n", nr, nz);
	fprintf(file_ptr, "# rmax x zmax = %f x %f microns\n", 1.0e6 * rmax, 1.0e6 * zmax);
	fprintf(file_ptr, "# dr x dz = %f x %f microns\n", 1.0e6 * dr, 1.0e6 * dz);
	fprintf(file_ptr, "# (sr, sz) = (%f, %f) microns\n", 1.0e6 * sr, 1.0e6 * sz);
	fprintf(file_ptr, "# (pr, pz) = (%f, %f) microns\n", 1.0e6 * pr, 1.0e6 * pz);
	fprintf(file_ptr, "# Electrode distance = %f microns\n", 
		1.0e6 * sqrt(SQR(pr-sr) + SQR(pz-sz)));
	fprintf(file_ptr, "# (iz1, iz2) = (%d, %d)\n", iz1, iz2);
	fprintf(file_ptr, "# (lz1, lz2) = (%f, %f) microns\n", 1.0e6 * lz1, 1.0e6 * lz2);
	fprintf(file_ptr, "# Layer thickness = %f microns\n", 1.0e6 * (lz2 - lz1));
	fprintf(file_ptr, "# Layer discrete steps = %d\n", iz2 - iz1);
	fprintf(file_ptr, "# Nolayer flag = %d\n" , nolayer);
	fprintf(file_ptr, "# dfree = %g m^2/s\n", dfree);
	fprintf(file_ptr, "# alpha_so = %.4f, theta_so = %.4f, "
			"lambda_so = %.4f, kappa_so = %.6f\n", 
			alpha_so, theta_so, 1.0/sqrt(theta_so), kappa_so);
	fprintf(file_ptr, "# Starting alpha_sp = %.4f, theta_sp = %.4f, "
			"lambda_sp = %.4f, kappa_sp = %.6f\n", 
			alpha_sp, theta_sp, 1.0/sqrt(theta_sp), kappa_sp);
	fprintf(file_ptr, "# Starting alpha_step = %.4f, theta_step = %.4f\n", 
			alpha_step, theta_step);
	fprintf(file_ptr, "# Constraints: minalpha = %.8f, maxalpha = %.8f\n", 
			minalpha, maxalpha);
	fprintf(file_ptr, "# Constraints: mintheta = %.8f, maxtheta = %.8f\n", 
			mintheta, maxtheta);
	fprintf(file_ptr, "# Constraints: minkappa = %.8f, maxkappa = %.8f\n", 
			minkappa, maxkappa);
	fprintf(file_ptr, "# Stopping criteria: simplex size < %g or # iterations = %d\n", 
			fit_tol, itermax);
	fprintf(file_ptr, "# alpha_sr = %.4f, theta_sr = %.4f, "
			"lambda_sr = %.4f, kappa_sr = %.6f\n", 
			alpha_sr, theta_sr, 1.0/sqrt(theta_sr), kappa_sr);
	if (opt_global_kappa) 
		fprintf(file_ptr, "# NOTE: kappa_sr and kappa_so set to kappa_sp (-g)\n");
	fprintf(file_ptr, "# nt = %d\n", nt);
	fprintf(file_ptr, "# tmax = %f s\n", tmax);
	fprintf(file_ptr, "# dt = %f ms\n", 1.0e3 * dt);
	fprintf(file_ptr, "# von Neumann dt / (dr^2/(6*dstar)) = %f\n", 
		dt * 6.0 * dstar_max / (dr*dr));
	fprintf(file_ptr, "# ns = %d\n", ns);
	fprintf(file_ptr, "# Source delay sd = %f s\n", sd);
	fprintf(file_ptr, "# Source duration st = %f s\n", st);
	fprintf(file_ptr, "# Current = %g nA\n", 1.0e9 * crnt);
	fprintf(file_ptr, "# Transport number = %f\n", trn);
	fprintf(file_ptr, "# Start time = %s", string); // ctime() added the \n 

	// Close the file 
	fclose(file_ptr);

	// Array of alpha values -- use 1D index to access elements of 
	// the 2D array, so a[i*(nr+1)+j] = a[i][j] (use macro INDEX(i,j) 
    alphas = create_array(nz*(nr+1), "alphas array");
	for (j=0; j<nr+1; j++) {
		for (i=0;     i<iz1+1; i++) alphas[INDEX(i,j)] = alpha_sr;
		for (i=iz1+1; i<iz2+1; i++) alphas[INDEX(i,j)] = alpha_sp;
		for (i=iz2+1; i<nz;    i++) alphas[INDEX(i,j)] = alpha_so;
	}

	// Array of 1/r values, except it is 0 for r=0 
	// - in layer.pro it's a 2D array with all column identical, 
	//   but it doesn't have to be; I'm doing 1D (i only) 
    param_struct.invr = create_array(nr+1, "param invr array");
	param_struct.invr[0] = 1.0 / dr;	
	param_struct.invr[1] = 0.0;
	for (j=2; j<nr+1; j++) {
		param_struct.invr[j] = 1.0 / ((j-1.)*dr);
	}


	// Source 
    param_struct.s = create_array(nz*(nr+1), "param s array");
	isource = lround(sz/dz);   	// index to z position of source 
	jsource = 1+lround(sr/dr);	// index to r position of source 
	param_struct.s[INDEX(isource,jsource)] 
		= (1.0 / alphas[INDEX(isource,jsource)]) * 
			sa * dt * 4.0 / (PI * SQR(dr) * dz);


	// Time and probe arrays 
	t = create_array(nt, "time");
	p = create_array(nt, "p");
	for (k=0; k<nt; k++) 
		t[k] = dt * k;

	iprobe = lround(pz/dz);   	// index to z position of probe 
	jprobe = 1+lround(pr/dr);	// index to r position of probe 



	// Fit parameters
	if (opt_verbose)
		printf("About to fit parameters\n");


	// Fit the model p[] to the data pdata[] 

	if (opt_verbose) {
		printf("\nSimplex fitting -- vertex changes:\n");
		printf("Iter\talpha_fit\ttheta_fit\tkappa_fit\tmse      \tfit size\n");
		printf("%d\t%f\t%f\t%f\n", 
			0, alpha_sp, theta_sp, kappa_sp);
	}

	if (opt_pathfile) {
		if ((pathfile_ptr = fopen(pathfilename,"w")) == NULL) {
			fprintf(stderr, "Error opening simplex path file %s\n", 
				pathfilename);
			exit(EXIT_FAILURE);
		}
		fprintf(pathfile_ptr, "Simplex fitting -- vertex changes:\n");
		fprintf(pathfile_ptr, 
			"Iter\talpha_fit\ttheta_fit\tkappa_fit\tmse      \tfit size\n");
		fprintf(pathfile_ptr, "%d\t%f\t%f\t%f\n", 
			0, alpha_sp, theta_sp, kappa_sp);
	}

	param_struct.nt = nt;
	param_struct.nd = nd;
	param_struct.nz = nz;
	param_struct.nr = nr;
	param_struct.iprobe = iprobe;
	param_struct.jprobe = jprobe;
	param_struct.iz1 = iz1;
	param_struct.iz2 = iz2;
	param_struct.nolayer = nolayer;
	param_struct.opt_global_kappa = opt_global_kappa;

	param_struct.dt = dt;
	param_struct.dr = dr;
	param_struct.sd = sd;
	param_struct.st = st;
	param_struct.alpha_so = alpha_so;
	param_struct.theta_so = theta_so;
	param_struct.kappa_so = kappa_so;
	param_struct.alpha_sp = alpha_sp;
	param_struct.theta_sp = theta_sp;
	param_struct.kappa_sp = kappa_sp;
	param_struct.alpha_sr = alpha_sr;
	param_struct.theta_sr = theta_sr;
	param_struct.kappa_sr = kappa_sr;
	param_struct.minalpha = minalpha;
	param_struct.maxalpha = maxalpha;
	param_struct.mintheta = mintheta;
	param_struct.maxtheta = maxtheta;
	param_struct.minkappa = minkappa;
	param_struct.maxkappa = maxkappa;
	param_struct.dfree = dfree;

    param_struct.t = create_array(nt, "param t array");

    param_struct.p = create_array(nt, "param p array");
    param_struct.t_data = create_array(nd, "param t_data array");
    param_struct.p_data = create_array(nd, "param p_data array");
	for (k=0; k<nt; k++) {
    	param_struct.t[k] = t[k];
    	param_struct.p[k] = p[k];
	}

	// Fill t_data and p_data 
	for (k=0; k<nd; k++) {
    	param_struct.t_data[k] = tdata[k];
    	param_struct.p_data[k] = pdata[k];
	}

/*****************************************
 Fit the model to determine the parameters
 *****************************************/
	// Initialize the simplex 
	simplex = gsl_vector_alloc(3);
	gsl_vector_set(simplex, 0, alpha_sp);
	gsl_vector_set(simplex, 1, theta_sp);
	gsl_vector_set(simplex, 2, kappa_sp);

	// Initialize step sizes 
	steps = gsl_vector_alloc(3);
	gsl_vector_set(steps, 0, alpha_step);
	gsl_vector_set(steps, 1, theta_step);
	gsl_vector_set(steps, 2, kappa_step);

	// Set up minimization method 
	fit_func.n = 3;  // 3 parameters to fit 
	fit_func.f = calc_mse_fit_layer;  // function to minimize 
	fit_func.params = &param_struct;  // extra parameters to function 

	fit_state = gsl_multimin_fminimizer_alloc(fit_algorithm, 3);
	gsl_multimin_fminimizer_set(fit_state, &fit_func, simplex, steps);

	// Run minimization
	do {
		fit_iter++;
		fit_status = gsl_multimin_fminimizer_iterate(fit_state);

		if (fit_status) break;

		fit_size = gsl_multimin_fminimizer_size(fit_state);
		fit_status = gsl_multimin_test_size(fit_size, fit_tol);

		if (opt_verbose)
			if (fit_status == GSL_SUCCESS) printf("Finished fit\n");

		alpha_fit = gsl_vector_get(fit_state->x, 0);
		theta_fit = gsl_vector_get(fit_state->x, 1);
		kappa_fit = gsl_vector_get(fit_state->x, 2);
		mse = fit_state->fval;

		if (opt_verbose)
			printf("%d\t%f\t%f\t%f\t%g\t%g\n", 
				(int) fit_iter, alpha_fit, theta_fit, kappa_fit, mse, fit_size);

		if (opt_pathfile) 
			fprintf(pathfile_ptr, "%d\t%f\t%f\t%f\t%g\t%g\n", 
				(int) fit_iter, alpha_fit, theta_fit, kappa_fit, mse, fit_size);

	} while (fit_status == GSL_CONTINUE && fit_iter < itermax);

	if (fit_status != GSL_SUCCESS) {
		printf("Warning: failed to converge, status = %d, "
			"# iterations = %zd\n", fit_status, fit_iter);
		if (opt_pathfile) 
			fprintf(pathfile_ptr, "Warning: failed to converge, "
			"status = %d, # iterations = %zd\n", fit_status, fit_iter);
	}

	if (opt_pathfile) 
		fclose(pathfile_ptr);

	// Output results
	double lambda_fit = 1./sqrt(theta_fit);
	if (opt_verbose) {
		printf("Fitted alpha = %f\n", alpha_fit);
		printf("Fitted theta = %f  (lambda = %f)\n", 
			theta_fit, lambda_fit);
		if (opt_global_kappa) 
			printf("Fitted kappa = %f s^-1 (in all layers)\n", kappa_fit);
		else
			printf("Fitted kappa = %f s^-1\n", kappa_fit);
	}


	// Get end time of program 
	end_time = time(NULL);
	strncpy(string, ctime(&end_time), MAX_LINELENGTH-1);
if (opt_verbose)
	printf("End time = %s", string);

	total_time = difftime(end_time, start_time);
if (opt_verbose)
	printf("Total time = %d seconds = %f minutes = %f hours\n", 
		(int) round(total_time), total_time/60., total_time/3600.);


	// Write results to output file. Two concentration columns,
	// the concentration from the layer model and also the data.
	// Each concentration column has its own time column preceding it. 
	// (The data was input with their own time values. 
	//  The concentration curve calculated from the model typically 
	//  has >> 1000 points and is downsampled to 1000. 
	//  The time columns are close, especially for high resolution.) 
	if ((file_ptr = fopen(outfilename,"a")) == NULL) {
		fprintf(stderr, "Error opening output file %s\n", outfilename);
		exit(EXIT_FAILURE);
	}
 
	fprintf(file_ptr, "# End time = %s", string); // ctime() added the \n 
	fprintf(file_ptr, "# Total time = %d seconds = %f minutes = %f hours\n", 
		(int) round(total_time), total_time/60., total_time/3600.);
	fprintf(file_ptr, "# --------------------------------------\n");
	fprintf(file_ptr, "# Results of fitting:\n");
	fprintf(file_ptr, "# Number of iterations = %d\n", (int)fit_iter);
	fprintf(file_ptr, "# Fitted alpha = %f\n", alpha_fit);
	fprintf(file_ptr, "# Fitted theta = %f  (lambda = %f)\n", 
		theta_fit, lambda_fit);
	if (opt_global_kappa) 
		fprintf(file_ptr, "# Fitted kappa = %f s^-1 (in all layers)\n", kappa_fit);
	else
		fprintf(file_ptr, "# Fitted kappa = %f s^-1\n", kappa_fit);
	fprintf(file_ptr, "# Final mean squared error = %g\n", mse);
	fprintf(file_ptr, "# Final simplex size = %g\n", fit_size);

	fprintf(file_ptr, "# Solution: alpha_sp\ttheta_sp\tlambda_sp\tkappa_sp"
	                  "\t     MSE\tsimplex size\t# iter.\tTime (s)"
	                  "\tTime (m)\tTime (h) \n");
	fprintf(file_ptr, "# Solution: %f\t%f\t%f\t%f\t%f\t%g  \t%7d\t%8d\t%f\t%f\n",
		alpha_fit, theta_fit, lambda_fit, kappa_fit, 
		mse, fit_size, (int) fit_iter, 
		(int) round(total_time), total_time/60., total_time/3600.);
	fprintf(file_ptr, "# --------------------------------------\n");
	fprintf(file_ptr, "# Probe concentration data:\n");
	fprintf(file_ptr, "#   time      \t  c (model) \t  t (data) "
		"\t    c (data) \n");


	// Print concentration arrays to output file 
	// If there are more than 1000 points in the concentration values 
	// from the model (normally nt >> 1000), downsample to 1000 points
	if (nt > 1000) 
		for (i=0; i<1000; i++) {
			k = (i * nt) / 1000;
			l = (i * nd) / 1000;
			fprintf(file_ptr, "%#12.8g\t%#12.8g\t%#12.8g\t%#12.8g\n", param_struct.t[k], param_struct.p[k], 
				param_struct.t_data[l], param_struct.p_data[l]);
		}
	else
		for (i=0; i<nt; i++) {
			fprintf(file_ptr, "%#12.8g\t%#12.8g\t%#12.8g\t%#12.8g\n", param_struct.t[i], param_struct.p[i], 
				param_struct.t_data[i], param_struct.p_data[i]);
		}
	fprintf(file_ptr, "\n\n\n");

	// Close the file 
	fclose(file_ptr);



	// Deallocate arrays 
	free(t);
	free(p);
//	free(s);
	free(tdata);
	free(pdata);
	free(alphas);
//	free(invr);

	free(param_struct.t);
	free(param_struct.s);
	free(param_struct.invr);
	free(param_struct.t_data);
	free(param_struct.p_data);
	free(param_struct.p);

	gsl_vector_free(simplex);
	gsl_vector_free(steps);
	gsl_multimin_fminimizer_free(fit_state);

if (opt_verbose)
	printf("All done\n");

	return(EXIT_SUCCESS);
}
