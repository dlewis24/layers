/**
  \file 3layer/3layer.c

  Calculates the extracellular concentration of a substance 
  diffusing in an environment consisting of 3 adjacent homogeneous, 
  isotropic layers.

  Usage:

    3layer [options] \<input_file\>

  where \<input_file\> is the name of the input parameter file.
  The output file will have the same basename as the input file
  but will have the extension '.dat'. The input file should not
  have the extension '.dat'.

  To get a list of options, run 3layer with no input filename.

 
  The input parameter file consists of:
  - Comment lines beginning with '#'
  - Lines assigning parameter values that look like 
    "parameter = value [<trailing text>]"
    - For example,
      duration = 50 s (source duration)
    - The program will assign the value 50 to the parameter 'duration'

  Notes:
  - Units are always m^2/s for dfree, nA for current, s for duration,
    and microns for distance -- no matter what appears in the optional
    trailing text
  - Parameter values specified on the command line override parameter
    values specified in the input file
  - Based on the IDL program layer.pro by Jan Hrabe, CABI, NKI

  \author David Lewis, CABI, NKI
  \copyright GNU Public License
  \date 2012-2013 
*/

// includes
#include "header.h"


/// Main program 
int main(int argc, char *argv[])
{
	/* Input parameters */
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
	int i, j, k; 
	i = j = k = -1;
	int linelength = -1;

	time_t start_time = -1;
	time_t end_time = -1;
	double total_time = -1.;
	double program_version = 0.2;

	/* struct for holding comments from the input parameter file */
	struct {
		int n;
		char line[MAXNUM_COMMENTLINES][MAX_LINELENGTH];
		char command[MAX_COMMAND_LENGTH];
	} comments;
	comments.n = 0;
	memset(comments.command, '\0', MAX_COMMAND_LENGTH);

	/* Geometry */
	double rmax = 1000.0e-6;	/* m */
	double zmax = 2000.0e-6;	/* m */
	int specified_zmax = FALSE;
	double lz1 = -1.;	/* Default depends on zmax */
	int specified_lz1 = FALSE;
	double lz2 = -1.;	/* Default depends on zmax */
	int specified_lz2 = FALSE;
	int iz1 = -1;
	int iz2 = -1;
	int nolayer = FALSE; /* layer flag (simplifies numerical 
	                        calculation of a homogeneous case) */
	double ez1 = -1.;	/* Location of lower edge of cylinder */
	int specified_ez1 = FALSE;
	double ez2 = -1.;	/* Location of upper edge of cylinder */
	int specified_ez2 = FALSE;

	/* Discretization steps in r, z, and t */
	int nr = 500;
	int nz = 1000;
	int nt = -1;	/* Default determined from von Neumann 
					   stability criterion */
	int specified_nt = FALSE;
	double nt_scale = -1.;
	int specified_nt_scale = FALSE;
	int ns = -1;
	int nds = -1;	/* Number of points in delay period before source */
	double dr = -1.;
	double dz = -1.;
	double dt = -1.;

	/* Source */
	double trn = 0.35;
	double crnt = 80.0e-9;	/* A */
	double samplitude = -1.;
	double tmax = 150.0; 	/* s */
	double sdelay = 10.0;	/* s */
	double sduration = 50.0;	/* s */
	double sr = 0.0;	/* Source defines origin in this program */
	double sz = -1.;	/* Default depends on zmax */
	int isource, jsource;
	isource = jsource = -1;

	/* Probe */
	double pr = 0.0;	/* m */
	double pz = -1.;	/* Default depends on sz */
	int specified_pz = FALSE;
	int iprobe, jprobe;
	iprobe = jprobe = -1;

	/* ECS parameters -- got defaults from paper -- p. 12 and Table 1 
	   of manuscript submitted in summer 2011 */
	int opt_global_kappa = FALSE;
	double alpha_so = 0.218;
	double theta_so = 0.447;
	double kappa_so = 0.0;
	double alpha_sp = 0.2;
	double theta_sp = 0.4;
	double kappa_sp = 0.0;
	double alpha_sr = 0.218;
	double theta_sr = 0.447;
	double kappa_sr = 0.0;
	double kappa_outside = 0.0;
	int specified_kappa_outside = FALSE;
	double dstar_so, dstar_sp, dstar_sr, dstar_max;
	dstar_so = dstar_sp = dstar_sr = dstar_max = -1.;

	double dfree = 1.24e-09;

	/* Arrays */
	double *t = NULL;	/* time */
	double *s = NULL;	/* source */
	double *p = NULL;	/* probe (calculated) */
	double *alphas = NULL;	/* alpha */
	double *invr = NULL;	/* 1/r */

	/* Parameters for output files */
	int opt_output_conc_image = FALSE;
	double image_spacing = 1.;
	char imagebasename[FILENAME_MAX-32];
	memset(imagebasename, '\0', FILENAME_MAX-32);

	/* Parameters for additional sources */
	int nsource;
	char *token;
	char additional_sources_string[ADDITIONAL_SOURCES_STRING_LENGTH];
	memset(additional_sources_string, '\0', ADDITIONAL_SOURCES_STRING_LENGTH);
	more_sources_struct_type more_sources;
	more_sources.n = 0;
	more_sources.source = NULL;
	source_struct_type new_source;

	/* Parameters to send to calc_mse_rti */
	double spdist = -1.;

	mse_rti_params_struct_type mse_rti_params;

	mse_rti_params.nt = -1;
	mse_rti_params.spdist = -1.;
	mse_rti_params.samplitude = -1.;
	mse_rti_params.sdelay = -1.;
	mse_rti_params.sduration = -1.;
	mse_rti_params.kappa = -1.;
	mse_rti_params.dfree = -1.;
	mse_rti_params.alpha = -1.;
	mse_rti_params.theta = -1.;

	mse_rti_params.t = NULL;
	mse_rti_params.p_model = NULL;
	mse_rti_params.p_theory = NULL;


	/* Parameters for curve fitting */
	double alpha_start = 0.2;  /* Sarting value for theoretical fit */
	double theta_start = 0.4;  /* Sarting value for theoretical fit */
	double alpha_fit = -1.; 
	double theta_fit = -1.; 
	double mse = -1.;
	gsl_vector *steps = NULL;
	gsl_vector *simplex = NULL;
	double alpha_step = 0.1;
	double theta_step = 0.2;
	const gsl_multimin_fminimizer_type *fit_algorithm = 
		gsl_multimin_fminimizer_nmsimplex;
	gsl_multimin_fminimizer *fit_state = NULL;
	gsl_multimin_function fit_func;
	size_t fit_iter = 0;
	int itermax = 100;
	int fit_status = -1;
	double fit_size = -1.;
	double fit_tol = 1.e-4;


	/* Get start time of program */
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
		print_usage(argv[0]);

/*
 * Read the parameters and the data from the input file
 */

	/* Get the names of the input and output files from the last 
	   argument on the command line. (The output filename determined 
	   here can be overridden with a command line argument.) */
	get_io_filenames(argv[argc-1], ".par", ".dat", infilename, outfilename);

	/* Open file */
	if ((file_ptr = fopen(infilename,"r")) == NULL) {
		fprintf(stderr, "Error opening input file %s\n", infilename);
		exit(EXIT_FAILURE);
	}

	/* Read input parameter file and parse the parameters */
	comments.n = 0;	/* counter for the comment lines in the parameter file */
	for (i=0; i<MAXNUM_LINES; i++) {
		if (fgets(string,MAX_LINELENGTH,file_ptr) == NULL) 
			break;  /* Assume that EOF was found, rather than an error */

		/* If the line that was read into 'string' was not EOF
		   and there was no error, continue; get length of 'string' */
		linelength = strlen(string);

		/* If the line starts with '#', it's a comment;
		   copy the comment, but don't parse it */
		if (string[0] == '#') {
			if (comments.n > MAXNUM_COMMENTLINES)
				fprintf(stderr, "Warning: Maximum # of comment lines "
					"exceeded.\nWill not copy more comment lines to "
					"the output file.\n");
			else
				strcpy(comments.line[comments.n], string);
			comments.n++;
		/* If the line is 1 character long, stop reading header */
		} else if (linelength < 3) {
			break;
		/* If the line is too long, just skip it */
		} else if (linelength >= MAX_LINELENGTH-1) {
			fprintf(stderr, "Warning: Line %d seems to be too long\n", i+1);
			/* Put a null character at end of string and feel good */
			string[MAX_LINELENGTH-1] = '\0'; 
			/* Ignore the rest of the line */
			if (fscanf(file_ptr, "%*[^\n] %[\n]", &junk_char) == EOF)
				error("fscanf returned EOF");
		/* Read parameters */
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
				crnt *= 1e-9;	/* Input in nA; convert to A */
			}
			if (STREQ(parameter, "delay")) sdelay = atof(value);
			if (STREQ(parameter, "duration")) sduration = atof(value);
			if (STREQ(parameter, "source_z")) {
				sz = atof(value);
				/* sz *= 1e-6; */	/* Input in microns; convert to m */
				if (! IS_ZERO(sz)) /* The source location should be 0. */
					error("source_z = %f microns but should be 0 "
						"(or not specified in the output file)", sz);
			}
			if (STREQ(parameter, "probe_z")) {
				pz = atof(value);
				specified_pz = TRUE;
				pz *= 1e-6;	/* Input in microns; convert to m */
			}
			if (STREQ(parameter, "probe_r")) {
				pr = atof(value);
				pr *= 1e-6;	/* Input in microns; convert to m */
			}
			if (STREQ(parameter, "nolayer")) nolayer = atoi(value);
			if (STREQ(parameter, "lz1")) {
				lz1 = atof(value);
				specified_lz1 = TRUE;
				lz1 *= 1e-6;	/* Input in microns; convert to m */
			}
			if (STREQ(parameter, "lz2")) {
				lz2 = atof(value);
				specified_lz2 = TRUE;
				lz2 *= 1e-6;	/* Input in microns; convert to m */
			}
			if (STREQ(parameter, "ez1")) {
				ez1 = atof(value);
				specified_ez1 = TRUE;
				ez1 *= 1e-6;	/* Input in microns; convert to m */
			}
			if (STREQ(parameter, "ez2")) {
				ez2 = atof(value);
				specified_ez2 = TRUE;
				ez2 *= 1e-6;	/* Input in microns; convert to m */
			}
			if (STREQ(parameter, "alpha_so")) alpha_so = atof(value);
			if (STREQ(parameter, "alpha_sp")) alpha_sp = atof(value);
			if (STREQ(parameter, "alpha_sr")) alpha_sr = atof(value);
			if (STREQ(parameter, "theta_so")) theta_so = atof(value);
			if (STREQ(parameter, "theta_sp")) theta_sp = atof(value);
			if (STREQ(parameter, "theta_sr")) theta_sr = atof(value);
			if (STREQ(parameter, "kappa_so")) kappa_so = atof(value);
			if (STREQ(parameter, "kappa_sp")) kappa_sp = atof(value);
			if (STREQ(parameter, "kappa_sr")) kappa_sr = atof(value);
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
				rmax *= 1e-6;	/* Input in microns; convert to m */
			}
			if (STREQ(parameter, "zmax")) {
				zmax = atof(value);
				zmax *= 1e-6;	/* Input in microns; convert to m */
				specified_zmax = TRUE;
			}
			if (STREQ(parameter, "tmax")) tmax = atof(value);
		}
	}


	/* Close the file */
	fclose(file_ptr);


/* 
 * Parse the command line options 
 */
	static struct option long_opts[] = {
		{"help", no_argument, NULL, 'h'},
		{"verbose", no_argument, NULL, 'v'},
		{"global_kappa", no_argument, NULL, 'g'},
		{"nr", required_argument, NULL, 0},
		{"nz", required_argument, NULL, 0},
		{"nt", required_argument, NULL, 0},
		{"nt_scale", required_argument, NULL, 0},
		{"probe_z", required_argument, NULL, 0},
		{"probe_r", required_argument, NULL, 0},
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
		{"alpha_start", required_argument, NULL, 0},
		{"theta_start", required_argument, NULL, 0},
		{"alpha_step", required_argument, NULL, 0},
		{"theta_step", required_argument, NULL, 0},
		{"tmax", required_argument, NULL, 0},
		{"fit_tol", required_argument, NULL, 0},
		{"itermax", required_argument, NULL, 0},
		{"outfile", required_argument, NULL, 0},
		{"pathfile", required_argument, NULL, 0},
		{"images", required_argument, NULL, 0},
		{"image_spacing", required_argument, NULL, 0},
		{"additional_sources", required_argument, NULL, 0},
		{NULL, no_argument, NULL, 0}
	};

	while ((opt = getopt_long(argc, argv, "hvg", 
	                          long_opts, &opt_index)) != -1) {
		switch (opt) {

		case 'h':   /* user needs help -- print usage */
			opt_help = TRUE;
			break;

		case 'v':   /* user wants verbose output */
			opt_verbose = TRUE;
			break;

		case 'g':   /* user wants kappa to be the same in all regions */
			opt_global_kappa = TRUE;
			break;

		case 0:   /* long option has no corresponding short option */
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
			} else if (STREQ("probe_z", long_opts[opt_index].name)) {
				pz = atof(optarg);
				specified_pz = TRUE;
				pz *= 1e-6;	/* Input in microns; convert to m */
			} else if (STREQ("probe_r", long_opts[opt_index].name)) {
				pr = atof(optarg);
				pr *= 1e-6;	/* Input in microns; convert to m */
			} else if (STREQ("ez1", long_opts[opt_index].name)) {
				ez1 = atof(optarg);
				specified_ez1 = TRUE;
				ez1 *= 1e-6;	/* Input in microns; convert to m */
			} else if (STREQ("ez2", long_opts[opt_index].name)) {
				ez2 = atof(optarg);
				specified_ez2 = TRUE;
				ez2 *= 1e-6;	/* Input in microns; convert to m */
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
			} else if (STREQ("alpha_start", long_opts[opt_index].name)) {
				alpha_start = atof(optarg);
			} else if (STREQ("theta_start", long_opts[opt_index].name)) {
				theta_start = atof(optarg);
			} else if (STREQ("alpha_step", long_opts[opt_index].name)) {
				alpha_step = atof(optarg);
			} else if (STREQ("theta_step", long_opts[opt_index].name)) {
				theta_step = atof(optarg);
			} else if (STREQ("tmax", long_opts[opt_index].name)) {
				tmax = atof(optarg);
			} else if (STREQ("fit_tol", long_opts[opt_index].name)) {
				fit_tol = atof(optarg);
			} else if (STREQ("itermax", long_opts[opt_index].name)) {
				itermax = atoi(optarg);
			} else if (STREQ("outfile", long_opts[opt_index].name)) {
				get_filename(optarg, outfilename);
			} else if (STREQ("pathfile", long_opts[opt_index].name)) {
				get_filename(optarg, pathfilename);
				opt_pathfile = TRUE;
			} else if (STREQ("images", long_opts[opt_index].name)) {
				get_filename(optarg, imagebasename);
				opt_output_conc_image = TRUE;
			} else if (STREQ("image_spacing", long_opts[opt_index].name)) {
				image_spacing = atof(optarg);
			} else if (STREQ("additional_sources", long_opts[opt_index].name)) {
				if (strlen (optarg) < ADDITIONAL_SOURCES_STRING_LENGTH)
					strcpy(additional_sources_string, optarg);
				else
					error("additional_sources_string is too long");

				token = strtok(additional_sources_string, " ,");
				more_sources.n = atoi(token);

				if (more_sources.n > 0) {
					more_sources.source = 
						(source_struct_type *)
							malloc( sizeof(source_struct_type) * more_sources.n );
					if (more_sources.source == NULL)
						error("Cannot allocate memory for more_sources");
					/* Read parameters for each additional source */
					for (nsource = 0; nsource < more_sources.n; nsource++) {
						more_sources.source[nsource].sz = 
							read_source_parameter("sz", nsource) * 1e-6; 
						more_sources.source[nsource].sr = 
							read_source_parameter("sr", nsource) * 1e-6; 
						more_sources.source[nsource].crnt = 
							read_source_parameter("crnt", nsource) * 1e-9; 
					}
				}
			}
			break;

		case ':':   /* missing argument for an option */
			error("Missing argument for a command-line option '-%c'", 
				optopt);
			break;

		case '?':   /* unknown option */
			error("Unknown command line option '-%c'", optopt);
			break;

		default:
			error("Problem in parsing command line options");
			break;
		}
	}


	/* See if any options are left */
	num_args_left = argc - optind;
	if ((num_args_left != 1) || opt_help)
		print_usage(argv[0]);


	if (opt_verbose) {
		printf("The name of the input file is %s\n", infilename);
		printf("The name of the output file will be %s\n", outfilename);
		if (opt_pathfile) 
			printf("The name of the simplex path file will be %s\n", 
				pathfilename);
	}

	/* Check for conflicts */
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


	/* Default values of parameters that depend on zmax */
	if (specified_pz == FALSE) {
		pz = 120.0e-6; /* m */
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

/* Set kappa in SR and SO to kappa_outside */
	if (specified_kappa_outside) {
		kappa_sr = kappa_outside;
		kappa_so = kappa_outside;
		if (opt_global_kappa) {
			error("You've specified both global kappa and kappa_outside.\n"
		"The global kappa option sets kappa_so and kappa_sr to kappa_sp.\n"
		"When you specify kappa_outside, kappa_so and kappa_sr are set \n"
		"to that value.");
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

/* If using a global value for kappa, set kappa_sr=kappa_so=kappa_sp */
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
		zmax = ez2 + (-ez1);	/* Calculate the cylinder length */
		coord_shift = (-ez1);
	} else {
		coord_shift = (zmax - (lz1+lz2))/2.;
	}
	sz = coord_shift;
	pz += coord_shift;
	lz1 += coord_shift;
	lz2 += coord_shift;

	/* Discretization intervals in r, z */
	dr = rmax / nr;
	dz = zmax / nz;
	if (fabs(dr - dz) > 1.0e-15) {
		dr = dz;
		rmax = dr * nr;
	}

	sz = round(sz / dz) * dz;
	pz = round(pz / dz) * dz;
	pr = round(pr / dr) * dr;

	/* Layer geometry */
	iz1 = (int) round((lz1 / dz));
	lz1 = iz1 * dz + dz / 2.0;

	iz2 = (int) round((lz2 / dz));
	lz2 = iz2 * dz + dz / 2.0;


	/* D* */
	dstar_so = theta_so * dfree;
	dstar_sp = theta_sp * dfree;
	dstar_sr = theta_sr * dfree;
	dstar_max = MAX(dstar_so, dstar_sp);
	dstar_max = MAX(dstar_max, dstar_sr);


	/* Check if layer thickness is numerically reasonable */
	if ( ((iz2 - iz1) < 2) && (nolayer == 0) ) 
		error("Layer has too few discrete steps to continue.");

	/* Calculate time step from nt or from von Neumann criterion */
	if (specified_nt == TRUE)
		dt = tmax / nt;
	else
		dt = 0.9 * dr*dr / (6.0 * dstar_max);

	/* Scale nt if the user specified to do so -- but do it via dt */
	if (specified_nt_scale == TRUE) {
		if (IS_ZERO(nt_scale)) error("nt_scale = 0");
		if (nt_scale < 0.) error("nt_scale < 0");
		dt /= nt_scale;
	}

	/* Round off tmax, sduration, and sdelay to multiples of dt */
	nt = lround(tmax / dt);
	tmax = dt * nt;
	ns = lround(sduration / dt);
	sduration = dt * ns;
	nds = lround(sdelay / dt);
	sdelay = dt * nds;

	if (sdelay >= tmax) 
		error("Source delay (%f) should be < tmax (%f)", sdelay, tmax);
	if (sduration >= tmax) 
		error("Source duration (%f) should be < tmax (%f)", sduration, tmax);
	if (sdelay+sduration >= tmax) 
		error("Source delay (%f) + duration (%f) should be < tmax (%f)", 
			sdelay, sduration, tmax);
		

	/* Calculate source amplitude in mol/s from current and transport number*/
	samplitude = crnt * trn / FARADAY; /* source strength in mol/s 
	                                      (not a concentration) */


	/* Assemble string with command that user input */
	i = assemble_command(argc, argv, comments.command);
	if (opt_verbose)
		printf("\nIn main(): The command used was\n\t%s\n(%d words)\n\n", 
			comments.command, i);


	/* Put start time in string */
	strncpy(string, ctime(&start_time), MAX_LINELENGTH-1);
 
	/* Feedback to user to check adjusted values */
	if (opt_verbose) {
		printf("Output from 3layer.c, version %.1f:\n", program_version);
		printf("Note that the z-values (sz, pz, lz1, and lz2) have been shifted \n");
		printf("by %f microns to center the SP layer in the volume.\n", 
			1.0e6*coord_shift);
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
		printf("alpha_sp = %.4f, theta_sp = %.4f, "
			"lambda_sp = %.4f, kappa_sp = %.6f\n",
			alpha_sp, theta_sp, 1.0/sqrt(theta_sp), kappa_sp);
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
		printf("source delay sdelay = %f s\n", sdelay);
		printf("source duration sduration = %f s\n", sduration);
		printf("Current = %g nA\n", 1.0e9 * crnt);
		printf("Transport number = %f\n", trn);
		printf("Start time = %s", string);
	}


	/* Open output data file */
	if ((file_ptr = fopen(outfilename,"w")) == NULL) {
		fprintf(stderr, "Error opening output file %s\n", outfilename);
		exit(EXIT_FAILURE);
	}

	/* Write out command line to output file */
	fprintf(file_ptr, "# 3layer Output File\n");
	fprintf(file_ptr, "# ~~~~~~~~~~~~~~~~~~\n");
	fprintf(file_ptr, "# Command used to run program:\n");
	fprintf(file_ptr, "# %s\n", comments.command);

	/* Write comments from input file to output file */
	if (comments.n > 0) {
		fprintf(file_ptr, "# --------------------------------------\n");
		fprintf(file_ptr, "# Comments from input parameter file:\n");
		for (j=0; j<comments.n; j++)
			fprintf(file_ptr, "%s", comments.line[j]);
		fprintf(file_ptr, "# --------------------------------------\n");
	}


	/* Write parameters to output file */
	fprintf(file_ptr, "# Output from 3layer.c, version %.1f:\n", program_version);
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
	fprintf(file_ptr, "# alpha_sp = %.4f, theta_sp = %.4f, "
		"lambda_sp = %.4f, kappa_sp = %.6f\n",
		alpha_sp, theta_sp, 1.0/sqrt(theta_sp), kappa_sp);
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
	fprintf(file_ptr, "# Source delay sdelay = %f s\n", sdelay);
	fprintf(file_ptr, "# Source duration sduration = %f s\n", sduration);
	fprintf(file_ptr, "# Current = %g nA\n", 1.0e9 * crnt);
	fprintf(file_ptr, "# Transport number = %f\n", trn);
	if (more_sources.n > 0) {
		fprintf(file_ptr, "# Number of extra sources = %d\n", more_sources.n);
		for (nsource = 0; nsource < more_sources.n; nsource++) 
			fprintf(file_ptr, "# Additional source #%d: \n" 
				"#\tsz = %lf microns, sr = %lf microns, crnt = %lf nA\n",
			nsource+1,
			1.0e6 * (more_sources.source[nsource].sz + coord_shift),
			1.0e6 * more_sources.source[nsource].sr,
			1.0e9 * more_sources.source[nsource].crnt);
	}
	fprintf(file_ptr, "# Start time = %s", string); /* ctime() added the \n */

	/* Close the file */
	fclose(file_ptr);

	/* Array of alpha values -- use 1D index to access elements of 
	   the 2D array, so a[i*(nr+1)+j] = a[i][j] (use macro INDEX(i,j) */
    alphas = create_array(nz*(nr+1), "alphas array");
	for (j=0; j<nr+1; j++) {
		for (i=0;     i<iz1+1; i++) alphas[INDEX(i,j)] = alpha_sr;
		for (i=iz1+1; i<iz2+1; i++) alphas[INDEX(i,j)] = alpha_sp;
		for (i=iz2+1; i<nz;    i++) alphas[INDEX(i,j)] = alpha_so;
	}

	/* Array of 1/r values, except it is 0 for r=0 */
	/* - in layer.pro it's a 2D array with all column identical, 
	     but it doesn't have to be; I'm doing 1D (i only) */
    invr = create_array(nr+1, "param invr array");
	invr[0] = 1.0 / dr;	
	invr[1] = 0.0;
	for (j=2; j<nr+1; j++) {
		invr[j] = 1.0 / ((j-1.)*dr);
	}


	/* Source */
    s = create_array(nz*(nr+1), "param s array");
	isource = lround(sz/dz);   	/* index to z position of source */
	jsource = 1+lround(sr/dr);	/* index to r position of source */
	s[INDEX(isource,jsource)] 
		= (1.0 / alphas[INDEX(isource,jsource)]) * 
			samplitude * dt * 4.0 / (PI * SQR(dr) * dz);

	/* If there are additional sources, add them to the s array */
	if (more_sources.n > 0) {
		for (nsource = 0; nsource < more_sources.n; nsource++) {
			new_source.sz = more_sources.source[nsource].sz;
			new_source.sr = more_sources.source[nsource].sr;
			new_source.crnt = more_sources.source[nsource].crnt;

			new_source.sz += coord_shift;	/* Convert to shifted coordinates */

			isource = lround((new_source.sz)/dz);
			jsource = 1+lround(new_source.sr/dr);
			samplitude = new_source.crnt * trn / FARADAY;

			if (isource < 0)
				error("adding additional source %d; isource = %d < 0", 
					nsource, isource);
			if (isource > nz-1)
				error("adding additional source %d; isource = %d > nz-1",
					nsource, isource);
			if (jsource < 0)
				error("adding additional source %d; jsource = %d < 0", 
					nsource, jsource);
			if (jsource > nr)
				error("adding additional source %d; jsource = %d > nr",
					nsource, jsource);

			s[INDEX(isource,jsource)] 
				+= (1.0 / alphas[INDEX(isource,jsource)]) * 
					samplitude * dt * 4.0 / (PI * SQR(dr) * dz);
		}
	}


	/* Time and probe arrays */
	t = create_array(nt, "time");
	p = create_array(nt, "p");
	for (k=0; k<nt; k++) 
		t[k] = dt * k;

	iprobe = lround(pz/dz);   	/* index to z position of probe */
	jprobe = 1+lround(pr/dr);	/* index to r position of probe */



	/* Calculate p[], the diffusion curve at the probe  */
	if (opt_verbose)
		printf("About to calculate diffusion curve\n");

	/* If the output concentration images will not be saved, 
	   set the image spacing time to -1 */
	if (! opt_output_conc_image)
		image_spacing = -1.;


	/* Calculate the concentration as a function of time and space; 
	   return the probe concentration as a function of time */
	calc_diffusion_curve_layer(nt, nz, nr, iprobe, jprobe, 
		iz1, iz2, nolayer, dt, dr, sdelay, sduration, 
		alpha_so, theta_so, kappa_so, 
		alpha_sp, theta_sp, kappa_sp, 
		alpha_sr, theta_sr, kappa_sr, 
		dfree, t, s, invr, 
		imagebasename, image_spacing, 
		p);


	/* Fit the traditional model (p_theory[]) to the concentration 
	   calculated with the multilayer model (p[]) to get the 
	   apparent parameters and characteristic curves. Surprisingly, 
	   you can often generate a characteristic curve with the 
	   traditional model that fits the calculated curve (from the 
	   multilayer model) very well.
	   This step is not really necessary -- the apparent parameters 
	   aren't physically meaningful -- but sometimes it's useful. */

	if (opt_verbose) {
		printf("\nFitting for apparent parameters/characteristic curve:\n");
		printf("Iter\talpha_fit\ttheta_fit\tmse      \tfit size\n");
		printf("%d\t%f\t%f\n", 
			0, alpha_start, theta_start);
	}

	if (opt_pathfile) {
		if ((pathfile_ptr = fopen(pathfilename,"w")) == NULL) {
			fprintf(stderr, "Error opening simplex path file %s\n", 
				pathfilename);
			exit(EXIT_FAILURE);
		}
		fprintf(pathfile_ptr, "\nFitting for apparent parameters/characteristic curve:\n");
		fprintf(pathfile_ptr, 
			"Iter\talpha_fit\ttheta_fit\tmse      \tfit size\n");
		fprintf(pathfile_ptr, "%d\t%f\t%f\n", 
			0, alpha_start, theta_start);
	}


	/* Fill struct with values for passing to mse function */
	spdist = sqrt(SQR(pr-sr) + SQR(pz-sz));
	mse_rti_params.nt = nt;
	mse_rti_params.spdist = spdist;
	mse_rti_params.samplitude = samplitude;
	mse_rti_params.sdelay = sdelay;
	mse_rti_params.sduration = sduration;
	mse_rti_params.kappa = 0.;	/* not using this currently */
	mse_rti_params.dfree = dfree;

    mse_rti_params.t = create_array(nt, "param t array");
    mse_rti_params.p_model = create_array(nt, "param p_model array");
    mse_rti_params.p_theory = create_array(nt, "param p_theory array");
	for (k=0; k<nt; k++) {
    	mse_rti_params.t[k] = t[k];
    	mse_rti_params.p_model[k] = p[k];
	}


	/* Initialize the simplex */
	simplex = gsl_vector_alloc(2);
	gsl_vector_set(simplex, 0, alpha_start);
	gsl_vector_set(simplex, 1, theta_start);

	/* Initialize step sizes */
	steps = gsl_vector_alloc(2);
	gsl_vector_set(steps, 0, alpha_step);
	gsl_vector_set(steps, 1, theta_step);

	/* Set up minimization method */
	fit_func.n = 2;  /* 2 parameters to fit */
	fit_func.f = calc_mse_rti;  /* function to minimize */
	fit_func.params = &mse_rti_params;  /* extra parameters to function */

	fit_state = gsl_multimin_fminimizer_alloc(fit_algorithm, 2);
	gsl_multimin_fminimizer_set(fit_state, &fit_func, simplex, steps);

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
		mse = fit_state->fval;

		if (opt_verbose)
			printf("%d\t%f\t%f\t%g\t%g\n", 
				(int) fit_iter, alpha_fit, theta_fit, mse, fit_size);

		if (opt_pathfile) 
			fprintf(pathfile_ptr, "%d\t%f\t%f\t%g\t%g\n", 
				(int) fit_iter, alpha_fit, theta_fit, mse, fit_size);

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

	double lambda_fit = 1./sqrt(theta_fit);
	if (opt_verbose) {
		printf("Fitted alpha = %f\n", alpha_fit);
		printf("Fitted theta = %f  (lambda = %f)\n", 
			theta_fit, lambda_fit);
	}


	/* Get end time of program */
	end_time = time(NULL);
	strncpy(string, ctime(&end_time), MAX_LINELENGTH-1);
if (opt_verbose)
	printf("End time = %s", string);

	total_time = difftime(end_time, start_time);
if (opt_verbose)
	printf("Total time = %d seconds = %f minutes = %f hours\n", 
		(int) round(total_time), total_time/60., total_time/3600.);


	/* Write concentration curves to output file. Two conc. columns
	   -- the data from the layer model and also the apparent curve */
	if ((file_ptr = fopen(outfilename,"a")) == NULL) {
		fprintf(stderr, "Error opening output file %s\n", outfilename);
		exit(EXIT_FAILURE);
	}
 
	fprintf(file_ptr, "# End time = %s", string); /* ctime() added the \n */
	fprintf(file_ptr, "# Total time = %d seconds = %f minutes = %f hours\n", 
		(int) round(total_time), total_time/60., total_time/3600.);
	fprintf(file_ptr, "# --------------------------------------\n");
	fprintf(file_ptr, "# Fit for characteristic curve:\n");
	fprintf(file_ptr, "# Number of iterations = %d\n", (int)fit_iter);
	fprintf(file_ptr, "# Fitted apparent alpha = %f\n", alpha_fit);
	fprintf(file_ptr, "# Fitted apparent theta = %f  (lambda = %f)\n", 
		theta_fit, lambda_fit);
	fprintf(file_ptr, "# Final mean squared error = %g\n", mse);
	fprintf(file_ptr, "# Final simplex size = %g\n", fit_size);
	fprintf(file_ptr, "# Solution: apparent alpha\tapparent theta\t"
	                  "apparent lambda\t     MSE\tsimplex size\t# iter."
	                  "\tTime (s)\tTime (m)\tTime (h) \n");
	fprintf(file_ptr, "# Solution: %f\t%f\t%f\t%f\t%g \t%7d\t%8d\t%f\t%f\n",
		alpha_fit, theta_fit, lambda_fit, mse, fit_size, (int) fit_iter, 
		(int) round(total_time), total_time/60., total_time/3600.);
	fprintf(file_ptr, "# --------------------------------------\n");
	fprintf(file_ptr, "# Probe concentration data:\n");
	fprintf(file_ptr, "#   time      \t  c (3-layer model) \t  c (characteristic curve) \n");


	/* Print concentration arrays to output file */
	if (nt > 1000) 
		for (i=0; i<1000; i++) {
			k = (i * nt) / 1000;
			fprintf(file_ptr, "%#12.8g\t%#12.8g\t%#12.8g\n", 
					t[k], p[k], mse_rti_params.p_theory[k]);
		}
	else
		for (i=0; i<nt; i++) {
			fprintf(file_ptr, "%#12.8g\t%#12.8g\t%#12.8g\n", 
					t[i], p[i], mse_rti_params.p_theory[i]);
		}
	fprintf(file_ptr, "\n");

	/* Close the file */
	fclose(file_ptr);



	/* Deallocate arrays */
	free(t);
	free(p);
	free(s);
	free(alphas);
	free(invr);

	free(mse_rti_params.t);
	free(mse_rti_params.p_model);
	free(mse_rti_params.p_theory);

	gsl_vector_free(simplex);
	gsl_vector_free(steps);
	gsl_multimin_fminimizer_free(fit_state);

if (opt_verbose)
	printf("All done\n");

	return(EXIT_SUCCESS);
}
