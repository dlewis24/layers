/** 
  \file fit-layer/header.h

  Header file for program fit-layer.

  \author David Lewis, CABI, NKI
  \copyright GNU Public License
  \date 2012-2013

 */

// Constants

/// Maximum number of lines in input file 
#define MAXNUM_LINES 10000

/// Maximum length of any line in input file
#define MAX_LINELENGTH 100

/// Maximum number of comment lines of input file to copy to output file
#define MAXNUM_COMMENTLINES 1000

/// Maximum number of characters of command to copy to output file
#define MAX_COMMAND_LENGTH 1000

/// FALSE assigned to 0
#define FALSE 0

/// TRUE assigned to 1
#define TRUE 1

/// Sets constant PI to M_PI if that is defined; otherwise, computes PI as \f$ 4 \: tan^{-1}(1) \f$
#ifdef M_PI
# define PI M_PI
#else
# define PI (4*atan(1.0))
#endif

/// Faraday constant in C/mol 
#define FARADAY 96485.3399

/// For comparing doubles with 0:  Sets SMALLNUM to DBL_EPSILON or to GSL_DBL_EPSILON if they're defined; otherwise set to another small number
#ifdef DBL_EPSILON
# define SMALLNUM DBL_EPSILON
#elif defined GSL_DBL_EPSILON
# define SMALLNUM GSL_DBL_EPSILON
#else
# define SMALLNUM 2.2204460492503131e-16
#endif


// Macros

/**
  \def IS_ZERO(x)
  Evaluates to TRUE if |\a x| < SMALLNUM, FALSE otherwise.
 */
#define IS_ZERO(x) (fabs(x) < (SMALLNUM) ? (TRUE) : (FALSE))

/** 
  \def STREQ(s1,s2)

  Compares strings \a s1 and \a s2; evaluates to TRUE if strings are equal, FALSE if not.

 */
#define STREQ(s1,s2) (strcmp(s1,s2) == 0)

/** 
  \def SQR(x)
  Computes the square of \a x.
 */
#define SQR(x) ((x) * (x))

/**
  \def MAX(x,y)
  Computes the maximum of \a x and \a y
 */
#define MAX(x,y) ((x) > (y) ? (x) : (y))

/**
  \def MIN(x,y)
  Computes the minimum of \a x and \a y
 */
#define MIN(x,y) ((x) < (y) ? (x) : (y))

/* Define a macro for calculating the 1D index of the arrays
 * from the 2D indices i and j, in order to make the code simpler.
 * The 2D indices are: i corresponds to z and j corresponds to r.
 * The number of rows in most of these matrices is (nr+1)
 */
/**
  \def INDEX(i,j)
  Computes the 1D index for concentration or alpha arrays, given pseudo-2D indices \a i (z index) and \a j (r index); uses column-major ordering.
 */
#define INDEX(i,j) ((i)*(nr+1)+(j))


// Function prototypes

// convo.c
void convolve3(int M, int N, double *a, double scale1, double scale2, double *invr, double *out);

// extras.c
void error(char *errorstring, ...);

void print_usage_fit_layer(char *program);

void check_filename(char *in, char *out);

double *create_array(int N, char *string);

int assemble_command(int argc, char *argv[], char *command);

// model.c
void calc_diffusion_curve_layer_fit_layer(int nt, int nz, int nr, int iprobe, int jprobe, int iz1, int iz2, int nolayer, double dt, double dr, double sd, double st, double alpha_so, double theta_so, double kappa_so, double alpha_sp, double theta_sp, double kappa_sp, double alpha_sr, double theta_sr, double kappa_sr, double dfree, double *t, double *s, double *invr, double *p);

