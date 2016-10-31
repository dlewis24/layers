/**
  \file 3layer/io.c

  Input/output functions that are useful for 3layer

  \author David Lewis, CABI, NKI
  \copyright GNU Public License
  \date 2012-2013

 */

// includes
#include "header.h"

/**
  \brief Make sure that the filename is not too long.

   If the input string is not to big, copy it to the output string. 
   Otherwise, exit with an error message.  This function is used 
   to check that the user input filename is not too long. 

  \param [in] in Input string = the user input filename (or base filename)
  \param [out] out Output string = the input string, if not too long
 */
void get_filename(char *in, char *out)
{
	// The input and output filenames are 4 characters longer than 
	// the base filename, hence the "- 4"
	if (strlen(in) < (FILENAME_MAX - 4)) 
		strcpy(out, in);
	else 
		error("Filename length is too long");
}


/**
  \brief Determine the input filename and the default output file name.

  From the last argument on the command line, determine 
  the input file name and the default output file name. 
  (If the output file name is specified on the command 
  line, it will be used instead of the default output 
  file name determined here.)

  If the argument has an extension, it is the input filename;
  remove the extension to get the basename and add the default
  output filename extension to get the default output filename.

  If the argument has no extension, it is the basename; then
  add the default extensions to get the input filename and the
  default output filename.

  \param [in] argstring Final command-line argument
  \param [in] inf_extension Default extension for input file
  \param [in] outf_extension Default extension for output file
  \param [out] infilename Name of input file
  \param [out] outfilename Name of output file

 */
void get_io_filenames(char *argstring, const char *inf_extension, const char *outf_extension, char *infilename, char *outfilename)
{
	char *strng;

	get_filename(argstring, infilename);
	strcpy(outfilename, infilename);

	strng = index(outfilename, '.');  /* Find the '.' in the filename */
	if (strng == NULL) {    /* If no '.', then add default extensions */
		strcat(infilename, inf_extension);
		strcat(outfilename, outf_extension);
	} else {
		strcpy(strng, outf_extension);  /* Replace ext. with default */
	}
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
 */
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


/*******************************************************************
	This function is for parsing a string which has several 
	doubles delimited by spaces or commas. It uses strtok. 
	The string has parameters for the extra sources. 

	For example, 
		--additional_sources "2 50.0 0.0 100.0 -50.0 0.0 100.0"
		-> 2 means there are 2 additional sources 
		-> +/-50.0, 0.0, and 100.0 mean the sz, sr, and crnt 
		   values for the sources

	Input: string describing the quantity (sz, sr, or crnt) 
		and # of the extra source -- both for diagnostics
	Return: value
*******************************************************************/
/**
  \brief Read a parameter from the string argument to the additional_sources option.

  The additional_sources option to 3layer takes several
  arguments. (The total number of arguments depends on the number
  of extra sources.) The getopt_long function (used in the main()
  function to parse 3layer's options) takes only one argument per
  option, so for additional_sources's argument use a string that
  contains the several arguments that additional_sources needs.

  This function uses strtok to parse the string to obtain the next
  parameter. The two arguments for this function are used only for 
  diagnostic purposes -- if there is an error, the function aborts 
  with an error message that lists these two arguments.

  The first element of the string (the number of extra sources) 
  was previously obtained with the first call to strtok, and 
  subsequent calls to strtok (via this function) use NULL rather 
  than the name of the string.

  \param [in] string String containing the type of parameter to be read
  \param [out] nsource The number of the source whose parameter is being read

  \return Value of the parameter that is read
 */
double read_source_parameter(char *string, int nsource)
{
	double value;
	char *token;

	token = strtok(NULL, " ,");
	if (token == NULL)
		error("Cannot read %s token; nsource = %d\n", string, nsource);

	value = atof(token);

	return (value);
}
