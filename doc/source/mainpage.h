// For main Doxygen page

/**
\mainpage Introduction 

<I> This documentation is generated automatically by Doxygen from the
source code comments found in the *.c and *.h files.  In particular,
every external (variable, macro, function, etc.) has a documentation
entry. </I>

This software implements the multilayer model of extracellular
diffusion described in (\ref saghyan "Saghyan et al., 2012") for
the case of 3 adjacent homogeneous layers of tissue.  The software
is useful in modeling or analyzing data from extracellular
diffusion measurements made with the Real-Time Iontophoretic
(RTI) method (\ref nicholson "Nicholson and Phillips, 1981")
of diffusion in 3 distinct layers that have different diffusion
parameters.  The two programs, 3layer and fit-layer, are
described below.

The implementation described in (Saghyan et al., 2012) was in 
IDL; this implementation is in C.  This implementation has some 
improvements over the IDL implementation, for example the 
addition of linear nonspecific clearance in the model.

Each layer is assumed to have a constant extracellular volume
fraction \f$\alpha\f$, a constant diffusion permeability
\f$\theta\f$, and a constant nonspecific clearance factor
\f$\kappa\f$.  The 3 layers we had in mind when writing this
software are stratum oriens (SO), stratum pyramidale (SP), and
stratum radiatum (SR) of the CA1 region of the hippocampus,
which is why the names of the diffusion parameters end in
_so, _sp, and _sr.

The environment modeled in the software consists of a cylinder of
tissue through the layers.  The axis of the cylinder goes through
the point source and is perpendicular to each layer.  The boundary
conditions are zero concentration (total absorption) at the
boundaries of the cylinder (top, bottom, and side).  Cylindrical
coordinates are used in order to take advantage of cylindrical
symmetry and reduce computation times.

The software uses the Forward Time Centered Space (FTCS) finite
difference scheme to solve the diffusion equation.  This simple
fully explicit scheme avoids having to deal with large sparse
matrices.  However, many time-steps must be used in order to meet 
the von Neumann stability criterion, \f$ dt < dr^2 / (6D) \f$ ,
which can result in long program run times, especially for 
the program fit-layer.

Note that the run time of each program depends heavily on the
size of the grid used for the finite difference calculations
(the \f$nr\f$ and \f$nz\f$ parameters).

References:

\anchor nicholson
  Nicholson C, Phillips JM (1981) Ion diffusion modified by 
  tortuosity and volume fraction in the extracellular 
  microenvironment of the rat cerebellum.  J Physiol 321:225–57.

\anchor saghyan 
  Saghyan A, Lewis DP, Hrabe J, Hrabetova S (2012) Extracellular 
  diffusion in laminar brain structures exemplified by hippocampus.  
  J Neurosci Methods 205:110-8.  Epub 2011 Dec 30.



\section license_sec License

Copyright (C) 2012-2013  
<a href="mailto:lewis@nki.rfmh.org">David Lewis</a>, CABI, NKI

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
<a href="http://www.gnu.org/licenses/old-licenses/gpl-2.0.html#SEC1"
>GNU General Public License</a> for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  
02110-1301, USA.


\section prereq_sec Prerequisites

This software runs on Linux. It was developed on CentOS 5.6
and tested on Centos 6.0, Debian 6.0.4, and Gentoo, and it
should compile and run on a recent version of any other Linux
distribution. In order to compile and run the software, the
following dependencies are needed:

  - A C compiler such as gcc which supports C99
  - The Make utility, in order to be able to use a makefile to compile 
    the software
  - The GNU Scientific Library, including header files (packages gsl 
    and gsl-devel on CentOS)

The software outputs data in ASCII format, so it can be easily
read into a plotting program like Gnuplot.


\section unpacking_sec Unpacking the Software

After the tar archive is downloaded, the software can be extracted
with WinZip (Windows) or the tar utility (Linux), if it hasn't
already been extracted during the download. For example in Linux,

<code>
$ <b>tar -xf layers-1.2.tar</b><br>
$ <b>cd layers-1.2</b>
</code>



\section install_sec Installation

There is no special installation procedure. Individual programs
are compiled separately. In order to avoid the user's having to
specify the absolute or relative path each time a program is run,
the programs can be copied to a directory in the user's path, or 
the PATH environment variable can be appropriately modified.



\section three_layer_sec 3layer

The program 3layer calculates the extracellular concentration
of a substance as a function of time from a point source
embedded in an environment comprised of 3 homogeneous layers.
The diffusion parameters in each layer, size of the cylinder,
thickness of each layer, etc., are given.

This program outputs the extracellular concentration at the probe
as a function of time, in order to model data from RTI diffusion
measurements.  It optionally outputs images showing the spatial
distribution of the concentration at regularly-spaced time points 
in order to model data from Integrative Optical Imaging (IOI) 
diffusion measurements.

\subsection three_layer_compiling_sec Compiling

To compile the program, run 'make' in the 3layer directory:

<code>
$ <b>cd 3layer</b><br>
$ <b>make 3layer</b><br>
...
</code>


\subsection three_layer_running_sec Running

When you run 3layer you specify an input file of parameters,
and when the program finishes it will write the calculated
concentration to an output file.  By default the output file
will have the same basename as but a different extension than
the input file.  You can get a usage statement and a list of
command-line options by running 3layer without specifying the
input file:

<code>
$ <b>./3layer</b><br>
Usage: 3layer [options] \<input_file\><br>
...
</code>

(On Linux/Unix, the characters "./" preceding the program name mean
"run the program from the current directory" and are not necessary 
if the current path is in your PATH environment variable.)

The program 3layer can take several minutes or hours to run,
depending on the size of the grid and the speed of your computer.


\subsection three_layer_input_sec Input File

Parameters to the program (such as the size of the environment and
the diffusion parameters) can be specified in the input file or on
the command line.  If any parameter is specified both in the input
file and on the command line, the parameter value specified on the
command line takes precedence.  If any parameter is not specified,
its default value from the program is used.  Parameter assignments
in the input file are of the form \<parameter name\> = \<value\>.

The 3layer directory has an example input file called "sample.par".
The comments in sample.par describe the input file format in 
more detail.

The location of the source defines the origin of the coordinate
system used for input.  As a result, positions specified in the
input file or on the command line (such as the probe position)
are relative to the source position; see the file coordinates.pdf.
This convention was chosen to be consistent with that of the
program fit-layer.

\subsection three_layer_output_sec Output File

The output file includes a list of parameters and concentration
data from the run.  Note that some of the input parameters get
adjusted by the program for a number of reasons:

- Shift of z-coordinates:  The program uses a coordinate system 
  where the bottom of the cylinder defines \f$z=0\f$ and the
  input file uses a coordinate system where the source defines
  the origin.

- Discretization:  The source and the probe positions are 
  adjusted to fall on grid points, and the layer boundaries 
  are adjusted to fall midway between the grid points.

- Equal resolutions in \f$r\f$ and \f$z\f$:  The cylinder radius 
  is adjusted if necessary to make the spatial resolutions in 
  both the \f$r\f$- and \f$z\f$-directions equal.

The output file includes the extracellular concentration
(not the tissue concentration) at the probe as a function
of time (columns 1 and 2) and a "characteristic curve" fit
to this data using the traditional model used in RTI analysis
(column 3).  The traditional model assumes that the environment
consists of one homogeneous layer.  Although the characteristic
curve often closely matches the curve from the 3-layer model,
the corresponding fit parameters (the "apparent parameters")
are not physically meaningful, because the assumption of one 
homogeneous region is incorrect.

The 3layer directory also has an example output file,
"sample.dat.orig".  It was generated from 3layer with
the command

<code>
$ ./3layer --outfile sample.dat.orig sample.par
</code>


\subsection three_layer_testing_sec Testing

You can run

<code>
$ <b>./3layer sample.par</b>
</code>

and compare your output file sample.dat with sample.dat.orig.
The program might take a few minutes to run, depending on the
speed of your computer.

The two output data files should differ in the following lines:
	-# The line listing the command used
	-# The lines listing the start time, end time, and total time
	-# The Solution lines will probably differ in the last 3 numbers, 
	   since they are the total time in seconds, minutes, and hours

If the parameters or the concentration data of these two output
files are substantially different for this test, please report it
to the program's author.


\subsection three_layer_graphing_sec Graphing

The output files are in a format that can be graphed with
gnuplot.  For example,

<code>
$ <b>gnuplot</b><br>
...<br>
gnuplot> <b>set style data lines</b><br>
gnuplot> <b>set xlabel "Time (s)"</b><br>
gnuplot> <b>set ylabel "Concentration (mM)"</b><br>
gnuplot> <b>set title "Diffusion curve calculated from multilayer model"</b><br>
gnuplot> <b>plot [0:150][0:1] "sample.dat" using 1:2 title "sample.dat"</b>
</code>

If you plot the diffusion curve from the test run with
the original data that was provided, the curves should be
indistinguishable:

<code>
gnuplot> <b>set title "Comparison of new data to original data"</b><br>
gnuplot> <b>plot [0:150][0:1] "sample.dat" using 1:2 title "sample.dat", "sample.dat.orig" using 1:2 title "sample.dat.orig"</b><br>
</code>

The diffusion curve looks very similar to a typical RTI diffusion
curve from a homogeneous region.  In fact, if you fit diffusion
data from multiple layers with the traditional model used in RTI,
which assumes that the diffusion occurs in one homogeneous region,
often the curve from the fit will match the data very well.
Such a curve is called a "characteristic curve" in (Saghyan et
al., 2012).  However, since the assumption of one homogeneous
region is incorrect, the ECS parameters determined in such a
fit are not physically meaningful.  The 3layer program fits the
curve from the multilayer model with such a characteristic curve
for comparison:

<code>
gnuplot> <b>set title "Comparison with characteristic curve"</b><br>
gnuplot> <b>plot [0:150][0:1] "sample.dat" using 1:2 title "3-Layer Model", "sample.dat" using 1:3 title "Characteristic Curve"</b><br>
gnuplot> <b>exit</b><br>
  <br>
$ <b>grep Fit sample.dat</b><br>
# Fit for characteristic curve:<br>
# Fitted apparent alpha = 0.264998<br>
# Fitted apparent theta = 0.340114  (lambda = 1.714698)
</code>

So if you fit data from a region with multiple layers using the
traditional model, the fit parameters can be very inaccurate.
A multilayer model is needed for an accurate fit.  For more
detail see (Saghyan et al., 2012).  The program fit-layer does
such a multilayer fit.


\section fit_layer_sec Fit-layer

The program fit-layer fits the multilayer model to RTI data
obtained in a region with 3 homogeneous layers in order
to determine the diffusion parameters in the middle layer.
The inverse problem is implemented with the downhill simplex
algorithm, using for the forward problem the same multilayer 
model that 3layer uses.

As mentioned above, traditional analysis of RTI data assumes
that the measurement region is homogeneous.  If the region is not
homogeneous, the diffusion parameters obtained from the fitting
procedure can be very inaccurate.  For example, in the CA1 region
of hippocampus we find that traditional analysis is adequate for
determining the parameters in SO and SR, since the measurements
in these layers can be made away from adjacent layers.  However,
RTI diffusion measurements in the thin SP layer are confounded
by diffusion in the adjacent SO and SR layers, so multilayer
analysis must be used to determine the diffusion parameters in SP.

\subsection fit_layer_compiling_sec Compiling

To compile the program, run 'make' in the fit-layer directory:

<code>
$ <b>cd ../fit-layer</b><br>
$ <b>make fit-layer</b><br>
...
</code>


\subsection fit_layer_running_sec Running

When you run fit-layer you specify an input file of parameters and
data, and when the program finishes it will write the fitted SP
parameters and the corresponding concentration to an output file.
By default the output file will have the same basename as but
a different extension than the input file.  You can get a usage
statement and a list of command-line options by running fit-layer
without specifying the input file:

<code>
$ <b>./fit-layer</b><br>
Usage: fit-layer [options] \<input_file\><br>
...<br>
</code>

The program fit-layer can take several minutes, hours, or even
weeks to run, depending on the size of the grid and the speed
of your computer.


\subsection fit_layer_input_sec Input File

The parameter assignment section of an input file for fit-layer
is similar to an input file for 3layer.  Parameters to the program
(such as the size of the environment and the diffusion parameters)
can be specified in the input file or on the command line.
If any parameter is specified both in the input file and on the
command line, the parameter value specified on the command line
takes precedence.  If any parameter is not specified, its default
value from the program is used.

The parameter assignment section is followed by two blank 
lines and then followed by the data section.  The data section 
has 2 columns:  time (s) and concentration (mM).  It can have 
other columns, which are ignored.  The first line of the data 
section is the heading.

The fit-layer directory has an example input file called 
"data.txt".  Data for this input file were taken from the output 
of 3layer, specifically from the file sample.dat.orig.  Parameters 
for this input file were taken from sample.par, except for the SP 
parameters, which are to be determined by fit-layer.  The comments 
in data.txt describe the input file format in more detail.

The location of the source defines the origin of the coordinate 
system used for input.  As a result, positions specified in the 
input file or on the command line (such as the probe position) 
are relative to the source position; see the file coordinates.pdf. 

This convention was chosen because in an RTI experiment 
distances from the source to the probe and to the layer 
boundaries are measured.  It is natural for an experimenter 
to enter these values in the input file, rather than to 
specify positions relative to a fictitious boundary, 
such as the bottom of the cylinder.

\subsection fit_layer_output_sec Output File

The output file includes a list of parameters used for the fit.
Note that some of the input parameters get adjusted by the
program for a number of reasons:

- Shift of z-coordinates:  The program uses a coordinate system 
  where the bottom of the cylinder defines \f$z=0\f$ and the
  input file uses a coordinate system where the source defines
  the origin.

- Discretization:  The source and the probe positions are 
  adjusted to fall on grid points, and the layer boundaries 
  are adjusted to fall midway between the grid points. 

- Equal resolutions in \f$r\f$ and \f$z\f$:  The cylinder radius 
  is adjusted if necessary to make the spatial resolutions in 
  both the \f$r\f$- and \f$z\f$-directions equal.

The output file also includes results of the fitting:  the fitted 
values of \f$\alpha\f$, \f$\theta\f$, and \f$\kappa\f$ of the SP 
layer.  Also included are details from the minimization algorithm 
used, such as the number of iterations, the size of the final 
simplex, and the total time for the program run.

The final section of the output file has concentration data 
for the probe location.  The first two columns are concentration 
(column 2) as a function of time (column 1).  This concentration 
was calculated using the multilayer model using the fitted 
parameters for the SP layer and the given parameters for the 
other layers.  The last two columns are for the input concentration 
(column 4) as a function of time (column 3) -- i.e., the 
concentration data taken from the input file.  

The fit-layer directory also has an example output file,
"data.dat.orig".  It was generated from fit-layer with
the command

<code>
$ ./fit-layer --nr 100 --nz 200 --outfile data.dat.orig data.txt
</code>

You can search for the fitted parameters with grep:

<code>
$ <b>grep Fitted data.dat.orig</b><br>
...
</code>


\subsection fit_layer_testing_sec Testing

You can run

<code>
$ <b>./fit-layer --nr 100 --nz 200 data.txt</b>
</code>

and check the fit

<code>
$ <b>grep Fitted data.dat</b><br>
# Fitted alpha = 0.102433<br>
# Fitted theta = 0.294797  (lambda = 1.841783)<br>
# Fitted kappa = 0.000178 s^-1<br>
</code>

and compare your output file data.dat with data.dat.orig.
The parameters from the fit correspond to the SP layer.  They are
close to the true values of \f$\alpha\f$ = 0.10, \f$\theta\f$ =
0.30, and \f$\kappa\f$ = 0.0 s^-1 for this test (values from
sample.par).  A small grid size (100 x 200) is used to reduce the
total runtime for this test.  Nevertheless, the program might take
several minutes to run, depending on the speed of your computer.

The two output data files should differ in the following lines:
	-# The line listing the command used
	-# The lines listing the start time, end time, and total time
	-# The Solution lines will probably differ in the last 3 numbers, 
	   since they are the total time in seconds, minutes, and hours

If the parameters or the concentration data of these two output
files are substantially different for this test, please report it
to the program's author.


\subsection fit_layer_graphing_sec Graphing

The output files are in a format that can be graphed with
gnuplot.  For example,

<code>
$ <b>gnuplot</b><br>
...<br>
gnuplot> <b>set style data lines</b><br>
gnuplot> <b>set xlabel "Time (s)"</b><br>
gnuplot> <b>set ylabel "Concentration (mM)"</b><br>
gnuplot> <b>set title "Fit of multilayer model to input data"</b><br>
gnuplot> <b>plot [0:150][0:1] "data.dat" using 1:2 title "Data from fit", "data.dat" using 3:4 title "Input data"</b>
</code>

If you plot the diffusion curve from the test run with
the original data that was provided, the curves should be
indistinguishable:

<code>
gnuplot> <b>set title "Comparison of new data to original data"</b><br>
gnuplot> <b>plot [0:150][0:1] "data.dat" using 1:2 title "data.dat", "data.dat.orig" using 1:2 title "data.dat.orig"</b>
</code>




\author <a href="mailto:lewis@nki.rfmh.org">David Lewis</a>, CABI, NKI
\author Based on IDL software written by Jan Hrabe and David Lewis, CABI, NKI
\date 2012-2013
*/