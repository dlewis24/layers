/** 
  \file fit-layer/convo.c

  Function for calculating the Laplacian of the concentration in 
  cylindrical coordinates via convolutions.

The change in the concentration with time depends on the Laplacian of 
the concentration.  In cylindrical coordinates, the Laplacian of 
\f$ c \f$ is 

\f[ \nabla^2{c} = \frac{\partial^2 c}{\partial z^2} 
                + \frac{\partial^2 c}{\partial r^2}
                + \frac{1}{r}   \frac{\partial c}  {\partial r}
                + \frac{1}{r^2} \frac{\partial^2 c}{\partial \phi^2} \f]

Because of circular symmetry, \f$ \partial^2c/\partial \phi^2 = 0 \f$ .

Note that in this program \f$ \Delta z = \Delta r \f$, which simplifies 
the kernels below.

The first 2 terms in the equation for the Laplacian are implemented 
numerically as if they were a 2D Laplacian in cartesian coordinates, 
by convolving the concentration array \f$ c \f$ with the Laplacian kernel
\f[ \mathbf{L} = 
\left( \begin{array}{ccc}
0 &  1 & 0 \\
1 & -4 & 1 \\
0 &  1 & 0
\end{array} \right) \f]

and dividing by \f$ \Delta r^2 \; (= \Delta z^2) \f$ 
(via a scaling factor).

The first derivative with respect to \f$ r \f$ in the third term is 
implemented numerically by convolving the \f$ c \f$ array 
with the derivative kernel

\f[ \mathbf{D} = 
\left( \begin{array}{ccc}
0 & 0 & 0 \\
-1 & 0 & 1 \\
0 & 0 & 0
\end{array} \right) \f]

and dividing by \f$ 2 \Delta r \f$ (via another scaling factor). 
Rows of the concentration array correspond to \f$ r \f$ and columns 
correspond to \f$ z \f$.

In the third term, \f$ 1/r \f$ is a problem numerically for \f$ r = 0 \f$ . 
But by L'Hopital's Rule, 
	\f[\lim_{r \to 0}  (\partial c/\partial r) / r = \partial ^2c/\partial r^2  \f]
So for the case \f$ r = 0 \f$ we use 

\f[ \nabla^2{c} = \partial ^2c/\partial z^2 + 2 \partial ^2c/\partial r^2  \qquad (r = 0) \f]

which is implemented by convolving \f$ c \f$ with the kernel
\f[ \mathbf{L0} = 
\left( \begin{array}{ccc}
0 &  1 & 0 \\
2 & -6 & 2 \\
0 &  1 & 0
\end{array} \right) \f]

Also the Laplacian is scaled by \f$ \theta D_{free} \f$
via the scaling factors.

  \author David Lewis, CABI, NKI
  \copyright GNU Public License
  \date 2012-2013

 */

#include <stdio.h>
#include <stdlib.h>

/**
  \def A(i,j)
  Retrieves the appropriate element from the 1D 'a' array given pseudo-2D indices \a i (z index) and \a j (r index), for the convolutions; uses column-major ordering.
 */
#define A(i,j) a[((i)*N+(j))]

/**
  \def OUT(i,j)
  Retrieves the appropriate element from the 1D 'out' array given pseudo-2D indices \a i (z index) and \a j (r index); uses column-major ordering.
 */
#define OUT(i,j) out[((i)*N+(j))]


/**
  \brief Calculates updates to the concentration in a layer 
         by applying the Laplacian in cylindrical coordinates; 
         also scales by \f$ \theta D_{free} \f$.

  Calculates the update terms

\f[
      \mathbf{out} = s_1 \mathbf{L} * \mathbf{a} + \frac{s_2}{r} \mathbf{D} * \mathbf{a} \qquad (r \neq 0)
\f]

  or

\f[
      \mathbf{out} = s_1 \mathbf{L0} * \mathbf{a} \qquad (r = 0)
\f]

where

      \f$ \mathbf{out} \f$ = output array

      \f$ s_1, s_2 \f$ = scaling factors

      \f$ \mathbf{a} \f$ = input array (concentration)

      \f$ * \f$ means convolution

      \f$ \mathbf{L} \f$ = 3x3 Laplace operator (see above)

      \f$ \mathbf{L0} \f$ = 3x3 Laplace operator modified for \f$ r = 0 \f$

      \f$ \mathbf{D} \f$ = 3x3 derivative operator for \f$ r \f$

Note that the index \f$ j = 1 \f$ corresponds to \f$ r = 0 \f$.

  \param [in] M Number of columns of input matrix (z)
  \param [in] N Number of rows of input matrix (r)
  \param [in] a Input matrix (concentration, c)
  \param [in] scale1 Scaling factor #1
  \param [in] scale2 Scaling factor #2
  \param [in] invr Vector of 1/r values
  \param [out] out Output matrix
 */

void convolve3(int M, int N, double *a, double scale1, double scale2, double *invr, double *out)
{
	int i, j;

	// Calculate the inner part of the array first, 
	// then the edges and corners later 
	for (i=1; i<M-1; i++) 
		for (j=2; j<N-1; j++) 
			OUT(i,j) = scale1 * (
				                  A(i-1,j) 
				+ A(i,j-1) - 4. * A(i,  j) + A(i,j+1)
				                + A(i+1,j)   )
			+ scale2 * ( (-A(i,j-1) + A(i,j+1))*invr[j]);

	/* j=1: This is the r=0 row, so use L0 rather than L */
	for (i=1; i<M-1; i++) 
		OUT(i,1) = scale1 * (
				                     A(i-1,1) 
				+ 2. * A(i,0) - 6. * A(i,  1) + 2. * A(i,2)
				                   + A(i+1,1)   );

	// j=1, i=0 
	OUT(0,1) = scale1 * (
				2. * A(0,0) - 6. * A(0,1) + 2. * A(0,2)
				                 + A(1,1)   );

	// j=1, i=M-1 
	OUT(M-1,1) = scale1 * (
				                     A(M-2,1) 
				+ 2. * A(M-1,0) - 6. * A(M-1,  1) + 2. * A(M-1,2)
				                   );
	/* End of j=1 section */


	// i=0 
	for (j=1; j<N-1; j++) 
		OUT(0,j) = scale1 * (
			  A(0,j-1) - 4. * A(0,j) + A(0,j+1)
			                + A(1,j)   )
			+ scale2 * ( (-A(0,j-1) + A(0,j+1))*invr[j] );

	// i=M-1 
	for (j=1; j<N-1; j++) 
		OUT(M-1,j) = scale1 * (
			                     A(M-2,j) 
			+ A(M-1,j-1) - 4. * A(M-1,j) + A(M-1,j+1)   )
			+ scale2 * ( (-A(M-1,j-1) + A(M-1,j+1))*invr[j] );

	// j=0 
	for (i=1; i<M-1; i++) 
		OUT(i,0) = scale1 * (
			       A(i-1,0) 
			- 4. * A(i,  0) + A(i,1)
			     + A(i+1,0)   )
			+ scale2 * ( A(i,1)*invr[0] );

	// j=N-1 
	for (i=1; i<M-1; i++) 
		OUT(i,N-1) = scale1 * (
			                   A(i-1,N-1) 
			+ A(i,N-2) - 4. * A(i,  N-1) 
			                 + A(i+1,N-1)   )
			+ scale2 * ( -A(i,N-2)*invr[N-1] );

	// i=0, j=0 
	OUT(0,0) = scale1 * (
			- 4. * A(0,0) + A(0,1)
			     + A(1,0)   )
			+ scale2 * ( A(0,1)*invr[0] );

	// i=0, j=N-1 
	OUT(0,N-1) = scale1 * (
			  A(0,N-2) - 4. * A(0,N-1) 
			                 + A(1,N-1)   )
			+ scale2 * ( -A(0,N-2)*invr[N-1] );

	// i=M-1, j=0 
	OUT(M-1,0) = scale1 * (
			       A(M-2,0) 
			- 4. * A(M-1,0) + A(M-1,1)   )
			+ scale2 * ( A(M-1,1)*invr[0] );

	// i=M-1, j=N-1 
	OUT(M-1,N-1) = scale1 * (
			                      A(M-2,N-1) 
			+ A(M-1,N-2) - 4. * A(M-1,N-1)   )
			+ scale2 * ( -A(M-1,N-2)*invr[N-1] );
}
