/*=================================================================
* createSparse.c
*=================================================================*/

/* $Revision: 1.5.6.3 $ */

#include <math.h> /* Needed for the ceil() prototype */
#include "mex.h"
#include "mkl.h"
/*#include "mkl_lapacke.h"*/

/* If you are using a compiler that equates NaN to be zero, you must
 * compile this example using the flag  -DNAN_EQUALS_ZERO. For example:
 *
 *     mex -DNAN_EQUALS_ZERO fulltosparse.c
 *
 * This will correctly define the IsNonZero macro for your C compiler.
 */

#if defined(NAN_EQUALS_ZERO)
#define IsNonZero(d) ((d)!=0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d)!=0.0)
#endif

void mexFunction(
		 int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]
		 )
{
    /* Declare variable */
    double *a, *b;
	char trans='N';
    MKL_INT m, n, nbrhs, lda, ldb, *ipvt,info=0;
    /* Check for proper number of input and output arguments */
    if (nrhs != 3) {
        mexErrMsgIdAndTxt( "MATLAB:fulltosparse:invalidNumInputs",
                "Two input arguments required.");
    }
    if(nlhs > 1){
        mexErrMsgIdAndTxt( "MATLAB:fulltosparse:maxlhs",
                "Too many output arguments.");
    }
    
        
    /* Get the size and pointers to input data */
    a = mxGetPr(prhs[0]);
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	ipvt = (MKL_INT*)mxGetData(prhs[1]);
	plhs[0] = mxDuplicateArray(prhs[2]);
	nbrhs = mxGetN(plhs[0]);
	b = mxGetPr(plhs[0]);
	lda = m;
	ldb = m;
	dgetrs_(&trans,&m,&nbrhs,a,&lda,ipvt,b,&ldb,&info);
    
    
    return;
}
