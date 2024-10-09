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
    /* Declare variable*/
    double *pr;
    MKL_INT m, n, lda, *ipvt,kl,ku,info=0;
    /* Check for proper number of input and output arguments */
    if (nrhs != 4) {
        mexErrMsgIdAndTxt( "MATLAB:fulltosparse:invalidNumInputs",
                "Four input arguments required.");
    }
    if(nlhs > 0){
        mexErrMsgIdAndTxt( "MATLAB:fulltosparse:maxlhs",
                "Too many output arguments.");
    }
    
        
    /* Get the size and pointers to input data */
    pr = mxGetPr(prhs[0]);
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	ipvt = (MKL_INT*)mxGetData(prhs[1]);
	kl = mxGetScalar(prhs[2]);
	ku = mxGetScalar(prhs[3]);
	lda = m;
	
	dgbtrf_(&n,&n,&kl,&ku,pr,&lda,ipvt,&info);
    
    return;
}
