#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include "gaussian_smooth.h"
#include "mex.h"

#include "gaussian_smooth.c"

/* A Matlab wrapper for the code of Anthony Gabrielson (see gaussian_smooth.c)*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  /* Declare a few variables. */
  int h, w;
  double *img, sigma;

  /* No flexibility here, we want both the image and sigma ! */
  if (nrhs < 2) {
    mexErrMsgTxt("Not enough input arguments (2 are required) !");
  } else if (!mxIsDouble(prhs[0])) {
    mexErrMsgTxt("Input array is not of type Double");
  }

  /* Get sigma. */
  sigma = mxGetScalar(prhs[1]);

  /* Get the size of the image. */
  h = mxGetM(prhs[0]);
  w = mxGetN(prhs[0]);

  /* Create the output image on which we'll work directly. */
  plhs[0] = mxCreateDoubleMatrix(h, w, mxREAL);
  img = mxGetPr(plhs[0]);

  /* Copy the input to the working image. */
  memcpy(img, mxGetPr(prhs[0]), h*w*sizeof(double)); 

  /* Verify that sigma is valid, and let's go ! */
  if (sigma <= 0) {
    mexWarnMsgTxt("Gaussian smoothing with invalid sigma !");
  } else {
    gaussian_smooth(img, w, h, sigma);
  }

  return;
}
