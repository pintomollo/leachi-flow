#include <math.h>
#include <string.h>
#include "ctmf.h"
#include "mex.h"

#include "ctmf.c"

/* We need ot provide the size of the available memory, here 3Gb. */
#define MEM_SIZE 3*1024*1024

/*
 * The Matlab wrapper for the C code from ctmf.c, implementing constant-time median
 * filtering (see median_mex.m).
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  /* Declaring the variables with their default values. */
  int h, w, nelem, i, niter = 1, radius = 1; 
  unsigned char *tmp_img, *median_img, *tmp_ptr;
  double *img, value, mymin, mymax, scaling_factor;

  /* We accept either 1, 2 or 3 input arguments, always in the same order.
   * 1. The image 2. the radius of the kernel 3. the number of iterative calls. */
  if (nrhs < 1) {
    mexErrMsgTxt("Not enough input arguments (1 is the minimum, 3 is the maximum) !");
  } else if (nrhs == 2) {
    radius = (int) mxGetScalar(prhs[1]);
  } else if (nrhs == 3) {
    radius = (int) mxGetScalar(prhs[1]);
    niter = (int) mxGetScalar(prhs[2]);
  }

  /* Get the dimensions of the image. */
  img = mxGetPr(prhs[0]);
  h = mxGetM(prhs[0]);
  w = mxGetN(prhs[0]);
  nelem = h*w;

  /* Allocate memory for the result and the computations. */
  if ((median_img = mxCalloc(nelem, sizeof(unsigned char))) == NULL) {
    mexErrMsgTxt("Memory allocation failed !");
  }
  if ((tmp_img = mxCalloc(nelem, sizeof(unsigned char))) == NULL) {
    mexErrMsgTxt("Memory allocation failed !");
  }

  /* We need to convert our most-lekely double precision to UINT8.
   * So we first need to find the range of values present. */
  mymax = mymin = img[0];
  for (i = 1; i < nelem; i++) {
    if (img[i] < mymin) {
      mymin = img[i];
    } else if (img[i] > mymax) {
      mymax = img[i];
    }
  }

  /* Then we compute the scaling factor. */
  scaling_factor = 255 / (mymax - mymin);

  /* And we convert the image, setting NaN to 0. */
  for (i = 0; i < nelem; i++){
    if (mxIsNaN(img[i])) {
      median_img[i] = 0;
    } else {
      value = ceil(scaling_factor*(img[i] - mymin));

      median_img[i] = (unsigned char) value;
    }
  }

  /* Now let's call this function as many times as requried. */
  for (i = 0; i < niter; i++) {
    tmp_ptr = tmp_img;
    tmp_img = median_img;
    median_img = tmp_ptr;

    /* Here we enforce single-channel analysis and hard coded the memory size. */
    ctmf(tmp_img, median_img, h, w, h, h, radius, 1, MEM_SIZE);
  }

  /* Clear the temporary image. */
  mxFree(tmp_img);

  /* Create the output variable. */
  plhs[0] = mxCreateDoubleMatrix(h, w, mxREAL);
  img = mxGetPr(plhs[0]);

  /* Copy the image, rescaling it properly. */
  scaling_factor = 1/scaling_factor;
  for (i=0;i < nelem; i++) {
    img[i] = ((double) ((double) median_img[i]) * scaling_factor) + mymin;
  }

  /* Free the last image. */
  mxFree(median_img);

  return;
}
