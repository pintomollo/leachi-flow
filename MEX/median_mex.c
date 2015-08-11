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
  mwSize h, w, c, dims[3];
  const mwSize *size;
  int nelem, i, j, niter = 1, radius = 1; 
  unsigned char *tmp_img, *median_img, *tmp_ptr;
  double *img, *img_out, value, mymin, mymax, scaling_factor;

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

  // The size of the image
  size = mxGetDimensions(prhs[0]);
  h = size[0];
  w = size[1];
  c = mxGetNumberOfElements(prhs[0]) / (h*w);

  // Prepare the output
  dims[0] = h;
  dims[1] = w;
  dims[2] = c;

  // Some temporary values
  nelem = h*w;

  // Create the output array
  plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  img_out = mxGetPr(plhs[0]);

  /* Get the input image. */
  img = mxGetPr(prhs[0]);

  /* Allocate memory for the result and the computations. */
  if ((median_img = mxCalloc(nelem, sizeof(unsigned char))) == NULL) {
    mexErrMsgTxt("Memory allocation failed !");
  }
  if ((tmp_img = mxCalloc(nelem, sizeof(unsigned char))) == NULL) {
    mexErrMsgTxt("Memory allocation failed !");
  }

  // Loop over the stack separatly
  for (j=0; j < c; j++) {

    /* We need to convert our most-likely double precision to UINT8.
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

    /* Copy the image, rescaling it properly. */
    scaling_factor = 1/scaling_factor;
    for (i=0;i < nelem; i++) {
      img_out[i] = ((double) ((double) median_img[i]) * scaling_factor) + mymin;
    }

    // Update the image pointers
    img += nelem;
    img_out += nelem;
  }

  /* Free the two image. */
  mxFree(tmp_img);
  mxFree(median_img);

  return;
}
