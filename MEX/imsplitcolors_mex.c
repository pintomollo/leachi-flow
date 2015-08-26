#include <math.h>
#include <string.h>
#include "mex.h"

// Bilinear interpolation, main interface
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // Declare variable
  int i, indx;
  mwSize w, h, npix, nbins, ndim;
  const mwSize *dims;
  double r = 0.0, g = 1.0/3.0, b = 2.0/3.0, val, t0, t1, t2, d0, d1, d2, dind;
  double *img, *table;

  // Check for proper number of input and output arguments
  if (nrhs != 2) {
    mexErrMsgIdAndTxt("MATLAB:imsplitcolors:invalidInputs",
        "Not the correct number of input arguments (2 required) !");
  }

  // Ensure the types of the two first arrays at least
  if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])) {
    mexErrMsgIdAndTxt("MATLAB:imsplitcolors:invalidInputs",
        "Input arguments must be of type double.");
  }

  // Get the number of pixels and the image
  npix = mxGetNumberOfElements(prhs[0]) / 3;

  // Get the number of bins and the counts table
  nbins = mxGetNumberOfElements(prhs[1]);

  // Then we need to split the colors
  if (nbins <= 3) {

    table = mxGetPr(prhs[1]);
    t0 = table[0];
    t1 = nbins > 1 ? table[1] : t0;
    t2 = nbins > 2 ? table[2] : t0;

    ndim = mxGetNumberOfDimensions(prhs[0]);
    dims = mxGetDimensions(prhs[0]);

    // Prepare the output, allocating the memory
    //if ((plhs[0] = mxCreateDoubleMatrix(npix, 3, mxREAL)) == NULL) {
    if ((plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL)) == NULL) {
      mexErrMsgIdAndTxt("MATLAB:imsplitcolors:memoryAllocation",
        "Memory allocation failed !");
    }
    img = mxGetPr(plhs[0]);

    memcpy(img, mxGetPr(prhs[0]), npix*3*sizeof(double)); 

    for (i=0; i < npix; i++) {
      val = img[i];

      d0 = val - t0;
      d1 = val - t1;
      d2 = val - t2;

      d0 = d0 > 0.5 ? d0 - 1 : d0;
      d1 = d1 > 0.5 ? d1 - 1 : d1;
      d2 = d2 > 0.5 ? d2 - 1 : d2;

      d0 = d0 < -0.5 ? d0 + 1 : d0;
      d1 = d1 < -0.5 ? d1 + 1 : d1;
      d2 = d2 < -0.5 ? d2 + 1 : d2;

      if (fabs(d0) <= fabs(d1)) {
        if (fabs(d0) <= fabs(d2)) {
          val = d0 + r;
        } else {
          val = d2 + b;
        }
      } else {
        if (fabs(d1) <= fabs(d2)) {
          val = d1 + g;
        } else {
          val = d2 + b;
        }
      }

      val = val > 1 ? val - 1 : val;
      val = val < 0 ? val + 1 : val;

      img[i] = val;
    }

  // Otherwise, we need to bin the colors
  } else {

    // The image
    img = mxGetPr(prhs[0]);

    /* Get the size of the image. */
    h = mxGetM(prhs[1]);
    w = mxGetN(prhs[1]);

    // Prepare the output, allocating the memory
    if ((plhs[0] = mxCreateDoubleMatrix(h, w, mxREAL)) == NULL) {
      mexErrMsgIdAndTxt("MATLAB:imsplitcolors:memoryAllocation",
        "Memory allocation failed !");
    }
    table = mxGetPr(plhs[0]);

    memcpy(table, mxGetPr(prhs[1]), nbins*sizeof(double)); 

    /// Need to to the 2D version if 2D table provided, sum only i+npix then
    /// If 3D, sum 1

    for (i=0; i < npix; i++) {
      dind = img[i]*nbins;
      indx = dind == nbins ? dind-1 : floor(dind);
      table[indx] += img[i + npix] * img [i + npix + npix];
    }
  }

  return;
}
