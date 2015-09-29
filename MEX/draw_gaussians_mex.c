#include <math.h>
#include "mex.h"

#ifndef __SQR__
#define __SQR__(A)        ((A) * (A))
#endif

/// Approximation of the exponential function from :
/// N. N. Schraudolph. A Fast, Compact Approximation of the Exponential Function. 
/// Neural Computation, 11(4):853–862, 1999.

/// 2x to 9x faster than exp(x)!
/// Can be off by about +-4% in the range -100 to 100.
double fast_exp(double y) {
  double d;
  y = (y < -700) ? -700 : (y > 700 ? 700 : y);
  *((int*)(&d) + 0) = 0;
  *((int*)(&d) + 1) = (int)(1512775 * y + 1072632447);
  return d;
}

// Drawing gaussian spots on a discretized image
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  mwSize ndim, n, nprops;
  int i, j, c, w, h;
  double *tmp, *cells, *img, *x, *y, *sigma, *ampl;
  bool free_memory = false;

  // Check for proper number of input and output arguments
  if (nrhs != 2) {
    mexErrMsgIdAndTxt("MATLAB:draw_gaussians:invalidInputs",
        "Not the proper number of input arguments (2 are expected) !");
  }

  // Accept only 2D data
  ndim = mxGetNumberOfElements(prhs[0]);
  if (ndim != 2) {
    mexErrMsgIdAndTxt("MATLAB:draw_gaussians:invalidDimensions",
        "Currently only 2D data are accepted !");

  } else {
    // Retrieve the dimensions of the image
    tmp = mxGetPr(prhs[0]);
    h = (int)tmp[0];
    w = (int)tmp[1];
  }

  // Prepare the output
  plhs[0] = mxCreateDoubleMatrix(h, w, mxREAL);
  img = mxGetPr(plhs[0]);

  // Get the dimensions of the gaussians
  cells = mxGetPr(prhs[1]);
  n = mxGetM(prhs[1]);
  nprops = mxGetN(prhs[1]) - ndim;

  // We might not need to do anything at all...
  if (n > 0) {
    // We need at least their diameters
    if (nprops < 1) {
      mexErrMsgIdAndTxt("MATLAB:draw_gaussians:invalidGaussians",
          "Not enough parameters of the Gaussian spots are provided !");
    } else if (nprops == 1) {
      // And thus need to allocate memory for the temporary amplitudes
      if ((ampl = mxCalloc(n, sizeof(double))) == NULL) {
        mexErrMsgIdAndTxt("MATLAB:draw_gaussians:memoryAllocation", 
          "Memory allocation failed !");
      }

      // Set all amplitudes to 1
      for (c=0; c < n; c++) {
        ampl[c] = 1.0;
      }

      // Need to free this memory at the end
      free_memory = true;

    // Just access the amplitude then
    } else {
      ampl = cells + 3*n;
    }

    // Get each column separately
    x = cells;
    y = x + n;
    tmp = y + n;

    // And allocate memory for the precomputed sigma
    if ((sigma = mxCalloc(n, sizeof(double))) == NULL) {
      mexErrMsgIdAndTxt("MATLAB:draw_gaussians:memoryAllocation", 
        "Memory allocation failed !");
    }

    // Precompute all sigmas
    for (c=0; c < n; c++) {
      sigma[c] = -1.0 / (2.0 * __SQR__(tmp[c]));
    }

    // Now go through every pixel of the image
    for (i=1; i <= w; i++) {
      for (j=1; j <= h; j++) {
        for (c=0; c < n; c++) {
          img[j-1] += ampl[c] * fast_exp(sigma[c] * (__SQR__((double)i - x[c]) + __SQR__((double)j - y[c])));
        }
      }

      img += h;
    }
  }

  // Free the allocated memory
  if (free_memory) {
    mxFree(ampl);
  }
  mxFree(sigma);

  return;
}
