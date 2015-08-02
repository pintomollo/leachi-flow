#include <math.h>
#include "mex.h"

#ifndef MIN
#define MIN(a,b) ((a) > (b) ? (b) : (a))
#endif

// Bilinear interpolation, main interface
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // Declare variable
  int i;
  mwSize npix, ndim;
  const mwSize *dims;
  double *hsv;
  double r, g, b, h, s, v, k, p, norm, norm2, maxval;
  unsigned char *rgb;
  bool rgb2hsv = true;

  // Check for proper number of input and output arguments
  if (nrhs != 1) {
    mexErrMsgIdAndTxt("MATLAB:imsplitcolors:invalidInputs",
        "Not the correct number of input arguments (1 required) !");
  }

  // Ensure the types of the two first arrays at least
  if (mxIsDouble(prhs[0])) {
    rgb2hsv = false;
    hsv = mxGetPr(prhs[0]);
  } else if (mxGetClassID(prhs[0]) == mxUINT8_CLASS) {
    rgb = (unsigned char *) mxGetData(prhs[0]);
  } else {
    mexErrMsgIdAndTxt("MATLAB:imsplitcolors:invalidInputs",
        "Input arguments must be of type double or uint8.");
  }

  // Get the number of pixels and the image
  npix = mxGetNumberOfElements(prhs[0]) / 3;
  ndim = mxGetNumberOfDimensions(prhs[0]);
  dims = mxGetDimensions(prhs[0]);

  if (rgb2hsv) {

    if ((plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL)) == NULL) {
      mexErrMsgIdAndTxt("MATLAB:imsplitcolors:memoryAllocation",
        "Memory allocation failed !");
    }
    hsv = mxGetPr(plhs[0]);
    norm = 1 / 255.0;
    norm2 = 1 / 6.0;

    for (i=0; i < npix; i++) {
      r = ((double) rgb[i]) * norm;
      g = ((double) rgb[i + npix]) * norm;
      b = ((double) rgb[i + npix + npix]) * norm;

      h = 0.0;

      if (r > g) {
        if (r > b) {
          v = r;
          s = v - MIN(g, b);

          if (s != 0) {
            h = (g - b) / s;
          }
        } else {
          v = b;
          s = v - MIN(r, g);

          if (s != 0) {
            h = 4 + (r - g) / s;
          }
        }
      } else {
        if (g > b) {
          v = g;
          s = v - MIN(r, b);

          if (s != 0) {
            h = 2 + (b - r) / s;
          }
        } else {
          v = b;
          s = v - MIN(r, g);

          if (s != 0) {
            h = 4 + (r - g) / s;
          }
        }
      }

      if (v > 0) {
        h = h * norm2;

        if (h < 0) {
          h += 1;
        }

        s /= v;
      }

      hsv[i] = h;
      hsv[i + npix] = s;
      hsv[i + npix + npix] = v;
    }
  } else {

    if ((plhs[0] = mxCreateNumericArray(ndim, dims, mxUINT8_CLASS, mxREAL)) == NULL) {
      mexErrMsgIdAndTxt("MATLAB:imsplitcolors:memoryAllocation",
        "Memory allocation failed !");
    }
    rgb = (unsigned char *) mxGetData(plhs[0]);

    maxval = 0;

    for (i=0; i < npix; i++) {
      h = 6 * hsv[i];
      s = hsv[i + npix];
      v = hsv[i + npix + npix] * 255;

      k = floor(h);
      p = h-k;

      switch ((int) k) {
        case 1 :
          r = 1 - s*p;
          g = 1;
          b = 1 - s;
          break;
        case 2 :
          r = 1 - s;
          g = 1;
          b = 1 - (s * (1-p));
          break;
        case 3 :
          r = 1 - s;
          g = 1 - s*p;
          b = 1;
          break;
        case 4 :
          r = 1 - (s * (1-p));
          g = 1 - s;
          b = 1;
          break;
        case 5 :
          r = 1;
          g = 1 - s;
          b = 1 - s*p;
          break;
        default :
          r = 1;
          g = 1 - (s * (1-p));
          b = 1 - s;
      }

      rgb[i] = (unsigned char) (r * v);
      rgb[i + npix] = (unsigned char) (g * v);
      rgb[i + npix + npix] = (unsigned char) (b * v);
    }
  }

  return;
}
