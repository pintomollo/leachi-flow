#include <math.h>
#include "mex.h"

// Bilinear interpolation, main interface
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // Declare variable
  int i, j, m, n, k, xf, yf, xc, yc, offset_img, offset_values;
  double ox, oy, ix, iy, dxf, dyf, dxc, dyc, x, y, z, nanval;
  mwSize w, h, c, dims[3];
  const mwSize *size;
  bool is_projective;
  double *proj, *tmp, *img, *values;

  // Check for proper number of input and output arguments
  if (nrhs != 5) {
    mexErrMsgIdAndTxt("MATLAB:imwarp:invalidInputs",
        "Not the correct number of input arguments (5 required) !");
  }

  // Ensure the types of the two first arrays at least
  if (!mxIsDouble(prhs[0])) {
    mexErrMsgIdAndTxt("MATLAB:bimwarpinvalidInputs",
        "Input arguments must be of type double.");
  }

  // The size of the image
  size = mxGetDimensions(prhs[0]);
  h = size[0];
  w = size[1];
  c = mxGetNumberOfElements(prhs[0]) / (h*w);
  img = mxGetPr(prhs[0]);

  // The project matrix
  if (mxGetNumberOfElements(prhs[1]) != 9) {
    mexErrMsgIdAndTxt("MATLAB:imwarp:invalidInputs",
      "The second argument must be a 3x3 transformation matrix");
  }
  proj = mxGetPr(prhs[1]);

  // The type of projection
  is_projective = (bool) mxGetScalar(prhs[2]);

  // The final size of the output
  tmp = mxGetPr(prhs[3]);
  n = (int) tmp[0];
  m = (int) tmp[1];

  // The offset
  tmp = mxGetPr(prhs[4]);
  ox = tmp[0];
  oy = tmp[1];

  // Prepare the output
  dims[0] = m;
  dims[1] = n;
  dims[2] = c;

  // Some temporary values
  offset_img = h*w;
  offset_values = m*n;

  // Create the array
  plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  values = mxGetPr(plhs[0]);

  // Precompute the value of NaN
  nanval = mxGetNaN();

  // The bilinear interpolation
  for (i=0; i < n; i++) {
    for (j=0; j < m; j++) {

      // The current indexes
      ix = i-ox+1;
      iy = j-oy+1;

      // Projecting the current pixel
      x = ix*proj[0] + iy*proj[1] + proj[2];
      y = ix*proj[3] + iy*proj[4] + proj[5];

      // Specifics to the projective method
      if (is_projective) {
        z = 1/(ix*proj[6] + iy*proj[7] + proj[8]);

        x *= z;
        y *= z;
      }

      // Adjust the indexes
      x -= 1;
      y -= 1;

      // Get the lower index
      xf = floor(x);
      yf = floor(y);

      // Its distance to the index
      dxf = x - xf;
      dyf = y - yf;

      // Avoid a singularity when the index are rounds, get the upper indexes
      if (dxf == 0) {
        xc = xf;
        dxc = 1;
      } else {
        xc = xf + 1;
        dxc = xc - x;
      }

      // Same for y
      if (dyf == 0) {
        yc = yf;
        dyc = 1;
      } else {
        yc = yf + 1;
        dyc = yc - y;
      }

      // Check whether all indexes are valid
      if (xf >= w || yf >= h || xc < 0 || yc < 0 || xc >= w || xf < 0 || yc >= h || yf < 0) {
        // All channels of the input image
        for (k=0; k < c; k++) {
          values[i*m + j + k*offset_values] = nanval;
        }

      // Compute the bilinear interpolation
      } else {
        // All channels of the input image
        for (k=0; k < c; k++) {
          values[i*m + j + k*offset_values] =
                      img[xf*h + yf + k*offset_img] * dxc * dyc +
                      img[xc*h + yf + k*offset_img] * dxf * dyc +
                      img[xf*h + yc + k*offset_img] * dxc * dyf +
                      img[xc*h + yc + k*offset_img] * dxf * dyf;
        }
      }
    }
  }

  return;
}
