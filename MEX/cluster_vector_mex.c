#include <math.h>
#include "mex.h"

// Bilinear interpolation, main interface
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // Declare variable
  int i, j, counter = 0;
  mwSize n_vals;
  double *tmp, *vect_x, *vect_y, *cluster, *sizes, thresh, tmp_x, tmp_y, tmp_c;

  // Check for proper number of input and output arguments
  if (nrhs != 3) {
    mexErrMsgIdAndTxt("MATLAB:cluster_vector:invalidInputs",
        "Not the correct number of input arguments (3 are required) !");
  }

  // Get the size of the vector
  n_vals = mxGetNumberOfElements(prhs[0]);

  // And the vectors themselves
  vect_x = mxGetPr(prhs[0]);
  vect_y = mxGetPr(prhs[1]);

  // Make sure they are identical vectors
  if (n_vals != mxGetNumberOfElements(prhs[1])) {
    mexErrMsgIdAndTxt("MATLAB:cluster_vector:invalidInputs",
        "The two first inputs shoulf have the same number of elements !");
  }

  // Get the threshold value for clustering
  thresh = mxGetScalar(prhs[2]);

  // Prepare the output, allocating the memory
  if ((plhs[0] = mxCreateDoubleMatrix(n_vals, 1, mxREAL)) == NULL) {
    mexErrMsgIdAndTxt("MATLAB:cluster_vector:memoryAllocation",
      "Memory allocation failed !");
  }
  cluster = mxGetPr(plhs[0]);

  // Loop over all vectors
  for (i = 0; i < n_vals; i++) {

    // Create a new cluster if need be
    if (cluster[i] == 0) {
      cluster[i] = ++counter;
    }

    // Store the current vector
    tmp_x = vect_x[i];
    tmp_y = vect_y[i];
    tmp_c = cluster[i];

    // Now look through the other vectors
    for (j = i+1; j < n_vals; j++) {

      // Consider only the unassigned ones
      if (cluster[j] == 0) {

        // If they meet the threshold, add them to our cluster
        if (fabs(tmp_x - vect_x[j]) <= thresh && fabs(tmp_y - vect_y[j]) <= thresh) {
          cluster[j] = tmp_c;
        }
      }
    }
  }

  // Measure un addition the size of the clusters
  if (nlhs > 1) {

    // Prepare the second output, allocating the memory
    if ((plhs[1] = mxCreateDoubleMatrix(counter, 1, mxREAL)) == NULL) {
      mexErrMsgIdAndTxt("MATLAB:cluster_vector:memoryAllocation",
        "Memory allocation failed !");
    }
    sizes = mxGetPr(plhs[1]);

    // Loop over all clusters and count them
    for (i = 0; i < n_vals; i++) {
      sizes[(int)(cluster[i]-1)]++;
    }
  }

  return;
}
