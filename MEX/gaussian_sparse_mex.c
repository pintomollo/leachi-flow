#include <math.h>
#include "mex.h"

// A few useful definitions
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef __MAX__
#define __MAX__(A, B)     ((A)>=(B)? (A) : (B))
#endif
#ifndef __MIN__
#define __MIN__(A, B)     ((A)<=(B)? (A) : (B))
#endif
#ifndef __SQR__
#define __SQR__(A)        ((A) * (A))
#endif

/* A Matlab wrapper for the code of Anthony Gabrielson (see gaussian_smooth.c)*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  /* Declare a few variables. */
  mwSize m, n;
  mwSize nzmax, nzstep, nzfull;
  mwIndex *irs, *jcs, *irs2, *jcs2;
  int i, nelem, count, curr_col, curr_row;
  int r, rr, pr, irr, prs, dr, windowsize, center;
  double sigma, thresh, x, fx, sum=0.0, dot, sqrttwopi;
  double *rs, *rs2, *kernel;

  /* No flexibility here, we want both the image and sigma ! */
  if (nrhs != 3) {
    mexErrMsgIdAndTxt( "Bleachi:gaussian_sparse_mex:invalidNumInputs",
        "Three input arguments required.");
  } else if (!(mxIsSparse(prhs[0]) && mxIsDouble(prhs[0]))) {
    mexErrMsgIdAndTxt( "Bleachi:gaussian_sparse_mex:inputNotSparse",
        "Input arguments must be of type double sparse.");
  }

  /* Get sigma and thresh. */
  sigma = mxGetScalar(prhs[1]);
  thresh = mxGetScalar(prhs[2]);

  // Get the size of the input data
  m  = mxGetM(prhs[0]);
  n  = mxGetN(prhs[0]);

  // Get the sparse data
  rs  = mxGetPr(prhs[0]);
  irs = mxGetIr(prhs[0]);
  jcs = mxGetJc(prhs[0]);

  // The total number of data
  nelem = jcs[n];
  nzmax = mxGetNzmax(prhs[0]);

  // Data required in case new elements are created
  nzstep = (mwSize)__MAX__(ceil(((double)nzmax)*0.1),1);
  nzfull = m*n;

  // Prepare the output
  plhs[0] = mxCreateSparse(m,n,nzmax,false);
  rs2  = mxGetPr(plhs[0]);
  irs2 = mxGetIr(plhs[0]);
  jcs2 = mxGetJc(plhs[0]);

  // Create the kernel
  sqrttwopi = __SQR__(M_PI * 2);
  windowsize = 1 + 2 * ceil(2.5 * sigma);
  center = windowsize / 2;

  // Allocate the space
  if ((kernel = mxCalloc(windowsize, sizeof(double))) == NULL){
    mexErrMsgIdAndTxt( "Bleachi:gaussian_sparse_mex:memoryAllocation",
        "Memory allocation failed.");
  }

  // Compute the kernel
  for (i=0; i<windowsize; i++){
    x = ((double)i) - ((double)center);
    fx = exp(-0.5*x*x/(sigma*sigma)) / (sigma * sqrttwopi);

    kernel[i] = fx;
    sum += fx;
  }

  // And normalize it
  for(i=0;i<windowsize;i++) kernel[i] /= sum;

  // Now loop over the indexes to compute the gaussian blur along the columns
  curr_col = 0;
  count = 0;
  for (i = 0; i < nelem; i++) {

    // Find in which column we currently are, updating the Jc vector at the same time
    while (jcs[curr_col] <= i) {
      jcs2[curr_col] = count;
      curr_col++;
      pr = -1;
    }

    // Store the current row
    curr_row = irs[i];

    // Now loop over the previous positions, in case they were empty and now
    // need to be computed
    for (r=curr_row-center;r<=curr_row;r++){

      // If we have already computed index pr, we can skip it
      if (r<=pr){
        continue;
      }
      pr = r;

      // Compute the gaussian sum on the previous pixels
      dot = 0.0;
      sum = 0.0;
      for(irr=i;irr>=__MAX__((i-center),0);irr--){

        prs = irs[irr];

        // Need to make sure that we are not too far away or in another column
        if((irr >= ((int)jcs[curr_col-1])) && ((prs-r)>=(-center))){

          dr = prs-r;
          dot += rs[irr] * kernel[center+dr];
          sum += kernel[center+dr];
        } else {
          break;
        }
      }

      // Same for the pixels afterwards
      for(irr=i+1;irr<=__MIN__((i+center),nelem-1);irr++){
        prs = irs[irr];

        if((irr < jcs[curr_col]) && (prs-r<=center)){
          dr = prs-r;
          dot += rs[irr] * kernel[center+dr];
          sum += kernel[center+dr];
        } else {
          break;
        }
      }
      //sum = dot/sum;
      sum = dot;

      // If the value satisfies the threshold, store it !
      if (sum > thresh) {

        // Here we might need to increase the number of elements in the matrix
        if (count >= nzmax){
          nzmax += nzstep;
          nzmax = __MIN__(nzmax, nzfull);

          mxSetPr(plhs[0], mxRealloc(rs2, nzmax*sizeof(double)));
          mxSetIr(plhs[0], mxRealloc(irs2, nzmax*sizeof(mwIndex)));
          mxSetNzmax(plhs[0], nzmax);

          rs2  = mxGetPr(plhs[0]);
          irs2 = mxGetIr(plhs[0]);
        }

        // Actually store the values
        rs2[count]  = sum;
        irs2[count] = (mwIndex)r;
        count++;
      }
    }

    // Now there are two exceptions where we need to compute positions after our current one:
    //  - we are at the last position of the current row
    //  - there is a large gap (larger than the kernel) before the next value
    if ((i+1==jcs[curr_col]) || ((((int)irs[i+1])-curr_row) > center+1)){

      // Anyway, do not go outside of our column
      for (r=curr_row+1;r<=__MIN__(curr_row+center, m-1);r++){
        pr = r;

        // Compute the gaussian sum on the previous values (we know for sure there are none after)
        dot = 0.0;
        sum = 0.0;

        for(irr=i;irr>=__MAX__((i-center),0);irr--){
          prs = irs[irr];

          if((irr >= jcs[curr_col-1]) && ((prs-r)>=(-center))){
            dr = prs-r;
            dot += rs[irr] * kernel[center+dr];
            sum += kernel[center+dr];
          } else {
            break;
          }
        }
        //sum = dot/sum;
        sum = dot;

        // Maybe store it ?
        if (sum > thresh) {

          // Here we might need to increase the number of elements in the matrix
          if (count >= nzmax){
            nzmax += nzstep;
            nzmax = __MIN__(nzmax, nzfull);

            mxSetPr(plhs[0], mxRealloc(rs2, nzmax*sizeof(double)));
            mxSetIr(plhs[0], mxRealloc(irs2, nzmax*sizeof(mwIndex)));
            mxSetNzmax(plhs[0], nzmax);

            rs2  = mxGetPr(plhs[0]);
            irs2 = mxGetIr(plhs[0]);
          }

          rs2[count]  = sum;
          irs2[count] = (mwIndex)r;
          count++;
        }
      }
    }
  }

  // Needed to create a valid spare matrix by filling the last Jc values with the total
  for (i = curr_col; i <= n; i++) jcs2[i] = count;

  // Free the kernel
  mxFree(kernel);

  return;
}
