#include "mex.h"
#include "svdcmp.cpp"
#include <unordered_map>
#include <math.h>
#include <string.h>

// Some useful macros
#ifndef __MAX__
#define __MAX__(A, B)     ((A)>=(B)? (A) : (B))
#endif
#ifndef __MIN__
#define __MIN__(A, B)     ((A)<=(B)? (A) : (B))
#endif
#ifndef __STEP__
#define __STEP__  1000
#endif

// Functions to read and write 3D vectors
template<class T> inline void get_shifted_vector(T *x, T *y, T *m) { x[0]=y[0]-m[0]; x[1]=y[1]-m[1]; x[2]=y[2]-m[2]; }
template<class T, class U> inline void get_cluster_id(T *x, U *y, U s) { x[0]=ceil(y[0]*s); x[1]=ceil(y[1]*s); x[2]=ceil(y[2]*s); }
template<class T, class U> inline void set_cluster_pos(T *x, U *y, T s) { x[0]=((double)y[0])*s; x[1]=((double)y[1])*s; x[2]=((double)y[2])*s; }

// Mathematical functions for 3D vectors
template<class T> inline void diff3(T *s,T *a,T *b) { s[0]=a[0]-b[0]; s[1]=a[1]-b[1]; s[2]=a[2]-b[2]; }
template<class T> inline void plus3(T *s,T *a,T *b) { s[0]=a[0]+b[0]; s[1]=a[1]+b[1]; s[2]=a[2]+b[2]; }
template<class T> inline void prod3(T *s,T *a,T b) { s[0]=a[0]*b; s[1]=a[1]*b; s[2]=a[2]*b; }

// And some linear algebra in 3D
template<class T> inline T norm3(T *a) { return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]); }
template<class T> inline T dot3(T *a,T *b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
template<class T> inline void cross3(T *s,T *a,T *b) { s[0]=a[1]*b[2]-a[2]*b[1]; s[1]=a[2]*b[0]-a[0]*b[2]; s[2]=a[0]*b[1]-a[1]*b[0]; }

// Addition in 4D
template<class T> inline void plus4(T *s,T *a,T *b) { s[0]=a[0]+b[0]; s[1]=a[1]+b[1]; s[2]=a[2]+b[2]; s[3]=a[3]+b[3]; }

// And some functions on matrices
// Loading the symmetric 3x3 matrix from the quadric
template<class T> inline void set33(T **x, T *Q) { 
  x[1][1] = Q[0];
  x[2][1] = x[1][2] = Q[1];
  x[2][2] = Q[2];
  x[3][1] = x[1][3] = Q[3];
  x[3][2] = x[2][3] = Q[4];
  x[3][3] = Q[5];
}

// Multiplication of a 3x3 matrix with a 3x1 vector
template<class T> inline void prod33(T *s, T **A, T *x) { 
  s[0] = A[1][1]*x[0] + A[1][2]*x[1] + A[1][3]*x[2];
  s[1] = A[2][1]*x[0] + A[2][2]*x[1] + A[2][3]*x[2];
  s[2] = A[3][1]*x[0] + A[3][2]*x[1] + A[3][3]*x[2];
}

// The same but with the transposed matrix
template<class T> inline void prodT33(T *s, T **A, T *x) { 
  s[0] = A[1][1]*x[0] + A[2][1]*x[1] + A[3][1]*x[2];
  s[1] = A[1][2]*x[0] + A[2][2]*x[1] + A[3][2]*x[2];
  s[2] = A[1][3]*x[0] + A[2][3]*x[1] + A[3][3]*x[2];
}

// Computing and extraction the unique part of a symmetrix n'*n matrix
template<class T> inline void upper_prod(T *Q,T *a) {
  Q[0]=a[0]*a[0];
  Q[1]=a[0]*a[1];
  Q[2]=a[1]*a[1];
  Q[3]=a[0]*a[2];
  Q[4]=a[1]*a[2];
  Q[5]=a[2]*a[2];
  Q[6]=a[0]*a[3];
  Q[7]=a[1]*a[3];
  Q[8]=a[2]*a[3];
}

// Storing the weighted sum of the two quadrics
template<class T> inline void store_quadrics(T *Q,T *S,T *B,T l) {
  Q[0] += S[0]*l + B[0]*(1-l);
  Q[1] += S[1]*l + B[1]*(1-l);
  Q[2] += S[2]*l + B[2]*(1-l);
  Q[3] += S[3]*l + B[3]*(1-l);
  Q[4] += S[4]*l + B[4]*(1-l);
  Q[5] += S[5]*l + B[5]*(1-l);
  Q[6] += S[6]*l + B[6]*(1-l);
  Q[7] += S[7]*l + B[7]*(1-l);
  Q[8] += S[8]*l + B[8]*(1-l);
}

// Out-of-Core simplification algorithm from :
//  [1] Lindstrom, P. (2000). Out-of-core simplification of large polygonal models. Proceedings of the 27th Annual Conference on Computer Graphics and Interactive Techniques - SIGGRAPH ’00, 259–262. http://doi.org/10.1145/344779.344912
//  [2] Lindstrom, P., & Silva, C. T. (2001). A memory insensitive technique for large model simplification. Proceedings Visualization, 2001. VIS ’01. http://doi.org/10.1109/VISUAL.2001.964502
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwSize nfaces, nvertices, nnewf, nnewv;
  int sub1, sub2, sub3, i, id, clust_id=0, curr_f=0;
  int c1[3], c2[3], c3[3], sizes[3], pad_size;
  double step = -100, istep, lambda = 0.5, thresh = 1e-3, val1, val2, val3;
  double v1[3], v2[3], v3[3], n[4], Qs[9], tmp1[4], tmp2[4], tmp3[4];
  double e1[3], e2[3], e3[3], m1[4], m2[4], m3[4], Qb1[9], Qb2[9], Qb3[9];
  double *face, *vertex, *tf, *tv, *tq, *tnf, *new_vertices, *new_faces, *quadrics;
  double **U, E[4], **V, b[3], *tmp;
  std::unordered_map<int,int> mymap;

  // Check for proper number of input and output arguments
  if (nrhs < 2) {
    mexErrMsgIdAndTxt("MATLAB:oocs_mex:invalidInputs",
        "Not the correct number of input arguments (2 minimum required) !");

  } else if (mxGetM(prhs[0]) != 3 || mxGetM(prhs[1]) != 3) {

    mexErrMsgIdAndTxt("MATLAB:oocs_mex:invalidInputs",
        "Not the correct size for the input arguments (3xN required) !");

  // Get the tables and the step
  } else if (nrhs == 2) {
    face = mxGetPr(prhs[0]);
    vertex = mxGetPr(prhs[1]);
  } else {
    face = mxGetPr(prhs[0]);
    vertex = mxGetPr(prhs[1]);
    step = mxGetScalar(prhs[2]);
  }

  // Get the number of points
  nfaces = mxGetN(prhs[0]);
  nvertices = mxGetN(prhs[1]);

  // The size of the increments
  nnewf = nnewv = __STEP__;
  pad_size = 9*__STEP__*sizeof(double);

  // Prepare the outputs, allocating the memory
  if ((plhs[0] = mxCreateDoubleMatrix(3, nnewf, mxREAL)) == NULL) {
    mexErrMsgIdAndTxt("MATLAB:oocs_mex:memoryAllocation",
      "Memory allocation failed !");
  }
  new_faces = mxGetPr(plhs[0]);

  if ((plhs[1] = mxCreateDoubleMatrix(3, nnewv, mxREAL)) == NULL) {
    mexErrMsgIdAndTxt("MATLAB:oocs_mex:memoryAllocation",
      "Memory allocation failed !");
  }
  new_vertices = mxGetPr(plhs[1]);

  if ((quadrics = (double*) mxCalloc(9*nnewv, sizeof(double))) == NULL) {
    mexErrMsgIdAndTxt("MATLAB:oocs_mex:memoryAllocation",
      "Memory allocation failed !");
  }

  // Initialize the range
  double mins[3] = {vertex[0], vertex[1], vertex[2]};
  double maxs[3] = {vertex[0], vertex[1], vertex[2]};

  // Compute the range of the coordinates
  tv = vertex + 3;
  for (i=0; i<nvertices-1; i++) {
    mins[0] = __MIN__(mins[0], tv[0]);
    mins[1] = __MIN__(mins[1], tv[1]);
    mins[2] = __MIN__(mins[2], tv[2]);

    maxs[0] = __MAX__(maxs[0], tv[0]);
    maxs[1] = __MAX__(maxs[1], tv[1]);
    maxs[2] = __MAX__(maxs[2], tv[2]);

    tv += 3;
  }

  // Compute the step if given as a negative input
  if (step < 0) {
    step = -__MAX__(__MAX__(maxs[0]-mins[0], maxs[1]-mins[1]), maxs[2]-maxs[2])/step;
  }
  istep = 1 / step;

  // Shift the whole by half a space
  mins[0] -= step/2;
  mins[1] -= step/2;
  mins[2] -= step/2;

  // Get the dimension of the mapping cube by getting the indexes of the max point
  get_shifted_vector(maxs, maxs, mins);
  get_cluster_id(sizes, maxs, istep);
  sizes[0] = sizes[1]*sizes[2];
  sizes[1] = sizes[2];

  // Get temporary pointers
  tf = face;
  tnf = new_faces;
  tv = new_vertices;

  // Loop over all faces
  for (i=0; i<nfaces; i++) {

    // Get the vertices of the current face
    get_shifted_vector(v1, vertex + 3*((int)tf[0]-1), mins);
    get_shifted_vector(v2, vertex + 3*((int)tf[1]-1), mins);
    get_shifted_vector(v3, vertex + 3*((int)tf[2]-1), mins);

    // Get the corresponding indexes
    get_cluster_id(c1, v1, istep);
    get_cluster_id(c2, v2, istep);
    get_cluster_id(c3, v3, istep);

    // And the corresponding single index
    sub1 = c1[0]*sizes[0] + c1[1]*sizes[1] + c1[2];
    sub2 = c2[0]*sizes[0] + c2[1]*sizes[1] + c2[2];
    sub3 = c3[0]*sizes[0] + c3[1]*sizes[1] + c3[2];

    // Process only the faces spanning three different clusters
    if (sub1 != sub2 && sub2 != sub3 && sub3 != sub1) {

      // Check if the cluster already exists
      id = mymap[sub1];

      // If not, create one, making sure there is enough space to store it
      if (id == 0) {
        if (clust_id >= nnewv) {
          nnewv += __STEP__;

          if((tmp = (double*) mxRealloc(new_vertices, 3*nnewv*sizeof(double))) == NULL) {
            mxFree(quadrics);
            mexErrMsgIdAndTxt("MATLAB:oocs_mex:memoryReallocation",
              "Memory reallocation failed !");
          }
          new_vertices = tmp;
          if ((tmp = (double*) mxRealloc(quadrics, 9*nnewv*sizeof(double))) == NULL) {
            mxFree(quadrics);
            mexErrMsgIdAndTxt("MATLAB:oocs_mex:memoryReallocation",
              "Memory reallocation failed !");
          }
          quadrics = tmp;
          memset(&quadrics[(nnewv-__STEP__)*9], 0, pad_size);
          tv = new_vertices + 3*clust_id;
        }
        set_cluster_pos(tv, c1, step);
        tv += 3;

        // Store the new cluster id
        clust_id++;
        mymap[sub1] = clust_id;
        sub1 = clust_id;

      // Otherwise, retrieve the id
      } else {
        sub1 = id;
      }

      // Same for the second cluster
      id = mymap[sub2];
      if (id == 0) {
        if (clust_id >= nnewv) {
          nnewv += __STEP__;

          if((tmp = (double*) mxRealloc(new_vertices, 3*nnewv*sizeof(double))) == NULL) {
            mxFree(quadrics);
            mexErrMsgIdAndTxt("MATLAB:oocs_mex:memoryReallocation",
              "Memory reallocation failed !");
          }
          new_vertices = tmp;
          if ((tmp = (double*) mxRealloc(quadrics, 9*nnewv*sizeof(double))) == NULL) {
            mxFree(quadrics);
            mexErrMsgIdAndTxt("MATLAB:oocs_mex:memoryReallocation",
              "Memory reallocation failed !");
          }
          quadrics = tmp;
          memset(&quadrics[(nnewv-__STEP__)*9], 0, pad_size);
          tv = new_vertices + 3*clust_id;
        }
        set_cluster_pos(tv, c2, step);
        tv += 3;

        clust_id++;
        mymap[sub2] = clust_id;
        sub2 = clust_id;
      } else {
        sub2 = id;
      }

      // And the third one
      id = mymap[sub3];
      if (id == 0) {
        if (clust_id >= nnewv) {
          nnewv += __STEP__;

          if((tmp = (double*) mxRealloc(new_vertices, 3*nnewv*sizeof(double))) == NULL) {
            mxFree(quadrics);
            mexErrMsgIdAndTxt("MATLAB:oocs_mex:memoryReallocation",
              "Memory reallocation failed !");
          }
          new_vertices = tmp;
          if ((tmp = (double*) mxRealloc(quadrics, 9*nnewv*sizeof(double))) == NULL) {
            mxFree(quadrics);
            mexErrMsgIdAndTxt("MATLAB:oocs_mex:memoryReallocation",
              "Memory reallocation failed !");
          }
          quadrics = tmp;
          memset(&quadrics[(nnewv-__STEP__)*9], 0, pad_size);
          tv = new_vertices + 3*clust_id;
        }
        set_cluster_pos(tv, c3, step);
        tv += 3;

        clust_id++;
        mymap[sub3] = clust_id;
        sub3 = clust_id;
      } else {
        sub3 = id;
      }

      // And now store the corresponding face, making sure there is enough space
      if (curr_f >= nnewf) {
        nnewf += __STEP__;

        if ((tmp = (double*) mxRealloc(new_faces, 3*nnewf*sizeof(double))) == NULL) {
          mxFree(quadrics);
          mexErrMsgIdAndTxt("MATLAB:oocs_mex:memoryReallocation",
            "Memory reallocation failed !");
        }
        new_faces = tmp;

        tnf = new_faces + 3*curr_f;
      }
      tnf[0] = sub1;
      tnf[1] = sub2;
      tnf[2] = sub3;
      tnf += 3;

      curr_f++;

      // Compute the first 4D quadric [1]:
      //   3D sum of cross-products:  x1 X x2 + x2 X x3 + x3 X x1
      //   1D scalar triple product: -[x1, x2, x3]
      cross3(tmp1, v1, v2);
      cross3(tmp2, v2, v3);
      cross3(tmp3, v3, v1);

      plus3(n, tmp1, tmp2);
      plus3(n, n, tmp3);

      n[3] = -dot3(v1, tmp2);

      // Extract the corresponding quadric
      upper_prod(Qs, n);

      // Then the second ones (the edge ones, [2]):
      //    3D normals to the edges: m = ||e||(e X n)
      //    1D dot product         : -0.5 (x1 + x2)' * m
      //
      //  where n is the normal to the surface, e one of the edge vectors (3 of them),
      //  x1 and x2 the corresponding ends of e, and the connected m are added

      // The normal to the face
      prod3(n, n, 1/norm3(n));

      // The edges
      diff3(e1, v2, v1);
      diff3(e2, v3, v2);
      diff3(e3, v1, v3);

      val1 = norm3(e1);
      val2 = norm3(e2);
      val3 = norm3(e3);

      // The normals to the edges
      cross3(tmp1, e1, n);
      cross3(tmp2, e2, n);
      cross3(tmp3, e3, n);

      // The scaled ms
      prod3(m1, tmp1, val1);
      prod3(m2, tmp2, val2);
      prod3(m3, tmp3, val3);

      // The middle points
      plus3(tmp1, v2, v1);
      plus3(tmp2, v3, v2);
      plus3(tmp3, v1, v3);

      // The dot product
      m1[3] = -0.5*dot3(tmp1, m1);
      m2[3] = -0.5*dot3(tmp2, m2);
      m3[3] = -0.5*dot3(tmp3, m3);

      // Add the ms together
      plus4(tmp1, m1, m3);
      plus4(tmp2, m2, m1);
      plus4(tmp3, m3, m2);

      // Compute the quadrics
      upper_prod(Qb1, tmp1);
      upper_prod(Qb2, tmp2);
      upper_prod(Qb3, tmp3);

      // And store them all, weighted as advised in OoCSx
      store_quadrics(&quadrics[(sub1-1)*9], Qs, Qb1, lambda);
      store_quadrics(&quadrics[(sub2-1)*9], Qs, Qb2, lambda);
      store_quadrics(&quadrics[(sub3-1)*9], Qs, Qb3, lambda);
    }

    tf += 3;
  }

  // Adjust the various matrices to the actually required spaces
  if ((tmp = (double*) mxRealloc(new_faces, 3*curr_f*sizeof(double))) == NULL) {
    mxFree(quadrics);
    mexErrMsgIdAndTxt("MATLAB:oocs_mex:memoryReallocation",
      "Memory reallocation failed !");
  }
  new_faces = tmp;
  if ((tmp = (double*) mxRealloc(new_vertices, 3*clust_id*sizeof(double))) == NULL) {
    mxFree(quadrics);
    mexErrMsgIdAndTxt("MATLAB:oocs_mex:memoryReallocation",
      "Memory reallocation failed !");
  }
  new_vertices = tmp;
  if ((tmp = (double*) mxRealloc(quadrics, 9*clust_id*sizeof(double))) == NULL) {
    mxFree(quadrics);
    mexErrMsgIdAndTxt("MATLAB:oocs_mex:memoryReallocation",
      "Memory reallocation failed !");
  }
  quadrics = tmp;

  // And adapt the Matlab structures
  mxSetPr(plhs[0], new_faces);
  mxSetN(plhs[0], curr_f);
  mxSetPr(plhs[1], new_vertices);
  mxSetN(plhs[1], clust_id);

  // Somehow the svdcmp uses matrices of size N+1, so need to allocate 4x4 ones !
  U = new double*[4];
  V = new double*[4];
  for(i = 0; i < 4; ++i) {
    U[i] = new double[4];
    V[i] = new double[4];
  }

  // Prepare temporary pointers
  tq = quadrics;
  tv = new_vertices;

  // Loop over all vertices
  for (i=0; i<clust_id; i++) {

    // The algorithm goes as follow:
    //   the quadric stores the solution of the linear system A*x=b for each vertex
    //   to solve the inversion, we decompose A=U*E*V' using SVD
    //   the centered vertex is then x = v + V*E*U'*(b-A*v)
    //   where v is the center of the current cluster

    // Get the current center of the cluster
    v1[0] = tv[0] - step/2;
    v1[1] = tv[1] - step/2;
    v1[2] = tv[2] - step/2;

    // The b part of the quadric
    b[0] = -tq[6];
    b[1] = -tq[7];
    b[2] = -tq[8];

    // And the symmetric A matrix from the quadric
    set33(U, tq);

    // Compute A*v already, as svdcmp will overwrite A
    prod33(tmp1, U, v1);

    // Decompose A
    svdcmp(U,3,3,E,V);

    // The robustness fix proposed by Lindstrom
    E[0] = 1/E[1];
    E[1] = E[2]*E[0] > thresh ? 1/E[2] : 0.0;
    E[2] = E[3]*E[0] > thresh ? 1/E[3] : 0.0;

    // b-A*v
    diff3(tmp2, b, tmp1);

    // U'*(b-A*v)
    prodT33(tmp3, U, tmp2);

    // E is diagonal, so: E*U'*(b-A*v)
    E[0] *= tmp3[0];
    E[1] *= tmp3[1];
    E[2] *= tmp3[2];

    // V*E*U'*(b-A*v)
    prod33(v2, V, E);

    // The shifted x !
    tv[0] = v1[0] + v2[0] + mins[0];
    tv[1] = v1[1] + v2[1] + mins[1];
    tv[2] = v1[2] + v2[2] + mins[2];

    // Move the pointers
    tv += 3;
    tq += 9;
  }

  // Free the unused memory
  mxFree(quadrics);

  return;
}
