//------------------------------------------------------------------
//
// Copyright (c) 2012 Jeffrey Byrne
// $Id: nesting_distance.cpp 158 2013-04-11 18:50:04Z jebyrne $
//
//------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h> 
#include <mex.h>
#include <matrix.h>
#include <math.h>
#include <vector>
#include <algorithm>

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))


//-----------------------------------------------------
float euclideanDistance(const float *P, const float *Q, int n)
{
  int k;
  float d=0, dpq;
  
  for (k=0; k<n; k++) {
    dpq = P[k]-Q[k];
    d += dpq*dpq;
  }
  return(d);
}

//-----------------------------------------------------
float mahalanobisDistance(const float *P, const float *Q, const float s, int n)
{
  int k;
  float d=0, dpq;
  
  for (k=0; k<n; k++) {
    dpq = P[k]-Q[k];
    d += dpq*s*dpq;
  }
  return(d);
}

//-----------------------------------------------------
float L1Distance(const float *P, const float *Q, int n)
{
  int k;
  float d, dpq;
  
  d = 0;
  for (k=0; k<n; k++) {
    d += fabs(P[k]-Q[k]);
  }
  return(d);
}



//-----------------------------------------------------
// Mexify!
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int i,j,k,n,m;
  int n_dim, n_lobes, n_scales, n_bands;
  const mxArray *mxP,*mxQ;
  mxArray *mxD;
  float *D;
  const float *P, *Q;
  int do_debug = 0;
  std::vector<float> v,p,b,S,v_trim;  
  float dpq,d,mu;
  int n_median;
  int k_outlier, k_inlier;
  float t_outlier,p_outlier,p_inlier;
  
  // Input check
  if (nrhs < 8) {
    mexErrMsgTxt("invalid mex input\n");
  }    
  
  // Inputs
  mxP = prhs[0];  
  mxQ = prhs[1];
  n_bands = (int) mxGetScalar(prhs[2]);
  n_lobes = (int) mxGetScalar(prhs[3]);
  n_scales = (int) mxGetScalar(prhs[4]);
  p_inlier = (float) mxGetScalar(prhs[5]);
  p_outlier = (float) mxGetScalar(prhs[6]);
  t_outlier = (float) mxGetScalar(prhs[7]);

  
  // Derived inputs
  P = (const float *)mxGetPr(mxP);
  Q = (const float *)mxGetPr(mxQ);
  n = mxGetN(mxP);  // points
  m = mxGetN(mxQ);  // points
  n_dim = mxGetM(mxP);  // dimension
  k_outlier = (int)(p_outlier*n_dim);
  k_inlier = (int)(p_inlier*n_dim);
  
  // Input verification
  if (mxGetM(mxP) != mxGetM(mxQ)) {
    mexErrMsgTxt("invalid mex input - dimensionality mismatch\n");
  }
  if (n_dim != (n_bands*n_lobes*n_scales)) {
    mexErrMsgTxt("invalid mex input - dimensionality mismatch\n");
  }
  
  // Outputs
  mxD = mxCreateNumericMatrix(n,m,mxSINGLE_CLASS,mxREAL);
  D = (float *) mxGetPr(mxD);  
    
  // Buffers
  v.resize(n_dim);  
  p.resize(n_dim);  
  b.resize(n_dim);  
  
  // Nesting distance
  for (j=0; j<m; j++) {          
    for (i=0; i<n; i++) {    // column-major

      
#if 0
      // Euclidean distance vector
      for (k=0; k<n_dim; k++) {
        dpq = P[i*n_dim + k]-Q[j*n_dim + k];
        v[k] = dpq*dpq;
        p[k] = v[k];
        b[k] = (float)((P[i*n_dim + k] > 0) != (Q[j*n_dim + k] > 0));  // hamming
      }
      std::nth_element(v.begin(), v.begin()+k_inlier, v.end());
      D[j*n+i] = 0;
      for (k=0; k<n_dim; k++) {
        if (p[k] <= v[k_inlier]) {
          D[j*n+i] += (float)b[k];
        }
      }
#endif

      
#if 0
      // Scaled pooling
      for (int k_band=0; k_band<n_bands; k_band++) {
        for (int k_lobe=0; k_lobe<n_lobes; k_lobe++) {
          for (int k_scale=0; k_scale<n_scales; k_scale++) {
            int k = (k_scale)*n_bands*n_lobes + k_lobe*n_bands + k_band;
            int kp = MIN((k_scale+1),(n_scales-1))*n_bands*n_lobes + k_lobe*n_bands + k_band;
            int km = MAX((k_scale-1),0)*n_bands*n_lobes + k_lobe*n_bands + k_band;
            
            dpq = (P[i*n_dim+k]-Q[j*n_dim+k])*(P[i*n_dim+k]-Q[j*n_dim+k]);
            float dpqp = (P[i*n_dim+k]-Q[j*n_dim+kp])*(P[i*n_dim+k]-Q[j*n_dim+kp]);
            float dpqm = (P[i*n_dim+k]-Q[j*n_dim+km])*(P[i*n_dim+k]-Q[j*n_dim+km]);
            v[k] = MIN(MIN(dpq,dpqp),dpqm);
          }          
        }
      }

    // Outlier?
    d = 0;
      if (t_outlier > 0) { // speedup
        std::nth_element(v.begin(), v.begin()+k_outlier, v.end());
        for (k=0; k<k_outlier; k++) {
          d += (float)v[k];
        }
      }
      // Inlier distance
      if (d <= t_outlier) {        
        std::nth_element(v.begin(), v.begin()+k_inlier, v.end());
        D[j*n+i] = 0;
        for (k=0; k<k_inlier; k++) {
          D[j*n+i] += (float)v[k];
        }
      }
      else {
        D[j*n+i] = -1;  // invalid distance
      }
#endif
      
#if 1
      // Euclidean distance vector
      for (k=0; k<n_dim; k++) {
        dpq = P[i*n_dim + k]-Q[j*n_dim + k];
        v[k] = dpq*dpq;
      }
      // Outlier?      
      d = 0;
//       if (t_outlier > 0) { // speedup
//         std::nth_element(v.begin(), v.begin()+k_outlier, v.end());
//         for (k=0; k<k_outlier; k++) {
//           d += (float)v[k];
//         }
//       }
      // Inlier distance
//      if (d <= t_outlier) {        
      if (1) {        
        std::nth_element(v.begin(), v.begin()+k_inlier, v.end());
        D[j*n+i] = 0;
        for (k=0; k<k_inlier; k++) {
          D[j*n+i] += (float)v[k];
        }
      }
      else {
        D[j*n+i] = -1;  // invalid distance
      }
#endif

      
#if 0
      // Euclidean distance       
      d=0;
      for (k=0; k<n_dim; k++) {
        dpq = P[i*n_dim + k]-Q[j*n_dim + k];
        v[k] = dpq*s_mahalanobis*dpq;
        d += v[k];
      }
      // Truncated and Trimmed
      if (d <= t_truncate) {
        int k_trim = (int)(((float)n_dim / ((float)ceil(d)+1)));
        //int k_trim = (int)(0.5*((float)n_dim / ((float)ceil(d)*(float)ceil(d))));  // TESTING TESTING TESTING        
        //mexPrintf("%d\n", k_trim);
        //int k_trim = (int)v_trim[(int)ceil(d)];
        std::nth_element(v.begin(), v.begin()+k_trim, v.end());
        D[j*n+i] = floor(d);
        for (k=0; k<k_trim; k++) {
          D[j*n+i] += (float)v[k];
        }
      }
      else {
        D[j*n+i] = d;
      }
#endif



#if 0
      // Signal
      int l,b;
      n_median = (int)((float)v_global.size()/(float)2);
      int n_support = n_lobes*n_bands;
      float d_signal;
      d_signal = 0;
      for (k=0; k<n_scales; k++) {
        for (l=0; l<n_lobes; l++) {
          for (b=0; b<n_bands; b++) {
            v_global[l*n_bands+b] = (P[i*n_dim + k*n_support + l*n_bands + b]-Q[j*n_dim + k*n_support + l*n_bands + b])*(P[i*n_dim +k*n_support + l*n_bands + b]-Q[j*n_dim +k*n_support + l*n_bands + b]);
          }        
        }
        // sum <= median 
        std::nth_element(v_global.begin(),v_global.begin()+n_median, v_global.end());
        for (l=0; l<n_median; l++) {
           d_signal += (float)v_global[l];
        }           
      }            

      // Clutter 
      if ((d_signal) < t_truncate) {
        for (k=0; k<n_dim; k++) {
          dpq = P[i*n_dim + k]-Q[j*n_dim + k];
          v[k] = dpq*dpq;
        }
        float d_clutter;
        d_clutter = 0;
        std::nth_element(v.begin(), v.begin()+n_trim, v.end());
        for (k=0; k<n_trim; k++) {
          d_clutter += (float)v[k];
        }
        D[j*n+i] = d_clutter;
      }
      else {
        D[j*n+i] = d_signal;
      }

      
#endif

#if 0
      // Signal
      float d_signal;
      d_signal=0;
      for (int s=0; s<n_scales; s++) {
        for (int b=0; b<n_lobes; b++) {
          dpq = P[i*n_dim + s*n_bands*n_lobes + b*n_bands + ((2*b)%n_bands)]-Q[j*n_dim + s*n_bands*n_lobes + b*n_bands + ((2*b)%n_bands)];
          d_signal += dpq*dpq;
        }
      }
      // Clutter 
      for (k=0; k<n_dim; k++) {
        dpq = P[i*n_dim + k]-Q[j*n_dim + k];
        v[k] = dpq*dpq;
      }
      float d_clutter;
      d_clutter = 0;
      std::nth_element(v.begin(), v.begin()+n_trim, v.end());
      for (k=0; k<n_trim; k++) {
        d_clutter += (float)v[k];
      }
      if ((d_signal+d_clutter) < t_truncate) {
        D[j*n+i] = d_clutter;
      }
      else {
        D[j*n+i] = d_signal+d_clutter;
      }
#endif

#if 0
      // L1 distance       
      d=0;
      for (k=0; k<n_dim; k++) {
        v[k] = fabs(P[i*n_dim + k]-Q[j*n_dim + k]);
        d += v[k];
      }
//       // Euclidean distance       
//       d=0;
//       for (k=0; k<n_dim; k++) {
//         dpq = P[i*n_dim + k]-Q[j*n_dim + k];
//         v[k] = dpq*dpq;
//         d += v[k];
//       }
      
      // Truncated and trimmed
      if (d <= t_truncate) {
         // sum <= median 
         std::nth_element(v.begin(), v.begin()+n_trim, v.end());
         D[j*n+i] = 0;
         for (k=0; k<n_trim; k++) {
           D[j*n+i] += (float)v[k];
         }        
      }
      else {
        D[j*n+i] = d;  
      }
#endif


    }
  }

  // Outputs  
  plhs[0] = mxD;  // Distance
} // end of mexFunction
