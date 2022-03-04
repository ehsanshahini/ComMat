/*
 * =============================================================
 * conn_c.c - calculate connectivity matrix for a dual lattice
 *
 * This is a MEX-file for MATLAB.
 * Copyright (c) 2008/11/28 2:43 byj lab, NTUCH
 * =============================================================
 */

/* $Revision: 1.1.4.15 $ */

#include "mex.h"
#include <math.h>

void connec(double *R, bool *A, int m, int ncols)
{
  int i,j,k,l,count = 0;
  double eps=0.00000001;
  
  if (ncols==2){
      for (k = 0; k < m; k++) {
          for (l = 0; l < k; l++){
              count = 0;
              for (i = 0; i < 3; i++) {
                  for (j = 0; j < 3; j++) {
                      if( fabs(*(R+i+3*k)-*(R+j+3*l))<eps && fabs(*(R+i+3*k+m*3)-*(R+j+3*l+m*3))<eps    ){
                          count++;
                      };
                  }
              }
              if (count==2){
                  *(A+k+m*l) = true;
                  *(A+l+m*k) = true;
              }
          }
      }
  }
  else{
      for (k = 0; k < m; k++) {
          for (l = 0; l < k; l++){
              count = 0;
              for (i = 0; i < 3; i++) {
                  for (j = 0; j < 3; j++) {
                      if(fabs(*(R+i+3*k)-*(R+j+3*l))<eps && fabs(*(R+i+3*k+m*3)-*(R+j+3*l+m*3))<eps && fabs(*(R+i+3*k+2*m*3)-*(R+j+3*l+2*m*3))<eps ){
                           count++;
                      };
                  }
              }
              if (count==2){
                  *(A+k+m*l) = true;
                  *(A+l+m*k) = true;
              }
          }
      }      
  }
}

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  double *R;
  double *A;
  int status,mrows,ncols;
  
  /*  Check for proper number of arguments. */
  /* NOTE: You do not need an else statement when using 
     mexErrMsgTxt within an if statement. It will never 
     get to the else statement if mexErrMsgTxt is executed. 
     (mexErrMsgTxt breaks you out of the MEX-file.) 
  */ 
  if (nrhs != 1) 
    mexErrMsgTxt("One input required.");
  if (nlhs != 1) 
    mexErrMsgTxt("One output required.");
  
 
  /* Check to make sure the first input argument is real. */
  if (mxIsComplex(prhs[0])) {
    mexErrMsgTxt("Input R must be a real matrix.");
  }
  
  /* Get the the input R. */
  R = mxGetPr(prhs[0]);
  
 
  /* Get the dimensions of the matrix input R. */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  
  if ( ncols != 3 && ncols !=2) {
    mexErrMsgTxt("The number of columns of input R must be 2 or 3.");
  }

  if ( mrows/3 != (double) mrows/3 ) {
    mexErrMsgTxt("The number of rows of input R must be multiple of 3.");
  }

  
  /* Set the output pointer to the output matrix. */
  plhs[0] = mxCreateLogicalMatrix(mrows/3, mrows/3);
  
  
  /* Create a C pointer to a copy of the output matrix. */
  A = mxGetPr(plhs[0]);
  
  /* Call the C subroutine. */
  connec(R,A,mrows/3,ncols);
}