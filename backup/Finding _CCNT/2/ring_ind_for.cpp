#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
void found_cycle(int **c,int *chain, int *nTT, int ***TT, int MaxCycle){
    int i,k;
    if (chain[0]==MaxCycle) return;
    chain[0]++;  /* chain[0] is the length of the chain */
    for(k=1;k<=c[chain[chain[0]]][0];k++){
        chain[chain[0]+1]=c[chain[chain[0]]][k];
//         for(i=1;i<=chain[0]+1;i++)printf("%d,",chain[i]);
        if (chain[chain[0]+1]<chain[1] || chain[chain[0]+1]==chain[chain[0]-1]) {
            
        /*  ^^^^^^^^^^^^^^^^^^^^^^^^^     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            the first atom must be min              going back*/
//             printf("2\n");
            continue;
        }
        //Check if chain[end]==chain[1:end-2]
        if (chain[chain[0]+1]==chain[1]){
            // found ring
            if(chain[chain[0]]>chain[2]){// the ring is in canonical order, so record it.
                for(i=0;i<chain[0];i++){
                    TT[chain[0]-3][nTT[chain[0]-3]][i]=chain[i+1];
                }
                nTT[chain[0]-3]++;
//                 printf("3\n");
                chain[0]--;
                return;
            }
//             printf("1\n");
            chain[0]--;
            return;
        }
//         printf("0\n");
        found_cycle(c,chain,nTT,TT,MaxCycle);    
    }
    chain[0]--;
}

void findring(int **c, int MaxCycle, int n, int *nTT, int ***TT){
	int k1,k2;
    int *chain;
    chain=(int*)calloc(MaxCycle+1,sizeof(int));  // chain is the searching path.
                                                 // chain[0] is the lenghth of the path.
 
    for(k1=1;k1<=n;k1++){
        chain[1]=k1; 
        for(k2=1;k2<=c[k1][0];k2++){
            chain[0]=1;
            chain[2]=c[chain[1]][k2];
            if (chain[2]<k1) continue; // not satisfy canonical order, the first node must be min.
            found_cycle(c,chain,nTT,TT,MaxCycle);
        }/*k2 end*/
    }/*k1 end*/  
    free(chain);
}

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){

    int i,j,k;
    
    int *T;
    int *nTT;
    int ***TT;
    int **c;
    
    int buffer=65536;    

    /*  Check for proper number of arguments. */
    if (nrhs != 2)
        mexErrMsgTxt("Two inputs are required.");
    if (nlhs != 1)
        mexErrMsgTxt("One output is required.");
    
    /* Get the input a*/
//     printf("mxIsSparse: %d\n",mxIsSparse(prhs[0]));
//     printf("mxIsLogical: %d\n",mxIsLogical(prhs[0]));
    if (!mxIsSparse(prhs[0])||mxGetM(prhs[0])!=mxGetN(prhs[0])) {
        mexErrMsgTxt("Input A must be a sparse logical square matrix.");
    }

    mwSize n=mxGetM(prhs[0]);
    bool *a;
    mwIndex  *ir, *jc;
    
    ir = mxGetIr(prhs[0]);
    jc = mxGetJc(prhs[0]);
    
    c=(int**)calloc(n+1,sizeof(int*));
    c[0]=NULL;  // c[0] is discarded
                // c[i][0] is the number of atom i
                // c[i][j] is the jth neithbor of atom i
    
    for(i=1;i<=n;i++){
        c[i]=(int*)calloc(jc[i]-jc[i-1]+1,sizeof(int));
        c[i][0]=jc[i]-jc[i-1];
        for(j=jc[i-1];j<jc[i];j++){
            c[i][j-jc[i-1]+1]=ir[j]+1;
        }
    }

//     for(i=1;i<=jc[n];i++){
//         for(j=0;j<=c[i][0];j++){
//             printf("c[%d][%d]=%d\n",i,j,c[i][j]);
//         }
//     }

    if (!mxIsNumeric(prhs[1]))
        mexErrMsgTxt("Argument must be numeric.");
    if (mxGetNumberOfElements(prhs[1]) != 1 || mxIsComplex(prhs[1]))
        mexErrMsgTxt("Argument must be non-complex scalar.");

    int MaxCycle=(int)mxGetScalar(prhs[1]);
    
    nTT=(int*)calloc(MaxCycle-2,sizeof(int));     // nTT[i] is the number of cycle with length i+3.
    TT=(int***)calloc(MaxCycle-2,sizeof(int**));  // T[i][j][k] is the kth node in jth cycle with length i+3.
    for(i=0;i<MaxCycle-2;i++){
        TT[i]=(int**)calloc(buffer,sizeof(int*));
        /*for(j=0;j<buffer;j++)TT[i][j]=calloc(i+3,sizeof(int));*/
        for(j=0;j<buffer;j++)TT[i][j]=(int *)calloc(i+3,sizeof(int));
    }
    /* Call the C subroutine. */
    findring(c,MaxCycle,n,nTT,TT);
    /* Set the output pointer to the output matrix. */
    plhs[0] = mxCreateCellMatrix(MaxCycle-2,1);
    for(i=0;i<MaxCycle-2;i++)mxSetCell(plhs[0],i,mxCreateNumericMatrix(nTT[i],i+3,mxINT32_CLASS,mxREAL));
    
    /* Create a C pointer to a copy of the output matrix. */
    for(i=0;i<MaxCycle-2;i++){
        /*T=mxGetData(mxGetCell(plhs[0],i));*/
        T=(int *)mxGetData(mxGetCell(plhs[0],i));
        for(k=0;k<i+3;k++){
                for(j=0;j<nTT[i];j++){
                    T[j+k*nTT[i]]=TT[i][j][k];
                }
        }
    }
    
// //     /* free variables (may be skipped?) */
//      for(i=1;i<=jc[n];i++){
//          free(c[i]);
//      }
//     free(c);
//     free(nTT);
//      for(i=0;i<MaxCycle-2;i++){
//          for(j=0;j<buffer;j++){
//              free(TT[i][j]);
//          }
//          free(TT[i]);
//     }
//     free(TT);
}
