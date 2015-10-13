// shareCalcCpp calculates shares

#include "armaMex.hpp"
void checkInput(int nrhs, const mxArray *prhs[]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    checkInput(nrhs, prhs);
    mat expmu = armaGetPr(prhs[0]);
    mat iweight = armaGetPr(prhs[1]);
    mat edelta = armaGetPr(prhs[2]);
        
    mat eg = expmu.each_col() % edelta; 
    mat id = ones<mat>(1, expmu.n_cols);
    mat si = eg.each_row() % (id / (id + sum(eg)) );
    mat s = si * iweight;

    if(nlhs > 0) {
        plhs[0] = armaCreateMxMatrix(s.n_rows, s.n_cols);
        armaSetPr(plhs[0], s);
    }
    if(nlhs == 2) {
        plhs[1] = armaCreateMxMatrix(si.n_rows, si.n_cols);
        armaSetPr(plhs[1], si);
    }
    
    return;
}

void checkInput(int nrhs, const mxArray *prhs[])
{
    const int args = 3;
    // Check the number of input arguments.
    if (nrhs != args)
        mexErrMsgTxt("Incorrect number of input arguments.");
    for(int i=0; i<args; i++) {
        // Check type of input.
        if ( (mxGetClassID(prhs[i]) != mxDOUBLE_CLASS) )
            mexErrMsgTxt("Input must be of type double.");
        // Check if input is real.
        if ( (mxIsComplex(prhs[i]))  )
            mexErrMsgTxt("Input must be real.");
    }
}
