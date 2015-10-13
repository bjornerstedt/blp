// findDeltaCpp calculates delta using the contraction for a single market

#include "armaMex.hpp"
void checkInput(int nrhs, const mxArray *prhs[]);

mat findDeltaCpp(mat s, mat expmu, mat iweight, mat edelta)
{
    int i = 0;
    double maxdev = 100;
    double tolerance = 1e-14;
    const int fpmaxit = 1000;
    mat newedel;
    mat id = ones<mat>(1, expmu.n_cols);
    while(maxdev > tolerance && i < fpmaxit) {
        mat eg = expmu.each_col() % edelta; 
        newedel = edelta % s / ( eg.each_row() % (id/(id + sum(eg))) * iweight);
        maxdev = as_scalar(max(abs(newedel - edelta)));
        edelta = newedel;
        i = i + 1;
    }
    return(newedel);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    checkInput(nrhs, prhs);
    mat s = armaGetPr(prhs[0]);
    mat expmu = armaGetPr(prhs[1]);
    mat iweight = armaGetPr(prhs[2]);
    mat edelta = armaGetPr(prhs[3]);
        
    mat C = findDeltaCpp(s, expmu, iweight, edelta);
    plhs[0] = armaCreateMxMatrix(C.n_rows, C.n_cols);
    armaSetPr(plhs[0], C);
    return;
}

void checkInput(int nrhs, const mxArray *prhs[])
{
    const int args = 4;
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
