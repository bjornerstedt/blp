The programs in this folder implement three dimensional Hessians. 

hess_integration_test.m - Compares the numerical with the analytic Hessian
    in MixedLogitDemandMarket. 
test_delta_hessian.m - compares the numerical and analytical delta Hessian

FOCHessianNLpk.m - find the derivative of the FOC of market equilibrium, first by
finding the share Hessian (in stages). Includes various array functions
implemented as anonymous functions.

uses:
    deltaDerivatives.m - analytical BLP delta Jacobian and Hessian
    share_hessians.m - analytical share Hessians of delta and theta
    shareJacobians.m -  analytical share Jacobians of delta and theta, extracted
                        from deltaJacobian method in MixedLogitDemand

basic functions:
    arraymult.m - is the function for three dimensional array mult
    diag_array.m - create diagonal array from matrix

    hessian.m - calculates numerical Hessian of f(x)
    crosshessian.m - calculates numerical Hessian of f(x,y)
    jacobian.m - numerical jacobian, to test analytical solutions
    jacobian_array.m - numerical jacobian of matrix valued function

Extracted functions:
    nlpart.m
    finddelta.m
    share.m

test_hess.m - tests Hessian on simple quadratic functions.

Todo: Compare mmult and arraymult