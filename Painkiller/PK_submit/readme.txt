Painkiller submit version, commit 121 of BLP code.

The Program folder contains the general estimation and simulation code. The Painkiller directory the code used in the analysis. 

The pk_bootstrap.m script is used for estimation of 
1) nested and random coefficient logit. 
2) CES and unit demand

The choice of estimation is made with boolean variables RC and ces at the top of the script. Estimation can consist of either

a) many different starting points and individual draws
b) bootstrap calculation of costs and merger effects from a selected estimation.

In both cases the calculations are performed by invoking either pkRC.m or pkNL.m that contain estimation settings for mixed and nested logit, respectively. Data input and estimation output is in the form of Matlab table objects. Invocation of pkRC.m is through a utility function batchrun, to facilitate parallel processing of multiple starting points.

Estimation and merger simulations are calculated using Matlab classes in the Program folder. Excel output of the results is in pk_bootstrap_ouptut.m.

The script pk_cost_calcs.m performs post merger cost calculations for mixed and nested logit. The script also invokes the pkRC.m and pkNL.m files mentioned above, using the estimation results in cost calculations.