The painkiller folder contains scripts to estimate and simulate with 
painkiller data. The important scripts are

pk_bootstrap.m
RC: Estimates either multiple starting points or selects one to do bootstrapped
merger simulation. 
Invokes:
    pkRC.m - RC invokes to estimate multiple/selected starting points
    pkNL.m - NL invokes to estimate 
    pkRC_output.m - to save RC multiple starting point output to Excel
    pk_bootstrap_output - to save bootstrap simulations to Excel

pk_cost_calcs.m
    To create product level costs pre- and post-merger for 
    with and without collusion, ces/unit demand, and RC/NL