To run test cases, use the command: 

runtests

Currently this command runs: 
    merger simulations:
        RCpkTest.m
        NLpkTest.m
    estimations
        testNLauto.m
        test_Mathias.m - implemented as one big test, despite having 5 estimations
Not implemented:
    simulated data


Overview of possible test cases (incomplete)

# Simulated market 

## Without instruments

* NL 

* RC 
 
- RC with endogenous price

- various draw methods

- merger after simulation

## with artificial instruments

- Same tests as above

## with count instruments

- Same tests as above

## Merger simulations

## Other output

* Elasticities



# Painkiller

## Nested Logit

* PK with merger
NLpkTest instead of testNLpk

## Mixed logit

* Mathias 
Fails due to either new dataset, new robust estimation or lack of FE support

* PK with merger
RCpkTest tests both CES and unit demand for markets under 1 year
Tests Halton2, halton, hypercube and quadrature

## Auto

testNLauto.m test merger using test.mat

# Other

* Estimate class

* Excel read and write tables
