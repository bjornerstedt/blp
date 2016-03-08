# Logit

This subproject is for testing aspects of logit demand. 

## Subproject 1: Simulation and estimation on individual choice data 

1. Use SimMarket to create individual data

To generate data:

* Generate X with individual and common characteristics
    * Use restaurant example
    * Individual obs corresponds to market
    * How do the characteristics of individuals identify demand? 
    Observed characteristics replace nests! How does this correspond to PC logit? 
* Generate shares from beta, alpha and sigma
* Stack shares (ie probabilities) to cover the interval [0,1]
* Take a uniform draw
* Select choice in which item falls.
 