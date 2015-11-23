The Development/market folder contains modifications to the demand and market classes to 
estimate and simulate market by market, rather than having all markets in the class.

MixedLogitDemandNew is new subclass of MixedLogitDemand
MixedLogitDemand2 is the old version
MixedLogitDemandMarket is the market demand class
std_tests tests.m that MixedLogitDemandNew does the same thing as MixedLogitDemand
test.m compares MixedLogitDemandNew and MixedLogitDemand, bottom up.
