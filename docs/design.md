---
title: "SimMarket Design"
author: "Jonas Björnerstedt"
date: ""
output: pdf_document
---

# Introduction

The purpose of this document is to collect documentation of various design
issues in SimMarket. It does not describe how to use the program or basic
programming issues that are best described by the code itself. The purpose is to
provide a background for more challenging issues that can be hard to understand
from reading the program code.

# Data structure and initialization

The demand code is used in three basic ways

1. Estimation: $p$ and $q$ known, parameters $\alpha$, $\sigma$, $d$, and $c$ unknown.

    (a) Estimation is done on all markets but in RC, finddelta calculates on
    each market at a time.
    
    (a) findCosts calculates c. Note that this has to be always tested in order
    for the test of the mapping to be complete

2. Simulation: $\alpha$, $\sigma$, $d$ and $c$ known, $p$ and $q$ unknown

    (a) In RC done on all markets at one time - share() splits per market
    
    (a) In NL done market by market

3. Generation: $\alpha$, $\sigma$, $d$ and $p$ known, $q$ and $c$ unknown. 

    (a) From $p$ and $d$ delta is generated, then shares

    (a) With market simulations, $p$ is unknown but $c$ is known. 

Initialization of the Estimate and derived classes (Market and 
Demand classes), 

1. Initialization is performed with the init() method. 

    (a) init(selection) is invoked by equilibrium(selection), 
    findCosts(selection), estimate(selection), estimateGMM()
    
    (a) init() is invoked by Market.init(). When is this necessary?:

        i. Estimation can be on selection. But simulation outside selection is
        problematic, it requires that $d$ is defined for new products.
        
        i. An alternative could be to test whether Demand.init() has been run.

2. initSimulation(market number) is used prior to each market simulation.

    (a) initSimulation is invoked by Market.initSimulation(), 

    (a) which in turn is invoked by findCosts() and equilibrium()
    
## Scenarios

Note that init is used both in initializing estimation and simulation. Three
scenarios have slightly different setup requirements.

1. Estimate on data

2. Simulate on existing data

3. Create new data

# Class parameters

## Estimate class

1. Estimation variables, optionally restricted to a subset of the data

    (a) X, Z and y matrices
    
    (a) Xorig and Zorig for non-demeaned calculations when fixed effects
    estimation is used.

2. Organization variables, defined for all data. Panelid is defined for all to
allow estimation and simulation to be performed on different samples of the
data.

    (a) obj.panelid - panel identifier (all data)

    (a) obj.vars - display vars

    (a) obj.lsdv - LSDV, persistent local var of init. Obsolete? (selection)

        i. Created if it does not exist or if there is no selection
        
        ii. Drops zero dummies (on selection if set), ie simulation cannot be
        done outside selection

## NLDemand class

1. obj.marketid and obj.dummarket - in Demand and Market classes this identifies
markets - selection

2. obj.nestcount - nest variable names. Obsolete?

    (a) Used in estimate(), shares(), shareJacobian(), elasticities(), generateShares() -
    to get number of nests

3. obj.nest - nesting vars (no selection)

    (a) Set in initAdditional

    (a) Used in initSimulation(), selects

4. obj.p - non-demeaned price (no log) (selection). 

    (a) used in shareJacobian on selection

    (b) used in generateShares(), demand calc.

    (c) used in estimate as log variable!

5. obj.ms - marketsize. Used in share calcs and demand calc. (selection)

6. obj.share - created shares. 

    (a) Created in initAdditional if depvar does not exist! (selection)

    (b) Used by NLDemand.estimate() to create obj.d

    (c) actualDemand returns obj.share.s

    (d) shareJacobian uses selection to calculate if isempty(P)

7. obj.G, H, GG, HH, GH

    (a) Created by initSimulation

    (b) Used by shares, shareJacobian and elasticities

8. obj.d

    (a) created by estimate()

    (b) initSimulation() sets obj.sim.d from obj.d on selection

    (c) obj.sim.d used in shares()

9. obj.alpha and obj.sigma obj.beta!

10. obj.sim - structure with selection, $d$ and market. Only market is really
needed. Remove?

11. obj.data

    (a) used in initAdditional to create obj.nest, obj.marketid, obj.shares, obj.q and obj.p
    
    (a) used in groupElasticities() for grouping

## RCDemand class

1. obj.draws - a Draws class instance with nonlinear draws. 

    (a) Only created once - why? 

        i. Market invokes RCDemand.init(), thus creating different draws.
        
        ii. Simulation should have the same draws as estimation.
        
        iii. Note that estimation will have the same draws as simulated market,
        if the same instance is used.

    (b) What happens in copy? What do we want to happen?! 

        i. Copying is done indirectly in a merger. There we want the same draws.
        
        ii. Copying demand and changing draw settings will have no effect,
        however.

2. obj.period - list of RCDemandMarket class instances (selection)

    (a) Created and initialized by initPeriods() called from init()
    
    (b) initSimulation calculates sets obj.d in each period. When is this used?
    
    (c) Used by methods to call corresponding method: findDelta(), shares(),
    shareJacobian(), deltaJacobian()
    
    (d) sharesAll() - where is this used? Not in SimMarket.

3. obj.x2 - nonlinear data matrix, with logprice in CES. 

    (a) Set in init()
  
    (b) Not demeaned in FE. (selection)

4. obj.nonlinparams - list of all nonlinear. 

    (a) Created in obj.randdraws() from init()
  
    (b) Used by methods: RCDemandMarket, 
  
    (c) Used by methods to get length : getSigma()

5. obj.vars2 - nonlin variable names. Replace with function.

6. obj.d

    (a) set in initSimulation() from edelta and alpha after estimation. Move to
    estimate?
    
    (b) Set manually in SimMarket for simulation
  
    (c) used by RCDemandMarket

7. obj.p

    (a) Set in init@NLDemand()
  
    (b) Used by RCDemandMarket

8. obj.Xorig(:, 1)

    (a) Used by initEstimation() - to create pHat and deltaHat (Xhat(:, 2:end))
  
    (b) initSimulation() - to calculate obj.d

9. obj.data

    (a) used in init to create obj.x2

## RCDemandMarket class

1. obj.alpha and obj.sigma

    (a) Set in constructor
  
    (b) Used in shares and shareJacobian()

2. obj.x2 - for estimation.

    (a) creates obj.vx in init()

3. obj.v and obj.iweight

4. obj.s - why?

    (a) findDelta(): to find delta giving shares

5. obj.p - why? :

    (a) in shareJacobian(): To calc on obj.p. Needs non-log price, x2 has log
    price in CES.

6. obj.d - why?:

    (a) in shares(): to do share calc

7. obj.nonlinprice - from RCDemand nonlinparams - strcmp(demand.var.price,
demand.nonlinparams)

    (a) shareJacobian() to update price col in data in simulation.

Check that markets are updated with all changes. 

* When is init invoked?

## Market class

1. obj.data - copy of demand.data

    (a) Set in estimate() and estimateGMM() and used directly in init()
  
1. obj.demand 

    (a) Set in constructor or directly
  
    (b) Used in init(), 
    
    (c) initSimulation() invokes obj.demand.initSimulation()
    
    (d) equilibrium() and fixedPoint() invoke obj.demand.shares(mkt)

    (d) foc() and margins() invoke obj.demand.shares(mkt) and obj.demand.shareJacobian(mkt)

    (d) findCosts() invokes obj.demand.shareJacobian(mkt)

    (d) summary() and compare() use obj.demand.data and obj.demand.var.market

1. obj.firm 

    (a) Set in init(), if not categorical transformed with unique(obj.firm)
    
        * Only set if not set before.

    (b) Can be directly set by user

    (b) Used in initSimulation(), summary() and compare()
    
1. obj.p

    (a) Set in init()
  
    (a) created in equilibrium() and fixedPoint()
  
    (b) Used in findCosts(), with obj.s

    (b) Used in getMarketShares(), with obj.q

    (b) Used in summary() and compare()

1. obj.c 

    (a) Created in findCosts()
  
    (b) Used in init(), equilibrium(), residuals(), estimateCosts(), fixedPoint()

    (b) Used in summary() and compare()

1. obj.RR 

    (a) Set in initSimulation()
  
    (b) Used in findCosts(), foc(), margins() 

1. obj.gamma 

    (a) Set by user and SimMarket
  
    (b) Used in foc()

1. obj.marketid 

    (a) Set in init()
  
    (b) Used in initSimulation(), equilibrium(), findCosts(), getMarketShares(), summary(), compare(),

## SimMarket class

# Requirements

### Estimate

For estimation, Estimate.var.panel has to be set unless
Estimate.settings.paneltype == 'none'.

### Demand classes

For NLDemand in init: 

1. obj.var.market and obj.var.price must be set

    (a) In RCDemand obj.var.nonlinear must be set

2. either 

    * obj.var.quantity has to be in the data - estimation 

        – obj.var.marketsize and obj.var.exog must be set

    * obj.alpha has to be set - simulation

To calculate shares() NLDemand needs obj.d, obj.alpha, obj.sigma and
corresponding obj.var.nests. The method initSimulation requires obj.nestlist to
be set

### Market class

The Market class has direct and indirect requirements. In init

1. obj.var.firm and obj.demand have to be set

2. Estimation requires obj.var.depvar to be set. Not the best 

3. From obj.demand, it copies

    (a) obj.demand.p
  
    (b) obj.demand.q
  
    (c) obj.demand.marketid
  
    (d) obj.demand.panelid - not currently used, but can be used for panel
    estimates of costs

From the obj.demand object, it requires

1. obj.demand.share.s (alt: obj.demand.actualDemand() )

2. obj.demand.shares(p) where market number has been set in 
obj.demand.initSimulation(marketnumber)

3. obj.demand.shareJacobian([]) and obj.s(sel) in findCosts. Invoke
obj.demand.initSimulation(market_number) first to set market number. This two
step process is used to avoid repeated initialization of market

4. obj.demand.shareJacobian( p ) and obj.demand.shares(pt)

