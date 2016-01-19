/**************************************************************************************************************
Bjornerstedt and Verboven, Does Merger Simulation Work, JAE: Applied
- Short version for Matlab replication
***************************************************************************************************************/

version 11
clear all


use painkillers

****************************************************************************
* PART 1. DATA MANAGEMENT
****************************************************************************


	replace Ptablets=Ptablets/(cpi/100)
	replace PX=PX/(cpi/100)
	egen totQ = sum(X) , by(time)
	sum totQ , meanonly 		
	gen BL0 = 2*r(mean)

	
// CREATE INSTRUMENTS
	egen num = count(product), by(time)
	egen numg = count(product), by(time form)
	egen numf = count(product), by(time firm)
	egen numfg = count(product), by(time firm form)
	egen numhg = count(product), by(time substance form)
	egen numfgh = count(product), by(time firm substance form)
	gen con=1
	
	global exogvar marketing1 sw sm time month2-month12 // used for nested logit
	mergersim init, nests(form substance) unit price(Ptablets) quantity(Xtablets) marketsize(BL0) firm(firm)
	xtivreg2 M_ls $exogvar (Ptablets M_lsjh M_lshg = num*) if year<2009, fe // robust


mergersim market if year==2008&month==12   

mergersim simulate if year == 2008&month==12, seller(1) buyer(2) sellereff(0) buyereff(0)

**** Save for Matlab
preserve
keep PX Ptablets BL0 M_costs M_price2 firm
export delimited using "test.csv", quote replace
restore
**** Save for Matlab

