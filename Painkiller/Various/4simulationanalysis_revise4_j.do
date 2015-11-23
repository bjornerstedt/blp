/*******************************************************************************************
This program performs merger simulations for different scenarios and descriptive table.
This is for revisions of JAE: Applied
*******************************************************************************************/

version 11
clear all
set more on
// Both product and time variables currently obtained through xtset


// CHOOSE MAIN SETTINGS
local ces 1	// 1 is ces 0 is unit
local cond 0 	// between 0 and 1 for conduct

****************************************************************************
* PART 0. READING DATA AND CREATING OF OUTSIDE GOOD AND PRICE/QUANTITIES
****************************************************************************

* READ OLD DATASET 9511

use painkillers9511main1
drop lP form2

replace GDPnom=770000 if year>2008 // no GDP info yet after 2008

// CREATE OUTSIDE GOOD
if `ces'==1 {
	egen totQ = sum(PX1), by(time)
	sum totQ if year<2009, meanonly
		local medincome = 649144.4 // time varying potential budget measure (very similar to budget0) (649144 is income in median period time=84)
	gen BL`ces' = 2*r(mean)*(GDPnom/`medincome'	)
	drop totQ
}
if `ces'==0 {
	replace Ptablets=Ptablets/(cpi/100)
	replace PX=PX/(cpi/100)
	*replace PX1=PX1/(cpi/100)
	egen totQ = sum(X) , by(time)
	sum totQ /*if year<2009*/, meanonly 		// sometimes be sure to drop post2009 data
	gen BL`ces' = 2*r(mean)
	*drop totQ
}


****************************************************************************
* PART 1. DATA MANAGEMENT
****************************************************************************

replace firm=15 if brand==9&year>2008 // Bamyl is brand=9: Ellem sold it to Meda on 1 April 2008 and this was misrecorded for first quarter after 2008


// PRODUCT DUMMIES AND TIME-INVARIANT CHARACTERISTICS

	replace product=200 if product==71 // we loose value label, but this product is otherwise prod44 and is missing for the premerger period, so it becomes more difficult to construct the delta and mcost
	quietly{
	tab product,gen(prod)
	}
	tab form,gen(form)
	tab substance,gen(substance)
	gen lpacksize=log(packsize)
	gen ldosage=log(dosage)


// CREATE INSTRUMENTS

	egen num = count(product), by(time)
	egen numg = count(product), by(time form)
	egen numf = count(product), by(time firm)
		egen numfg = count(product), by(time firm form)
	egen numhg = count(product), by(time substance form)
		egen numfgh = count(product), by(time firm substance form)

	gen con=1
	global charact con ldosage lpacksize form2 substance2 substance3

	foreach var of varlist $charact {
		bysort time: egen sum1_`var' = sum(`var')
		bysort time firm: egen sum2_`var' = sum(`var')
		gen i1_`var' = sum2_`var' - `var'
		gen i2_`var' = sum1_`var' - sum2_`var'
		}

global exogvar marketing1 sw sm time month2-month12 // used for nested logit
global exogvar1 marketing1 sw sm month2-month12 // used for logit and BLP

gen time_p=time*(substance==1)
gen time_i=time*(substance==3)
gen time_a=time*(substance==5)

global exogvar2 marketing1 sw sm month2-month12 time_p time_i time_a

global iv i1_ldosage i1_lpacksize i2_ldosage i2_lpacksize 



****************************************************************************
* PART 1.1 LOGIT DEMAND ANALYSIS WITH BLP INSTRUMENTS (AS INTRO TO BLP MODEL)
****************************************************************************

*logit with BLP instruments i1_* and i2_*

if `ces'==1 {
mergersim init, ces price(Ptablets) revenue(PX1) marketsize(BL`ces') firm(firm)
	xtivreg2 M_ls $exogvar1 (M_lp = i1_* i2_* ),fe robust
}
if `ces'==0 {
mergersim init, unit price(Ptablets) quantity(Xtablets) marketsize(BL`ces') firm(firm)
	xtivreg2 M_ls $exogvar1 (Ptablets = i1_* i2_* ),fe robust
}
mergersim market if year<2009
mergersim simulate if year == 2008 & month==12, seller(1) buyer(2) sellereff(0) buyereff(0) conduct(0) 


/*
****************************************************************************
* PART 1.2 	BASIC EXPLORATORY DEMAND ANALYSIS AND MERGER SIMULATION
****************************************************************************

if `ces'==1 {
mergersim init, nests(form substance) ces price(Ptablets) revenue(PX1) marketsize(BL`ces') firm(firm)
	xtivreg2 M_ls $exogvar (M_lp M_lsjh M_lshg = num*),fe robust 
	*xtivreg2 M_ls $exogvar1 (M_lp M_lsjh M_lshg = num* ), fe robust  // gives lower elasticity
	*xtivreg2 M_ls $exogvar2 (M_lp M_lsjh M_lshg = num* ), fe robust 

	*xtivreg2 M_ls $exogvar (M_lp M_lsjh M_lshg = num* i1_* i2_* ), fe robust  // gives higher elasticity
	*xtivreg2 M_ls $exogvar1 (M_lp M_lsjh M_lshg = num* i1_* i2_* ), fe robust  // gives very low elasticity
	*xtivreg2 M_ls $exogvar2 (M_lp M_lsjh M_lshg = num* i1_* i2_* ), fe robust
	
	***xtivreg2 M_ls $exogvar (M_lp M_lsjh M_lshg = num* $iv ), fe robust  // may use this as alternative!
	*xtivreg2 M_ls $exogvar1 (M_lp M_lsjh M_lshg = num* $iv ), fe robust  
	*xtivreg2 M_ls $exogvar2 (M_lp M_lsjh M_lshg = num* $iv ), fe robust
}

if `ces'==0 {
mergersim init, nests(form substance) unit price(Ptablets) quantity(Xtablets) marketsize(BL`ces') firm(firm)
	xtivreg2 M_ls $exogvar (Ptablets M_lsjh M_lshg = num*),fe robust
	*xtivreg2 M_ls $exogvar1 (Ptablets M_lsjh M_lshg = num*),fe robust
	*xtivreg2 M_ls $exogvar2 (Ptablets M_lsjh M_lshg = num*),fe robust

	*xtivreg2 M_ls $exogvar (Ptablets M_lsjh M_lshg = num* i1_* i2_*),fe robust
	*xtivreg2 M_ls $exogvar1 (Ptablets M_lsjh M_lshg = num* i1_* i2_*),fe robust
	*xtivreg2 M_ls $exogvar2 (Ptablets M_lsjh M_lshg = num* i1_* i2_*),fe robust

	***xtivreg2 M_ls $exogvar (Ptablets M_lsjh M_lshg = num* $iv),fe robust // may use this as alternative!
	*xtivreg2 M_ls $exogvar1 (Ptablets M_lsjh M_lshg = num* $iv),fe robust
	*xtivreg2 M_ls $exogvar2 (Ptablets M_lsjh M_lshg = num* $iv),fe robust
}

mergersim market if year==2008 & month==12, conduct(`cond')   // to know elasticities in december 2008

mergersim simulate if year == 2008 & month==12, seller(1) buyer(2) sellereff(0) buyereff(0) conduct(`cond') detail
gen M_lerner=(Ptablets-M_costs)/Ptablets
table substance if year==2008 & month==12,contents(count M_price2 mean M_lerner mean M_price2 mean M_price_ch sum PX1) 
table firm substance if year==2008 & month==12,contents(mean M_price_ch) 


/* CONCLUSION (on 8/12/2014)
Results are sensitive to instruments and trend.
Ideally do not include trend to avoid the issue, but then elasticity can be very low. 
Results somewhat less sensitive if interact trend with substance, but still not satisfying.

Adding i1_* and i2_* can give very high elasticities or very low elasticities
Adding only iv (i.e. limited BLP instruments, non-form and non-substance) seems to give comparable results to specification from case

Seems best to stay with the simplest specification from the case, or add limited instruments from BLP (i.e.iv)
*/


*/

****************************************************************************
* PART 2. 	SECOND STAGE ANALYSIS FOR DEMAND AND COST
****************************************************************************

// 2.1 DEMAND ANALYSIS (PICK SPECIFICATION FROM PART 1.2)

* first stage (use product dummies instead of fixed effects)

if `ces'==1 {
mergersim init, nests(form substance) ces price(Ptablets) revenue(PX1) marketsize(BL`ces') firm(firm)
	ivreg2 M_ls $exogvar prod1-prod56 (M_lp M_lsjh M_lshg = num*), noconstant robust
}

if `ces'==0 {
mergersim init, nests(form substance) unit price(Ptablets) quantity(Xtablets) marketsize(BL`ces') firm(firm)
	ivreg2 M_ls $exogvar prod1-prod56 (Ptablets M_lsjh M_lshg = num*), noconstant robust
}
mergersim market if year<2009, conduct(`cond')  // do this to compute the marginal costs for later


* second stage: Nevo fixed effects on time-invariant characteristics

**keep if e(sample)
matrix beta=e(b)
matrix delta=beta[1,"prod1".."prod56"]'

preserve
clear
svmat delta
save temp, replace
restore

preserve
sort product
gsort product -year -month
by product: keep if _n==1
merge 1:1 _n using temp
reg delta ldosage lpacksize form2 substance2 substance3
 
restore

* new delta after the merger

gen delta_new=M_delta+(log(packsize_new)-lpacksize)*_b[lpacksize]

// 2.2 MARGINAL COSTS AND CHANGE BECAUSE OF PACKAGING

gen lcost=log(M_costs)

* first stage

reg lcost time month2-month12 prod1-prod56 if year<2009, noconstant robust

* second stage: Nevo fixed effects on time-invariant characteristics

**keep if e(sample)
matrix beta=e(b)
matrix lcost_fe=beta[1,"prod1".."prod56"]'

preserve
clear
svmat lcost_fe
save temp, replace
restore

preserve
sort product
gsort product -year -month
by product: keep if _n==1
merge 1:1 _n using temp
reg lcost_fe ldosage lpacksize form2 substance2 substance3 

restore
scalar ps= _b[lpacksize]

* new marginal cost after the merger

gen costs_new=M_costs*(packsize_new/packsize)^_b[lpacksize]
gen cost_ch=(costs_new-M_costs)/M_costs



****************************************************************************
* PART 3. MERGER SIMULATION - BASE AND EXTENDED SCENARIOS
****************************************************************************

gen firmsubst=firm
replace firmsubst=15 if firm==5 // assign Ellem to Meda
replace firmsubst=15 if firm==6 // assign Recip to Meda
replace firmsubst=99 if brand==15 // make Ipren of McNeil separate "firm", i.e. McNeil-Ibu
replace firmsubst=199 if brand==5 // make Alindrin of Meda separate "firm", i.e. Meda-Ibu
label define firmsubstlabel 1 "AZT" 2 "GSK" 3 "McNeil" 4 "Nycomed" 13 "Bayer" 14 "McNeil-ASA"  15 "Meda-ASA" 16 "Searle" 99 "McNeil-Ibu" 199 "Meda-Ibu"
label values firmsubst firmsubstlabel
tab firm if year==2008
tab firmsubst if year==2008


* must repeat the above demand analysis first to get the parameters again
if `ces'==1 {
mergersim init, nests(form substance) ces price(Ptablets) revenue(PX1) marketsize(BL`ces') firm(firm)
	xtivreg2 M_ls $exogvar (M_lp M_lsjh M_lshg = num*),fe robust /*first*/
}
if `ces'==0 {
mergersim init, nests(form substance) unit price(Ptablets) quantity(Xtablets) marketsize(BL`ces') firm(firm)
	xtivreg2 M_ls $exogvar (Ptablets M_lsjh M_lshg = num*),fe robust /*first*/
}

mergersim market if year==2008&month==12, conduct(`cond')   // to know elasticities


* 3.1 MERGER SIMULATION WITHOUT COST CHANGES

mergersim simulate if year == 2008&month==12, seller(1) buyer(2) sellereff(0) buyereff(0) conduct(`cond') detail
gen M_lerner=(Ptablets-M_costs)/Ptablets
egen M_tot=sum(M_quantity), by(year month)
egen M_tot2=sum(M_quantity2), by(year month)
gen M_sh=M_quantity/M_tot
gen M_sh2=M_quantity2/M_tot2
gen M_shdif=M_sh2-M_sh

*replication of above to check if correct
table firm if year==2008&month==12, contents(mean M_price_ch) 
table firm if year==2008&month==12, contents(sum M_sh sum M_sh2)
*by substance
table substance if year==2008&month==12, contents(mean M_price_ch)
table substance if year==2008&month==12 [fw=Xtablets], contents(mean M_price_ch)
table substance if year==2008&month==12, contents(sum M_sh sum M_sh2 sum M_shdif)

*by firm*substance
table firmsubst if year==2008&month==12, contents(mean M_price_ch)
table firmsubst if year==2008&month==12 [fw=Xtablets], contents(mean M_price_ch)
table firmsubst if year==2008&month==12, contents(sum M_sh sum M_sh2)


drop M_lerner M_tot M_tot2 M_sh M_sh2 M_shdif


* 3.2 MERGER SIMULATION WITH COST CHANGES BECAUSE OF PACKAGE SIZE
* packagesize reduction for Paracetamol + product removals of high packages sizes of Treo and Banyl, who maintain twin products with low package size

mergersim simulate if year == 2008&month==12, seller(1) buyer(2) newcosts(costs_new) conduct(`cond') //method(fixedpoint)
gen M_lerner=(Ptablets-M_costs)/Ptablets
egen M_tot=sum(M_quantity), by(year month)
egen M_tot2=sum(M_quantity2), by(year month)
gen M_sh=M_quantity/M_tot
gen M_sh2=M_quantity2/M_tot2
gen M_shdif=M_sh2-M_sh

*by substance
table substance if year==2008&month==12, contents(mean M_price_ch)
table substance if year==2008&month==12 [fw=Xtablets], contents(mean M_price_ch)
table substance if year==2008&month==12, contents(sum M_sh sum M_sh2 sum M_shdif)
*by firm*substance
table firmsubst if year==2008&month==12, contents(mean M_price_ch)
table firmsubst if year==2008&month==12 [fw=Xtablets], contents(mean M_price_ch)
table firmsubst if year==2008&month==12, contents(sum M_sh sum M_sh2 sum M_shdif)

drop M_lerner M_tot M_tot2 M_sh M_sh2 M_shdif

* cost change from packsize reduction
table substance if year==2008&month==12,contents(count cost_ch mean cost_ch sd cost_ch) 
table substance if year==2008& month==12 [fw=Xtablets],contents(count cost_ch mean cost_ch sd cost_ch) 
table firmsubst if year==2008& month==12,contents(count cost_ch mean cost_ch sd cost_ch) 
table firmsubst if year==2008& month==12 [fw=Xtablets],contents(count cost_ch mean cost_ch sd cost_ch) 


* 3.3 AS IN 3.2 PLUS COST CHANGE TO RATIONALIZE AVERAGE PRICE INCREASE

* (ONLY FOR CES=1)

*drop M_lerner M_tot M_tot2 M_sh M_sh2 costs_new1

gen costs_new1=costs_new
if `cond'==0 {
*replace costs_new1=costs_new*0.9944 if firmsubst==1
*replace costs_new1=costs_new*0.7724 if firmsubst==2
replace costs_new1=costs_new*0.9912 if firmsubst==1
replace costs_new1=costs_new*0.7702 if firmsubst==2
replace costs_new1=costs_new*1.0052 if firmsubst==4
replace costs_new1=costs_new*1.1082 if firmsubst==13
replace costs_new1=costs_new*1.1237 if firmsubst==14
replace costs_new1=costs_new*0.9880 if firmsubst==15
replace costs_new1=costs_new*1.0069 if firmsubst==99
replace costs_new1=costs_new*0.9990 if firmsubst==199
}


if `cond'==0.75 {
*replace costs_new1=costs_new*1.0246 if firmsubst==1
*replace costs_new1=costs_new*0.8712 if firmsubst==2
replace costs_new1=costs_new*1.0291 if firmsubst==1
replace costs_new1=costs_new*0.8750 if firmsubst==2
replace costs_new1=costs_new*0.9760 if firmsubst==4
replace costs_new1=costs_new*1.0810 if firmsubst==13
replace costs_new1=costs_new*1.1025 if firmsubst==14
replace costs_new1=costs_new*0.9770 if firmsubst==15
replace costs_new1=costs_new*0.9780 if firmsubst==99
replace costs_new1=costs_new*0.9730 if firmsubst==199
}

mergersim simulate if year==2008&month==12, seller(1) buyer(2) newcosts(costs_new1) conduct(`cond') //method(fixedpoint)
gen M_lerner=(Ptablets-M_costs)/Ptablets
egen M_tot=sum(M_quantity), by(year month)
egen M_tot2=sum(M_quantity2), by(year month)
gen M_sh=M_quantity/M_tot
gen M_sh2=M_quantity2/M_tot2
gen M_shdif=M_sh2-M_sh

*by substance
table substance if year==2008&month==12, contents(mean M_price_ch)
table substance if year==2008&month==12 [fw=Xtablets], contents(mean M_price_ch)
table substance if year==2008&month==12, contents(sum M_sh sum M_sh2 sum M_shdif)
*by firm*substance
table firmsubst if year==2008&month==12, contents(mean M_price_ch)
table firmsubst if year==2008&month==12 [fw=Xtablets], contents(mean M_price_ch)
table firmsubst if year==2008&month==12, contents(sum M_sh sum M_sh2 sum M_shdif)


drop M_lerner M_tot M_tot2 M_sh M_sh2 M_shdif


mmmmmmmmmmmmmmmmmmmmmmmmmmmm


* CONSIDER PRICE INCREASE DIRECTLY

generate newp = Ptablets
replace newp = 1.428*Ptablets if firmsubst == 1
replace newp = 1.461*Ptablets if firmsubst == 2
replace newp = 1.012*Ptablets if firmsubst == 4
replace newp = 1.108*Ptablets if firmsubst == 13
replace newp = 1.153*Ptablets if firmsubst == 14
replace newp = 1.070*Ptablets if firmsubst == 15
replace newp = 1.024*Ptablets if firmsubst == 99
replace newp = 1.001*Ptablets if firmsubst == 199

/*
generate newp = Ptablets
replace newp = 1.436*Ptablets if substance == 1
replace newp = 1.020*Ptablets if substance == 3
replace newp = 1.141*Ptablets if substance == 5
*/

*mergersim init, nests(form substance) ces price(Ptablets) revenue(PX1) marketsize(BL`ces') firm(firm) alpha(-0.3037296) sigmas(0.8350244 0.6670308)
mergersim init, nests(form substance) ces price(Ptablets) revenue(PX1) marketsize(BL`ces') firm(firm) alpha(-0.3037296) sigmas(0.8350244 0.6670308)
mergersim market if year==2008&month==12
mergersim demand if year == 2008 & month == 12, price(newp)


gen M_price_ch=(newp-Ptablets)/Ptablets

egen M_tot=sum(M_quantity), by(year month)
egen M_tot3=sum(M_quantity3), by(year month)
gen M_sh=M_quantity/M_tot
gen M_sh3=M_quantity3/M_tot3
gen M_shdif=M_sh3-M_sh

*price changes and market share changes
table substance if year==2008&month==12 [fw=Xtablets], contents(mean M_price_ch)
table substance if year==2008&month==12, contents(sum M_shdif)

table firmsubst if year==2008&month==12 [fw=Xtablets], contents(mean M_price_ch)
table firmsubst if year==2008&month==12, contents(sum M_shdif)

drop M_tot M_tot3 M_sh M_sh3 M_shdif M_price_ch



dddddddddddddddddddddd



****************************************************************************
* PART 5. EX POST ANALYSIS
****************************************************************************

* SUMSTATS IN TABLE 1, 2 AND 3
egen tot=sum(PX1),by(year)
gen PX2=PX1/tot

table form substance if year==2008,c(sum PX2) row col
table brand substance if year==2008,c(sum PX2) row col
sum PX1  Xtablets Xddd Xnormal Ptablets Pddd Pnormal marketing sickwomen sickmen GDPnom popwomen popmen if year<2009
sum num numg numhg numf numfg numfgh $iv



// define time period and merger period
gen short=0
replace short=1 if date>566
gen merger=date>592
gen merger1=date>591


* 5.0 COMPUTE MARGINAL COST AGAIN, BUT NOW FOR FULL PERIOD

replace firm =1 if firm==2&date>=592  //date=592 in May 2009 (assume in section 5.2 that AZT/GSK are one firm to compute costs)
if `ces'==1 {
mergersim init, nests(form substance) ces price(Ptablets) revenue(PX1) marketsize(BL`ces') firm(firm)
	ivreg2 M_ls $exogvar prod1-prod56 (M_lp M_lsjh M_lshg = num*), noconstant robust
}
if `ces'==0 {
mergersim init, nests(form substance) unit price(Ptablets) quantity(Xtablets) marketsize(BL`ces') firm(firm)
	ivreg2 M_ls $exogvar prod1-prod56 (Ptablets M_lsjh M_lshg = num*), noconstant robust
}
mergersim market , conduct(0)  
gen costs=M_costs //gen costs=M_costs*(packsize_new/packsize)^-0.3279819
gen CX=costs*Xtablets


* 5.1 EX POST REGRESSIONS, TABLE 7

replace firm=2 if brand==18 	//make Panodil back a GSK firm (so as before start of section 5.0)


// BY FIRM

preserve
keep if short

*adjustment 1: only keep main brands (available before and after
tab brand merger
keep if brand==1|brand==4|brand==5|brand==6|brand==8|	brand==9|brand==14|brand==15|brand==16|brand==18|	brand==21|brand==22
		//Alvedon,Albyl,Alindrin,Alka-seltzer,Aspirin,	Bamyl,Ibumetin,Ipren,Magnecyl,Panodyl,				Reliv,Treo
tab brand merger

*adjustment 2: assign Recip and Ellem brands to Meda (since both sold to Meda)
tab brand firm
	// Albyl (brand=4) sold by Recip (firm=6) and Meda (firm=15) --> Recip sold to Meda in 01-20008
	// Alindrin (brand=5) sold by Recip (firm=6) and Meda (firm=15) --> Recip sold to Meda in 01-20008
	// Bamyl (brand=9) sold by Ellem (firm=5) and Meda (firm=15) --> Ellem sold to Meda in 05-20008
tab firm merger
replace firm=15 if firm==5
replace firm=15 if firm==6
tab firm merger

*adjustment 3: redefine firm as firm*substance
	//redefine McNeil(Pfizer) and Meda as McNeil-ASA, McNeil-Ibu, Meda-ASA and Meda-Ibu
	//no other redefinitions needed, since other firms are only in one segment
replace firm=99 if brand==15 // make Ipren of McNeil separate "firm", i.e. McNeil-Ibu
replace firm=199 if brand==5 // make Alindrin of Meda separate "firm", i.e. Meda-Ibu
label define firmlabel 14 "McNeil-ASA"  15 "Meda-ASA", modify
label define firmlabel 99 "McNeil-Ibu" 199 "Meda-Ibu", add

local num 8
tab brand firm
	// the above table shows which brands are sold by the 8 "firms"

*do the collapse and regressions
	collapse (sum) Xtablets PX CX (mean) merger,by(year month firm)
	gen date = ym(year, month)
	format date %tm
	gen Ptab=PX/Xtablets
	gen lPtab=log(Ptab)
	gen Ctab=CX/Xtablets
	gen lCtab=log(Ctab)
	tab firm,gen(firm)
	forvalues i = 1/`num' {
		generate mergerf`i' = merger*firm`i'
	}
	egen tX=sum(Xtablets),by(year month)
	gen sX=Xtablets/tX
	gen lnX=log(Xtablets)

reg lPtab firm2-firm`num' mergerf1-mergerf`num',r
*reg lCtab firm2-firm`num' mergerf1-mergerf`num',r
*reg sX firm2-firm`num' mergerf1-mergerf`num',r
reg lnX firm2-firm`num' mergerf1-mergerf`num',r

gen PX1=PX/1000000	
table firm merger ,c(mean PX1)	
table firm merger ,c(mean Ptab)	

restore

// BY MONTH

preserve
keep if short

*do the collapse and regressions
	collapse (sum) Xtablets PX CX (mean) merger,by(year month)
	gen date = ym(year, month)
	format date %tm
	gen Ptab=PX/Xtablets
	gen lPtab=log(Ptab)
	gen Ctab=CX/Xtablets
	gen lCtab=log(Ctab)

	gen lnX=log(Xtablets)
	
reg lPtab merger,r
reg lCtab merger,r
reg lnX merger,r

restore

ddfdf

// BY SUBSTANCE

preserve
keep if short

*do the collapse and regressions
	collapse (sum) Xtablets PX CX (mean) merger,by(year month substance)
	gen date = ym(year, month)
	format date %tm
	gen Ptab=PX/Xtablets
	gen lPtab=log(Ptab)
	gen Ctab=CX/Xtablets
	gen lCtab=log(Ctab)
	tab substance,gen(substance)
	forvalues i = 1/3 {
		generate mergers`i' = merger*substance`i'
	}
	egen tX=sum(Xtablets),by(year month)
	gen sX=Xtablets/tX
	gen lnX=log(Xtablets)

reg lPtab substance2-substance3 mergers1-mergers3,r
*reg lCtab substance2-substance3 mergers1-mergers3,r
*reg sX substance2-substance3 mergers1-mergers3,r
reg lnX substance2-substance3 mergers1-mergers3,r

restore
reerer
// package size changes

preserve
*adjustment 1: only keep main brands (available before and after
keep if brand==1|brand==4|brand==5|brand==6|brand==8|	brand==9|brand==14|brand==15|brand==16|brand==18|	brand==21|brand==22
*adjustment 2: assign Recip and Ellem brands to Meda (since both sold to Meda)
replace firm=15 if firm==5
replace firm=15 if firm==6

table firm substance if year==2008 [pweight=PX1], c(median packsize median packsize_new)
table firm substance if year==2008 [pweight=PX1], c(mean packsize mean packsize_new)
table brand substance if year==2008 [pweight=PX1], c(mean packsize mean packsize_new)

* Bamyl example (Meda): product drop
*table product if brand==9 & year==2009&month==5, c(mean packsize sum PX1 mean Ptablets)
*table product if brand==9 & year==2009&month==7, c(mean packsize sum PX1 mean Ptablets)
* Treo example (McNeil): product drop (fully visible in month==10)
*table product if brand==22 & year==2009&month==5, c(mean packsize sum PX1 mean Ptablets)
*table product if brand==22 & year==2009&month==7, c(mean packsize sum PX1 mean Ptablets)
* Aspirin example (Bayer): no change in packsize, but price increase
*table product if brand==8 & year==2009&month==5, c(mean packsize sum PX1 mean Ptablets)
*table product if brand==8 & year==2009&month==7, c(mean packsize sum PX1 mean Ptablets)
restore

zzzzzzzzzzzz
* 5.2 GRAPHS

* STEP 1 CREATE VARIABLES

*prices (only tablets for now!)
generate PtabPara = Ptablets if substance ==1
generate PtabIbu = Ptablets if substance ==3
generate PtabAsa = Ptablets if substance ==5

*unit sales
generate XtabPara = Xtablets if substance ==1
generate XtabIbu = Xtablets if substance ==3
generate XtabAsa = Xtablets if substance ==5

*value sales
generate PXpara = PX if substance==1
generate PXibu = PX if substance==3
generate PXasa = PX if substance==5


* STEP 2 DO DESCRIPTIVES

keep if PX>100000


* STEP 3 COLLAPSE THE DATA TO SHOW EVOLUTION OF PRICES, SALES, #PRODUCTS AND #FIRMS IN GRAPHS

collapse (mean) PtabPara PtabAsa PtabIbu firm (sum) XtabPara XtabAsa XtabIbu PXpara PXasa PXibu, by(year month)


* 3.1 Aggregate regressions

gen date = ym(year, month)
format date %tm


gen PtabPara1=PXpara/XtabPara
gen PtabAsa1=PXasa/XtabAsa
gen PtabIbu1=PXibu/XtabIbu

gen lPtabPara1=log(PtabPara1)
gen lPtabAsa1=log(PtabAsa1)
gen lPtabIbu1=log(PtabIbu1)

gen sh_Para=PXpara/(PXpara+PXasa+PXibu)
gen sh_Asa=PXasa/(PXpara+PXasa+PXibu)
gen sh_Ibu=PXibu/(PXpara+PXasa+PXibu)

gen shX_Para=XtabPara/(XtabPara+XtabAsa+XtabIbu)
gen shX_Asa=XtabAsa/(XtabPara+XtabAsa+XtabIbu)
gen shX_Ibu=XtabIbu/(XtabPara+XtabAsa+XtabIbu)

gen short=0
replace short=1 if date>566

gen merger=date>592
gen merger1=date>591

tsset date

reg lPtabPara1 merger if short==1
reg lPtabAsa1 merger if short==1
reg lPtabIbu1 merger if short==1

reg shX_Para merger if short==1
reg shX_Asa merger if short==1
reg shX_Ibu merger if short==1


* 3.2 Graphs
	*date=591 refers to 04/2009
	*date=567 refers to 04/2007
	
*main figures for paper
label variable PtabPara "Paracetamol"
label variable PtabIbu "Ibuprofen"
label variable PtabAsa "ASA"
label variable date "Month"

gen shX_Para1=shX_Para*100
gen shX_Ibu1=shX_Ibu*100
gen shX_Asa1=shX_Asa*100
label variable shX_Para1 "Paracetamol"
label variable shX_Ibu1 "Ibuprofen"
label variable shX_Asa1 "ASA"

gen sh_Para1=sh_Para*100
gen sh_Ibu1=sh_Ibu*100
gen sh_Asa1=sh_Asa*100
label variable sh_Para1 "Paracetamol"
label variable sh_Ibu1 "Ibuprofen"
label variable sh_Asa1 "ASA"


* FIGURE 1 IN PAPER (PRICES)
line PtabPara PtabAsa PtabIbu date if short==1, title("Price evolution analgesics (April 2007-April 2011)") ytitle("Price (in Krone)") xline(591) note("Note: vertical line refers to the month of merger (April 2009)") saving(Figure2,replace)
* FIGURE 2 IN PAPER (SHARES IN VOLUMES)
line shX_Para1 shX_Asa1 shX_Ibu1 date if short==1, title("Market share evolution analgesics (April 2007-April 2011)") ytitle("Market share (in percent)") xline(591) note("Note: vertical line refers to the month of merger (April 2009)") saving(Figure3,replace)
* FIGURE NOT IN PAPER (SHARES IN VALUES)
line sh_Para1 sh_Asa1 sh_Ibu1 date if short==1, title("Market share evolution analgesics (April 2007-April 2011)") ytitle("Market share (in percent)") xline(591) note("Note: vertical line refers to the month of merger (April 2009)") saving(Figure3,replace)


overview firms & brands
2007
2008
2009
2010

AZT:
Alvedon, Curadon, Reliv
Alvedon, Curadon, Reliv
Alvedon, Reliv, Panodil
Alvedon, Reliv, Panodil

GSK:
Panodil
Panodil
Panodil
-

Nycomed:
Ibumetin
Ibumetin
Ibumetin
Ibumetin

Ellem (after correction):
Bamyl
Bamyl
-
-

Recip:
Albyl Alindrin
-
-
-

Bayer:
Alka Aspirin
Alka Aspirin
Alka Aspirin
Alka Aspirin

Pfizer:
Ipren Magnecyn Treo
Ipren Magnecyn Treo
Ipren Magnecyn Treo
Ipren Magnecyn Treo

Meda:
-
Albyl Alindrin Bamyl
Albyl Alindrin Bamyl
Albyl Alindrin Bamyl

CONCLUSION
- AZT&GSK merge in 2009, AZT removes Curadon in 2009
- Nycomed: no change
- Bayer: no change
- Pfizer: no change
- Ellem: sells Bamyl, but also done by Meda --> corrected (=takeover on April 1)
- Recip: sells Albyl&Alindrin, but from 2008 done by Meda (=takeover)














