**** JORDAGK.DO

*** May 31, 2016
*** This program estimates IRFs and both 2-step and 1-step multipliers using the Jorda method

*** Uses Gordon-Krenn (GK) transformation (i.e. variables are divided by potential)

***
*** Requires:
***     rzdat.xlsx, updated April 7, 2016
********************************************************************************

 #delimit;

drop _all;
clear all;

set more 1;
set matsize 800;

capture log close;
log using jordagkirfs_results.log, replace;


/*******************************************************************************
  SET PARAMETERS THAT GOVERN SPECIFICATION
*******************************************************************************/

local sample = 1;  /*1 = full sample, 2 = post-WWII */

local omit nomit;  /*either nomit(don't omit subsample) or wwii (omit WWII) */

local state fstate;  /* slack or zlb or recession or ag*/

local shock newsy; /* shock identification: either newsy or bp */

local p = 4; /*number of lags of control variables*/

local trends = 0; /*0 = no trends, 1 = trends */

local tax = 0; /*0 = exclude taxes, 1 = include taxes */


*******************************************************************************;
** RAW DATA IMPORTATION AND DATA SETUP;
*******************************************************************************;

*import excel rzdat.xlsx, sheet("rzdat") firstrow;

insheet using rz_dat_updated.csv;

drop if quarter<1889;

gen qdate = q(1889q1) + _n-1;
tsset qdate, q;

/* World War II rationing sample.  Start in 1941q3 because of constraints on
     auto production.  Most rationing ended in Aug/Sept 1945, a few items in
	 Nov/Dec 1945 */

gen wwii = quarter>=1941.5 & quarter<1946;

gen nomit = 0;  /* indicator for no omit */


*** DEFINE QUARTIC TREND;

gen t = _n;
gen t2 = t^2;
gen t3 = t^3;
gen t4 = t^4;

*** DEFINE STATE VARIABLE;

gen slack = unemp >= 6.5;  /* unemployment state with fixed threshold */
gen slack8 = unemp>=8;
gen slackhp = unemp>=hpunemp_split;  /* unemployment state with hp threshold */

gen zlb = zlb_dummy;  /*  zlb state */
gen zlb5 = tbill<=0.5; /*tbill rate less than 0.5*/
gen both =unemp >= 6.5 & zlb_dummy>=1; /* consider slack and ZLB state together*/

*** NORMALIZATION;

/* choice of potential GDP for normalization:

   rgdp_potCBO (cubic trend early, CBO late) or rgdp_pott6 (6th degree for full),
   both fitted excluding Great Depression, WWII:  quarter>=1930 & quarter<1947*/

local ynorm rgdp_pott6; /* rgdp_pott6 or rgdp_potcbo */

* BASIC VARIABLES;

gen newsy = news/(L.`ynorm'*L.pgdp);
gen rgov = ngov/pgdp;
gen rtax = nfedcurrreceipts_nipa/pgdp;
gen taxy = nfedcurrreceipts_nipa/ngdp;
gen debty = pubfeddebt_treas/L.ngdp;
gen lpgdp = ln(pgdp);
gen ly = ln(rgdp);

gen infl = 400*D.lpgdp;

* normalize variables and shorten names;

gen y = rgdp/`ynorm';
gen g = rgov/`ynorm';
 
gen bp = g; /* Blanchard-Perotti shock is just orthogonalized current g */

*** AG DEFINITION OF STATE:  ag = 1 is extreme recession, ag = 0 is extreme expansion;

gen z = 100*(F3.ly - L4.ly)/7;  /* AG definition of state */

*The mean of z is approx. 0.8 and std is 0.5.  AG specically use those numbers so we do too;

gen znorm = (z - 0.8)/0.5;

gen fznorm = exp(-1.5*znorm)/(1 + exp(-1.5*znorm));

gen ag = fznorm;

*******************************************************************************;
** CUMULATIVE VARIABLES;
*******************************************************************************;

gen cumuly = 0;
gen cumulg = 0;
 
forvalues i = 0/20 {;

   gen f`i'cumuly = F`i'.y + cumuly;
   gen f`i'cumulg = F`i'.g + cumulg;
   
   gen recf`i'cumulg = f`i'cumulg*L.`state';
   gen expf`i'cumulg = f`i'cumulg*(1-L.`state');
   
   replace cumuly = f`i'cumuly;
   replace cumulg = f`i'cumulg;
   
};


*******************************************************************************;
**  INTERACTION OF SHOCKS WITH STATE;
*******************************************************************************;

 foreach var in newsy bp { ;
 
   gen rec0`var' = `var'*L.`state';
   gen exp0`var' = `var'*(1-L.`state');
 
 };

*******************************************************************************;
** CREATE LISTS;
*******************************************************************************;

   if `sample'==1 {;

         gen h = t - 1;  /* h is the horizon for the irfs */
         global trendlist t t2 t3 t4;
   };

    else {;
         drop if quarter<1947;
         gen h = t - 232 - 1;
         global trendlist t t2;
     };
	 
forvalues i = 1/`p' {; 

  foreach var in newsy y g taxy debty infl{;

    gen rec`var'`i' = L`i'.`var'*L.`state';
    gen exp`var'`i' = L`i'.`var'*(1-L.`state');
 
  };
};

  if `trends'==0 {;
  
    if `tax'==0 {;
  
      global newsylinxlist L(1/`p').newsy L(1/`p').y L(1/`p').g ;
      global bplinxlist L(1/`p').y L(1/`p').g ;
	  global newsynlxlist L.`state' recnewsy? expnewsy? recy? expy? recg? expg? ;
      global bpnlxlist L.`state' recy? expy? recg? expg? ;
	
	};
	
	else {;
	
	  global newsylinxlist L(1/`p').newsy L(1/`p').y L(1/`p').g L(1/`p').taxy; /*L(1/`p').infl;*/
      global bplinxlist L(1/`p').y L(1/`p').g L(1/`p').taxy; /* L(1/`p').infl;*/
	  global newsynlxlist L.`state' recnewsy? expnewsy? recy? expy? recg? expg? rectaxy? exptaxy?; /* expinfl? recinfl?;*/
      global bpnlxlist L.`state' recy? expy? recg? expg? rectaxy? exptaxy?; /* expinfl? recinfl?;*/
	
    };
  };
  
  else {;
  
    if `tax'==0 {;
	
      global newsylinxlist L(1/`p').newsy L(1/`p').y L(1/`p').g $trendlist;
      global bplinxlist L(1/`p').y L(1/`p').g $trendlist;
      global newsynlxlist L.`state' recnewsy? expnewsy? recy? expy? recg? expg? $trendlist;
      global bpnlxlist L.`state' recy? expy? recg? expg? $trendlist;
	
	};
	
	else {;
	
	global newsylinxlist L(1/`p').newsy L(1/`p').y L(1/`p').g L(1/`p').taxy $trendlist;
    global bplinxlist L(1/`p').y L(1/`p').g L(1/`p').taxy $trendlist;
	global newsynlxlist L.`state' recnewsy? expnewsy? recy? expy? recg? expg? rectaxy? exptaxy? $trendlist;
    global bpnlxlist L.`state' recy? expy? recg? expg? rectaxy? exptaxy? $trendlist;
	
    };
	
};


global newsylinshock newsy;
global newsynlshock rec0newsy exp0newsy;

global bplinshock bp;
global bpnlshock rec0bp exp0bp;


** INITIALIZE SUM OF EFFECTS TO 0 AND PARAMETERS SERIES TO MISSING;

gen sumliny = 0; gen sumling = 0;
gen sumexpy = 0; gen sumexpg = 0;
gen sumrecy = 0; gen sumrecg = 0;

foreach var in bylin byexp byrec bglin bgexp bgrec up95bylin up95byexp up95byrec up95bglin up95bgexp up95bgrec
  lo95bylin lo95byexp lo95byrec lo95bglin lo95bgexp lo95bgrec seylin seyexp seyrec seglin segexp segrec
  multlin multexp multrec {;
  
  quietly gen `var' = .;
  
}; 


*******************************************************************************;
** ESTIMATION OF IRFS
*******************************************************************************;

forvalues i = 0/0 {; /*Must treat horizon = 0 different in case the shock is BP*/

  ivreg2 F`i'.y $`shock'linshock $`shock'linxlist if L`p'.`omit'==0 & F`i'.`omit'==0 & `omit'==0, robust bw(auto);  

  gen bylinh`i' = _b[$`shock'linshock];
  
  gen seylinh`i' = _se[$`shock'linshock];
  
  ivreg2 F`i'.y exp0`shock' rec0`shock' $`shock'nlxlist if L`p'.`omit'==0 & F`i'.`omit'==0 & `omit'==0, robust bw(auto);	

  gen byexph`i' = _b[exp0`shock'];
  gen byrech`i' = _b[rec0`shock'];
  
  gen seyexph`i' = _se[exp0`shock'];
  gen seyrech`i' = _se[rec0`shock']; 
  
  
  if "`shock'" == "bp" {; 
  
  gen bglinh`i' = 1; 
  gen seglinh`i' = 0; 
  gen bgexph`i' = 1; 
  gen bgrech`i' = 1;
  gen segexph`i' = 0;
  gen segrech`i' = 0;
};

else {;

  *ivreg2 F`i'.g $`shock'linshock $`shock'linxlist if L`p'.`omit'==0 & F`i'.`omit'==0 & `omit'==0, robust bw(auto);
  ivreg2 F`i'.g $`shock'linshock $`shock'linxlist if L`p'.`omit'==0 & F`i'.`omit'==0 & `omit'==0, robust bw(auto);

  gen bglinh`i' = _b[$`shock'linshock];  
  gen seglinh`i' = _se[$`shock'linshock]; 

  ivreg2 F`i'.g exp0`shock' rec0`shock' $`shock'nlxlist if L`p'.`omit'==0 & F`i'.`omit'==0 & `omit'==0, robust bw(auto);	

  gen bgexph`i' = _b[exp0`shock']; 
  gen bgrech`i' = _b[rec0`shock'];
  
  gen segexph`i' = _se[exp0`shock'];
  gen segrech`i' = _se[rec0`shock']; 

};

replace sumliny = bylinh`i' + sumliny;
  replace sumling = bglinh`i' + sumling;
  
  replace sumexpy = byexph`i' + sumexpy;
  replace sumexpg = bgexph`i' + sumexpg;
  
  replace sumrecy = byrech`i' + sumrecy;
  replace sumrecg = bgrech`i' + sumrecg;
  
   gen multlinh`i' = sumliny/sumling;
   gen multexph`i' = sumexpy/sumexpg;
   gen multrech`i' = sumrecy/sumrecg;
  
  foreach var in bylin byexp byrec bglin bgexp bgrec multlin multexp multrec {;
  
    quietly replace `var' = `var'h`i' if h==`i';
	
  };
  
  foreach var in ylin glin yexp gexp yrec grec {;
  
    quietly replace up95b`var' = b`var'h`i' + 1.96*se`var'h`i' if h==`i';
	quietly replace lo95b`var' = b`var'h`i' - 1.96*se`var'h`i' if h==`i';
	
  };

};



forvalues i = 1/20 {;


ivreg2 F`i'.y $`shock'linshock $`shock'linxlist if L`p'.`omit'==0 & F`i'.`omit'==0 & `omit'==0, robust bw(auto);

  gen bylinh`i' = _b[$`shock'linshock];  
  gen seylinh`i' = _se[$`shock'linshock];
  
ivreg2 F`i'.g $`shock'linshock $`shock'linxlist if L`p'.`omit'==0 & F`i'.`omit'==0 & `omit'==0, robust bw(auto);

  gen bglinh`i' = _b[$`shock'linshock];  
  gen seglinh`i' = _se[$`shock'linshock]; 

ivreg2 F`i'.y exp0`shock' rec0`shock' $`shock'nlxlist if L`p'.`omit'==0 & F`i'.`omit'==0 & `omit'==0, robust bw(auto);	

  gen byexph`i' = _b[exp0`shock'];
  gen byrech`i' = _b[rec0`shock'];
  
  gen seyexph`i' = _se[exp0`shock'];
  gen seyrech`i' = _se[rec0`shock']; 

ivreg2 F`i'.g exp0`shock' rec0`shock' $`shock'nlxlist if L`p'.`omit'==0 & F`i'.`omit'==0 & `omit'==0, robust bw(auto);	

  gen bgexph`i' = _b[exp0`shock']; 
  gen bgrech`i' = _b[rec0`shock'];
  
  gen segexph`i' = _se[exp0`shock'];
  gen segrech`i' = _se[rec0`shock'];
  
  
  replace sumliny = bylinh`i' + sumliny;
  replace sumling = bglinh`i' + sumling;
  
  replace sumexpy = byexph`i' + sumexpy;
  replace sumexpg = bgexph`i' + sumexpg;
  
  replace sumrecy = byrech`i' + sumrecy;
  replace sumrecg = bgrech`i' + sumrecg;
  
   gen multlinh`i' = sumliny/sumling;
   gen multexph`i' = sumexpy/sumexpg;
   gen multrech`i' = sumrecy/sumrecg;
  
  foreach var in bylin byexp byrec bglin bgexp bgrec multlin multexp multrec seyexp seyrec segexp segrec {;
  
    quietly replace `var' = `var'h`i' if h==`i';
	
  };
  
  foreach var in ylin glin yexp gexp yrec grec {;
  
    quietly replace up95b`var' = b`var'h`i' + 1.96*se`var'h`i' if h==`i';
	quietly replace lo95b`var' = b`var'h`i' - 1.96*se`var'h`i' if h==`i';
	
  };

  
};


display as text "MULTIPLIERS:  2 STEP";

rename multlin multlin2;
rename multexp multexp2;
rename multrec multrec2;

outsheet h multlin2 multexp2 multrec2 using junk2step.csv if h<=20, comma replace;

outsheet h byexp byrec bgexp bgrec seyexp seyrec segexp segrec using junk2stepirfs.csv if h<=20, comma replace;

label var bglin "Gov, linear model";
label var bylin "GDP, linear model";
label var bgexp "GOV, expansion";
label var byexp "GDP, expansion";
label var bgrec "Gov, recession";
label var byrec "GDP, recession";


tw (rarea up95bglin lo95bglin h, bcolor(gs12) clw(medthin medthin)) 
  (scatter bglin h, c(l ) clp(l ) ms(i ) clc(black) mc(black) clw(medthick)) if h<=20,
  saving(junkg.gph,replace);

tw (rarea up95bylin lo95bylin h, bcolor(gs12) clw(medthin medthin)) 
  (scatter bylin h, c(l ) clp(l ) ms(i ) clc(black) mc(black) clw(medthick)) if h<=20,
  saving(junky.gph,replace);
  
  
 tw (rarea up95bgrec lo95bgrec h, bcolor(gs12) clw(medthin medthin))
    (scatter up95bgexp bgexp lo95bgexp bgrec h, clw(medthin medthick medthin medthick)
  c(l l l l l) clp(- l - l) clc(red red red blue ) ms(i o i i i i) mc(red red red blue)) if h<=20,
  saving(junkgnl.gph,replace);

tw (rarea up95byrec lo95byrec h, bcolor(gs12) clw(medthin medthin))
    (scatter up95byexp byexp lo95byexp byrec h,  clw(medthin medthick medthin medthick)
  c(l l l l l) clp(- l - l) clc(red red red blue ) ms(i o i i i i) mc(red red red blue)) if h<=20,
  saving(junkynl.gph,replace);
  
graph combine junkg.gph junky.gph junkgnl.gph junkynl.gph, col(2) iscale(0.5);

*****************************************************************************;

drop mult???h*;

drop sey*;


foreach var in multlin1 multexp1 multrec1 Fkplin Fkpexp Fkprec seylin seyexp seyrec ptestdiff Fdifflin Fdiffexp Fdiffrec{;
  
  quietly gen `var' = .;
  
}; 
*******************************************************************************;
** ESTIMATION OF CUMULATIVE;
*******************************************************************************;

forvalues i = 0/20 {; 

  ivreg2 f`i'cumuly (f`i'cumulg = $`shock'linshock) $`shock'linxlist if L`p'.`omit'==0 & F`i'.`omit'==0 & `omit'==0, robust bw(auto);
  gen Fkplinh`i'= e(widstat); /* Kleibergen-Paap rk Wald F statistic*/
  gen Fdifflinh`i'= Fkplinh`i'- 23.1085; 
  gen multlinh`i' = _b[f`i'cumulg];
  gen seylinh`i' = _se[f`i'cumulg]; /* HAC robust standard error*/
  
  ivreg2 f`i'cumuly (expf`i'cumulg = exp0`shock') $`shock'nlxlist if L`p'.`omit'==0 & F`i'.`omit'==0 & `omit'==0, robust bw(auto);
  gen Fkpexph`i'= e(widstat);
  gen Fdiffexph`i'= Fkpexph`i'- 23.1085; 
  gen multexph`i' = _b[expf`i'cumulg];
  gen seyexph`i' = _se[expf`i'cumulg];
  
  ivreg2 f`i'cumuly (recf`i'cumulg = rec0`shock') $`shock'nlxlist if L`p'.`omit'==0 & F`i'.`omit'==0 & `omit'==0, robust bw(auto);
  gen Fkprech`i'= e(widstat);  
  gen Fdiffrech`i'= Fkprech`i'- 23.1085; 
  gen multrech`i' = _b[recf`i'cumulg];
  gen seyrech`i' = _se[recf`i'cumulg];
  
  ivreg2 f`i'cumuly (expf`i'cumulg recf`i'cumulg = exp0`shock' rec0`shock') $`shock'nlxlist if L`p'.`omit'==0 & F`i'.`omit'==0 & `omit'==0, robust bw(auto);
  test expf`i'cumulg=recf`i'cumulg;	
  gen ptestdiffh`i' = r(p);
  
  	
 foreach var in multlin multexp multrec {;
  
    quietly replace `var'1 = `var'h`i' if h==`i';
	
  };
  
  foreach var in seylin seyexp seyrec ptestdiff Fkplin Fkpexp Fkprec  {;
  
    quietly replace `var' = `var'h`i' if h==`i';
	
  };
  
 foreach var in  Fdifflin Fdiffexp Fdiffrec {;
  
    quietly replace `var' = `var'h`i' if h==`i';
	quietly replace `var' = 30 if `var'>30;
	
  };
};


display as text "First stage F-statistic (Kleibergen-Paap rk Wald F statistic): Linear, Expansion, Recession";

list h Fkplin Fkpexp Fkprec if h<=20;
outsheet h Fkplin Fkpexp Fkprec using junk.csv if h<=20, comma replace ;
outsheet h Fdifflin Fdiffexp Fdiffrec using junkfdiff.csv if h<=20, comma replace ;

display as text "Multipliers from 1 step and 3 step approaches: Linear, Expansion, Recession";

list h multlin1 multlin2 multexp1 multexp2 multrec1 multrec2 if h<=20;
outsheet h multlin1 multexp1 multrec1 using junk1step.csv if h<=20, comma replace ;

display as text "Multipliers and corresponding standrad errors: Linear, Expansion, Recession";

list h multlin1 seylin multexp1 seyexp multrec1 seyrec ptestdiff if h<=20;
outsheet h multlin1 seylin multexp1 seyexp multrec1 seyrec ptestdiff using junkmultse.csv if h<=20, comma replace ;

/* used for understanding why 2-step different from 1-step - we now know the answer is different samples
foreach var in bglin bylin bgexp byexp bgrec byrec {;
  gen sum`var' = sum(`var');
};

list h bgrec byrec sumbgrec sumbyrec multrec2 multrec1 if h<=20;*/


capture log close;

 
