/* Import the dataset */
FILENAME REFFILE '/home/u63745816/sasuser.v94/BIST P8110 Applied Regression II/Midterm Project/midtermdata.csv';

PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=midterm;
	GETNAMES=YES;
RUN;

PROC CONTENTS DATA=midterm; RUN;


/*  Preparation: Combine those categories for race, IV drug use, and Karnofsky score. */
data midtermCombine;
	set midterm;
	if raceth in (3, 4, 5, 6) then racethCombined = 3; 
	else racethCombined = raceth;
	
	if ivdrug in (2, 3) then ivdrugCombined = 2; 
	else ivdrugCombined = ivdrug;
	
	if karnof in (70, 80) then karnofCombined = 80; 
	else karnofCombined = karnof;
run;

/* Define formats for categorical variables */
proc format;
	value censor_fmt
		0 = 'Otherwise'
		1 = 'AIDS defining diagnosis or death';
		
	value censord_fmt
		0 = 'Otherwise'
		1 = 'Death';
		
    value tx_fmt
        0 = 'Two-Drug Regimen'
        1 = 'Three-Drug Regimen';
    
    value sex_fmt
        1 = 'Male'
        2 = 'Female';
    
    value raceth_fmt
        1 = 'White Non-Hispanic'
        2 = 'Black Non-Hispanic'
        3 = 'Other';
    
    value ivdrug_fmt
        1 = 'Never'
        2 = 'Current or Previous';
    
    value hemophil_fmt
        1 = 'Yes'
        0 = 'No';
    
    value karnof_fmt
        100 = 'Normal'
        90 = 'Minor symptoms'
        80 = 'Limited normal activity';
    
    value strat2_fmt
        0 = 'CD4 ≤ 50'
        1 = 'CD4 > 50';
run;

/* Apply formats to the variables */
data midtermCombined;
    set midtermCombine;
    format censor censor_fmt. 
    		censor_d censord_fmt.
    		tx tx_fmt.
    		sex sex_fmt.
    		racethCombined raceth_fmt.
    		ivdrugCombined ivdrug_fmt.
    		hemophil hemophil_fmt.
    		karnofCombined karnof_fmt.
    		strat2 strat2_fmt.;
run;

/* Step 1: Summarize the baseline characteristics */
	* Continuous variables for all patients;
proc means data=midtermCombine n mean median std Q1 Q3;
	var age cd4 priorzdv time time_d;
	title "Baseline Characteristics Summary - Continuous Variables (All Patients)";
run;
	* Categorical variables summary for all patients;
proc freq data=midtermCombine;
	tables sex racethCombined ivdrugCombined hemophil karnofCombined strat2 censor censor_d/ missing;
    title "Baseline Characteristics Summary - Categorical Variables (All Patients)";
run;

	* Summarize the baseline characteristics by treatment regimen;
proc means data=midtermCombine n mean median std Q1 Q3;
    class tx;
    var age cd4 priorzdv time time_d;
    title "Baseline Characteristics Summary - Continuous Variables by Treatment Group";
run;
	*Categorical variables summary by treatment group;
proc freq data=midtermCombine;
    tables tx*(sex racethCombined ivdrugCombined hemophil karnofCombined strat2 censor censor_d) / missing;
    title "Baseline Characteristics Summary - Categorical Variables by Treatment Group";
run;

	* Perform statistics tests;
	** Wilcoxon rank sum test for continuous variables by treatment group;
proc npar1way data=midtermCombine wilcoxon;
    class tx;
    var age cd4 priorzdv time time_d;
    title "Wilcoxon Rank Sum Test for Continuous Variables by Treatment Group";
run;

	** Fisher’s exact test for smaller categorical variables;
proc freq data=midtermCombine;
    tables tx*sex / fisher;
    tables tx*hemophil / fisher;
    tables tx*censor / fisher;
    tables tx*censor_d / fisher;
    title "Fisher’s Exact Test for Categorical Variables by Treatment Group";
run;

	** Pearson’s Chi-Square test for larger categorical variables;
proc freq data=midtermCombine;
    tables tx*(racethCombined ivdrugCombined karnofCombined strat2) / chisq;
    title "Pearson’s Chi-Square Test for Categorical Variables by Treatment Group";
run;

/* Step 2: Kaplan-Meier survival curves and log-rank tests to compare PFS and OS 
between the two treatment regimens */
ods graphics on;
proc lifetest data=midtermcombine plots=survival(atrisk cb) notable;
	time time*censor(0);
	strata tx / test = logrank;
	survival out = km_pfs;
	format tx tx_fmt.;
	title "Kaplan-Meier Survival Curves and Log-Rank Test for PFS";
run;
ods graphics off;

ods graphics on;
proc lifetest data=midtermCombine plots=survival(atrisk cb) notable;
    time time_d*censor_d(0);
    strata tx / test = logrank;
    survival out = km_os;
    format tx tx_fmt.;
    title "Kaplan-Meier Survival Curves and Log-Rank Test for OS";
run;
ods graphics off;

/* Univariable analysis for PFS (Table 2) */
proc phreg data=midtermCombine;
    class tx(ref='0') / param=ref;
    model time*censor(0) = tx / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for PFS (Treatment)";
run;

proc phreg data=midtermCombine;
    class strat2(ref='0') / param=ref;
    model time*censor(0) = strat2 / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for PFS (CD4 stratum at screening)";
run;

proc phreg data=midtermCombine;
    class sex(ref='1') / param=ref;
    model time*censor(0) = sex / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for PFS (Gender)";
run;

proc phreg data=midtermCombine;
    class racethCombined(ref='3') / param=ref;
    model time*censor(0) = racethCombined / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for PFS (Race)";
run;

proc phreg data=midtermCombine;
    class ivdrugCombined(ref='1') / param=ref;
    model time*censor(0) = ivdrugCombined / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for PFS (IV Drug Use)";
run;

proc phreg data=midtermCombine;
    class hemophil(ref='0') / param=ref;
    model time*censor(0) = hemophil / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for PFS (Hemophilia)";
run;

proc phreg data=midtermCombine;
    class karnofCombined(ref='80') / param=ref;
    model time*censor(0) = karnofCombined / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for PFS (Karnofsky Score)";
run;

proc phreg data=midtermCombine;
    model time*censor(0) = age / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for PFS (Age)";
run;

proc phreg data=midtermCombine;
    model time*censor(0) = cd4 / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for PFS (CD4)";
run;

proc phreg data=midtermCombine;
    model time*censor(0) = priorzdv / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for PFS (priorzdv)";
run;

/* Univariable analysis for OS (Table 3) */
proc phreg data=midtermCombine;
    class tx(ref='0') / param=ref;
    model time_d*censor_d(0) = tx / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for OS (Treatment)";
run;

proc phreg data=midtermCombine;
    class strat2(ref='0') / param=ref;
    model time_d*censor_d(0) = strat2 / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for OS (CD4 stratum at screening)";
run;

proc phreg data=midtermCombine;
    class sex(ref='1') / param=ref;
    model time_d*censor_d(0) = sex / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for OS (Gender)";
run;

proc phreg data=midtermCombine;
    class racethCombined(ref='3') / param=ref;
    model time_d*censor_d(0) = racethCombined / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for OS (Race)";
run;

proc phreg data=midtermCombine;
    class ivdrugCombined(ref='1') / param=ref;
    model time_d*censor_d(0) = ivdrugCombined / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for OS (IV Drug Use)";
run;

proc phreg data=midtermCombine;
    class hemophil(ref='0') / param=ref;
    model time_d*censor_d(0) = hemophil / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for OS (Hemophilia)";
run;

proc phreg data=midtermCombine;
    class karnofCombined(ref='80') / param=ref;
    model time_d*censor_d(0) = karnofCombined / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for OS (Karnofsky Score)";
run;

proc phreg data=midtermCombine;
    model time_d*censor_d(0) = age / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for OS (Age)";
run;

proc phreg data=midtermCombine;
    model time_d*censor_d(0) = cd4 / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for OS (CD4)";
run;

proc phreg data=midtermCombine;
    model time_d*censor_d(0) = priorzdv / risklimits covb ties=efron;
    title "Univariate Cox Proportional Hazards for OS (priorzdv)";
run;

/* Multivariable analysis for PFS (Table 4) */
proc phreg data=midtermCombine;
    class tx(ref='0') sex(ref='1') karnofCombined(ref='80') strat2(ref='0') / param=ref;
    model time*censor(0) = tx age sex karnofCombined strat2 / risklimits covb ties=efron;
    title "Multivariable Analysis for PFS(2) (Table 4)";
run;

/* Multivariable analysis for OS (Table 5) */
proc phreg data=midtermCombine;
    class tx(ref='0') sex(ref='1') karnofCombined(ref='80') strat2(ref='0') / param=ref;
    model time_d*censor_d(0) = tx age sex karnofCombined strat2 / risklimits covb ties=efron;
    title "Multivariable Analysis for OS(2) (Table 5)";
run;

/* Test effect modifier between between the treatment regimens and PFS */
proc phreg data=midtermCombine;
    class tx(ref='0') sex(ref='1') karnofCombined(ref='80') strat2(ref='0') / param=ref;
    model time*censor(0) = tx age sex karnofCombined strat2 tx*age/ risklimits covb ties=efron;
    title "Effect Modifier 'Age' between treatment regimens and PFS";
run;

proc phreg data=midtermCombine;
    class tx(ref='0') sex(ref='1') karnofCombined(ref='80') strat2(ref='0') / param=ref;
    model time*censor(0) = tx age sex karnofCombined strat2 tx*sex/ risklimits covb ties=efron;
    title "Effect Modifier 'Sex' between treatment regimens and PFS";
run;

/* Test effect modifier between between the treatment regimens and OS */
proc phreg data=midtermcombine;
	class tx(ref = '0') sex(ref='1') karnofCombined(ref='80') strat2(ref='0') / param=ref;
    model time_d*censor_d(0) = tx age sex karnofCombined strat2 tx*age / risklimits covb ties=efron;
    title "Effect Modifier 'Age' between treatment regimens and PFS";
run;

proc phreg data=midtermcombine;
	class tx(ref = '0') sex(ref='1') karnofCombined(ref='80') strat2(ref='0') / param=ref;
    model time_d*censor_d(0) = tx age sex karnofCombined strat2 tx*sex / risklimits covb ties=efron;
    title "Effect Modifier 'Sex' between treatment regimens and OS";
run;


ods graphics on;
proc phreg data=midtermCombine;
	class tx(ref = '0') sex(ref='1') karnofCombined(ref='80') strat2(ref='0') / param=ref;
    model time*censor(0) = tx age sex karnofCombined strat2 / risklimits covb ties=efron;
    assess PH / resample;
    title "Test the PH assumption of treatment regimens and PFS";
run;
ods graphics off;

ods graphics on;
proc phreg data=midtermCombine;
	class tx(ref = '0') sex(ref='1') karnofCombined(ref='80') strat2(ref='0') / param=ref;
    model time_d*censor_d(0) = tx age sex karnofCombined strat2 / risklimits covb ties=efron;
    assess PH / resample;
    title "Test the PH assumption of treatment regimens and OS";
run;
ods graphics off;