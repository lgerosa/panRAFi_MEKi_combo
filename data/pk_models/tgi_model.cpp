[PROB]
Author: Logan

Time unit: day (generic)
Volume units: L (generic)
Validated: Yes

[CMT] @annotated
DUMMY : dummy compartment for mrgsolve to run

[PARAM] @annotated
TVKG     :   0.0130  : Tumor growth rate (1/week)
TVKS     :   0.0346  : Tumor shrinkage rate (1/week)
TVBSL    :   65.1    : Baseline tumor burden (mm)
EPSVAL   :   42.3    : Residual (mm)

[OMEGA] @annotated 
ETA_KG    :  0                      : 1
ETA_KS    :  0              : 2
ETA_BSL   :  0   : 3

[SIGMA] @annotated
ADD : 1 : Additive error

[MAIN]


// TGI parameters
double KG    = TVKG *  exp(ETA_KG);
double KS    = TVKS  * exp(ETA_KS);
double BSL   = TVBSL * exp(ETA_BSL);

double TT  = 0;  // treatment beginning

double TTG    = (log(KS) - log(KG))/(KG+KS) + TT;
double TR6    = exp(KG*(6-TT))  + exp(-KS*(6-TT))-1;
double TR12   = exp(KG*(12-TT)) + exp(-KS*(12-TT))-1;
double TR18   = exp(KG*(18-TT)) + exp(-KS*(18-TT))-1;
double TR24   = exp(KG*(24-TT)) + exp(-KS*(24-TT))-1;

[ODE]
dxdt_DUMMY = 0;

[TABLE]
double IPRED = 0;

if (TIME <= TT) {
	IPRED = BSL * exp(KG*(TIME-TT));
}

if (TIME > TT) {
	IPRED = BSL * (exp(KG*(TIME-TT)) + exp(-KS*(TIME-TT)) - 1);
}

double DV = IPRED + sqrt(EPSVAL) * ADD; 

[CAPTURE] @annotated
 IPRED : Concentration without residual variability
 DV    : Concentration with residual variability
 KG : Tumor growth rate (1/week)
 KS : Tumor shrinkage rate (1/week)  
 BSL : Baseline tumor burden (mm)