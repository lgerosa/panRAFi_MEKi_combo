[PROB]
Author: Logan

Validated: Yes 

[PKMODEL]
cmt = "DEPOT CENT PERIPH", depot = TRUE

[GLOBAL]
#define CONC1  (CENT/(V2/1000))  

[PARAM] @annotated
 TVCL  :  3.28   :1 Clearance (L/day)
TVV2  :  174   :2 Volume of central compartment (L)
TVQ   :  4.6   :3 Intercompartmental clearance (L/day)
TVV3  :  69.2   :4 Volume of peripheral compartment (L)
TVKA  :  0.443  :5 First-order absorption rate constant (L/day)
TVF1  : 0.82     : typical bioavailablity 
DMAX  : -0.8  : 6 Reduction in bioavailablity dose  
D50   :  127   :7 D50 
TALAG : 0.502  : 8 Lag time (day)
PROP   : -0.203  :9 Prop error
ADD     :  109  :9 ADD error
ADOSE :   200   : 
BWT :   70 : 
BWTref : 70 : 
CLBWT  : 0.75: 
V2BWT  :  1   : 


E_V2i: 0 : 
E_CLi: 0 : 
E_Qi : 0 : 
E_V3i: 0 : 
E_KAi: 0 : 


[OMEGA] @annotated @block
ETA_CL   : 0.423        : 
ETA_V2   : 0.321  0.701 :

[OMEGA] @annotated
ETA_V3   : 2.73     :
ETA_KA   : 1.22  :
ETA_Q    : 0     :


[SIGMA] @annotated
EPS_1 : 1 : error term

[MAIN]

double CLCOV = pow((BWT/BWTref), CLBWT); 
double V2COV = pow((BWT/BWTref), V2BWT); 
double DOSEF =  1 + ((DMAX*ADOSE)/(D50+ADOSE)); 

// PK parameters
double CL   = TVCL * CLCOV * exp(ETA_CL + E_CLi);
double V2   = TVV2 * V2COV * exp(ETA_V2 + E_V2i);
double Q    = TVQ * exp(ETA_Q + E_Qi);
double V3   = TVV3 * exp(ETA_V3 + E_V3i);
double KA   = TVKA * exp(ETA_KA + E_KAi);
double ALAG = TALAG ; 
double F1   = DOSEF ;

// absorption lag time in depot compartment
ALAG_DEPOT = ALAG;
F_DEPOT = F1; 

[TABLE]
double IPRED = CONC1; 
double W = sqrt(pow(PROP,2)*pow(IPRED,2) + pow(ADD,2));  
double DV = IPRED + EPS_1*W;

[CAPTURE] @annotated
DV    : Concentration (uM) with residual variability
IPRED : Concentration (uM) without residual variability