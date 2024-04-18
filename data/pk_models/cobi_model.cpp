[PROB]
Author: Jenny Nguyen, 2018
Source: 101.ctl

Time unit: day
Dose units: mg
Volume units: mL
Validated: Yes (Logan 2019)

[PKMODEL]
cmt = "DEPOT CENT PERIPH", depot = TRUE

[PARAM] @annotated
TVKA :  35.6    : First-order absorption rate constant (L/day)
TVCL : 322      : Clearance (L/day)
TVV2 : 511      : Volume of central compartment (L)
TVV3 : 295      : Volume of peripheral compartment (L)
TVQ  : 210      : Intercompartmental clearance (L/day)
ALAG :   0.0199 : Lag time (day)
TAGE :  -0.217  : Effect of age on TVCL
TBWT :   0.795  : Effect of body weight on TVV2
TVF1 :   1      : Relative bioavailability
AGE  :  57      : Typical individual age (yrs)
BWT  :  80      : Typical individual body weight (kg)
OCC  :   1      : Occasion (1 - 3)

[OMEGA] @annotated @block @correlation @name IIV 
// off-diagonal elements as correlations; diagnal elements as variances
ETA_KA : 2.76                             : ETA on KA 
ETA_CL :  0   0.342                       : ETA on CL
ETA_V2 :  0   0.882   0.237               : ETA on V2
ETA_V3 :  0   0.381  -0.066  0.631        : ETA on V3
ETA_Q  :  0     0       0      0    0.808 : ETA on Q

[OMEGA] @annotated @name IOV 
// diagonal matrix is assumed
ETA_IOV1 : 0.211 : IOV for OCC1
ETA_IOV2 : 0.211 : IOV for OCC2
ETA_IOV3 : 0.211 : IOV for OCC3

[SIGMA] @annotated
// additive SD = 0.423 --> var = 0.1866
ADD : 0.187 : Additive error on log concentration 

[MAIN]
// effect of covariates on parameters
double CLAGE = pow((AGE/57), TAGE);
double V2BWT = pow((BWT/80), TBWT); 

// PK parameters

// define interoccasions variability
double OCC1 = 0;
double OCC2 = 0;
double OCC3 = 0; 

if (OCC == 1) OCC1 = 1;
if (OCC == 2) OCC2 = 1; 
if (OCC == 3) OCC3 = 1;

double ETA_IOV = ETA_IOV1*OCC1 + ETA_IOV2*OCC2 + ETA_IOV3*OCC3;

double KA    = TVKA * exp(ETA_KA);
double CL    = TVCL * CLAGE * exp(ETA_CL); 
double V2    = TVV2 * V2BWT * exp(ETA_V2);
double V3    = TVV3 * exp(ETA_V3);
double Q     = TVQ  * exp(ETA_Q);
double F1    = TVF1 * exp(ETA_IOV); 

// depot compartment lag time
F_DEPOT    = F1;
ALAG_DEPOT = ALAG;

[TABLE]
double IPRED = log(CENT / (V2 / 1000)); // dose = mg, V2 = L, conc = ng/mL on log scale
double DV    = IPRED + ADD;
double IPREDnormal = exp(IPRED);
double DVnormal = exp(DV);

[CAPTURE] @annotated
IPRED : Log concentration (ng/mL) without residual variability
DV    : Log concentration (ng/mL) with residual variability
IPREDnormal : concentration (ng/mL) without residual variability
DVnormal :  concentration (ng/mL) with residual variability