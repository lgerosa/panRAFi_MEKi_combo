# exported from PySB model 'RAS_RAF_RAFi_rules_compartments_energy'

from pysb import Model, Monomer, Parameter, Expression, Compartment, Rule, Observable, Initial, MatchOnce, EnergyPattern, Annotation, MultiState, Tag, ANY, WILD, as_complex_pattern
from sympy import log

Model()

Monomer('Ras', ['raf', 'state'], {'state': ['GTP']})
Monomer('Raf', ['ras', 'raf', 'i'])
Monomer('I', ['raf'])

Parameter('M_vol', 0.001)
Parameter('C_vol', 1.0)
Parameter('RasGTP_0', 0.0)
Parameter('Raf_0', 0.01)
Parameter('I_0', 0.0)
Parameter('kr_GR', 10.0)
Parameter('kf_GR', 1.0)
Parameter('kr_RR', 10.0)
Parameter('kf_RR', 1.0)
Parameter('phi_RR', 1.0)
Parameter('kr_RI', 0.1)
Parameter('kf_RI', 1.0)
Parameter('phi_RI', 1.0)
Parameter('f', 1.0)
Parameter('g', 1.0)

Expression('Gf_RR', log(kr_RR/kf_RR))
Expression('Ea0_RR', -phi_RR*log(kr_RR/kf_RR) - log(kf_RR))
Expression('Gf_RI', log(kr_RI/kf_RI))
Expression('Ea0_RI', -phi_RI*log(kr_RI/kf_RI) - log(kf_RI))
Expression('Gf_f', log(f))
Expression('Gf_g', log(g))
Expression('Gf_RRI', Gf_f)
Expression('Gf_IRRI', Gf_f + Gf_g)

Compartment(name='M', parent=None, dimension=2, size=M_vol)
Compartment(name='C', parent=M, dimension=3, size=C_vol)

Observable('Rtot_obs', Raf())
Observable('Itot_obs', I())
Observable('R_BRAFmut_active_obs', Raf(i=None))
Observable('R_RAFwt_active_obs', Raf(raf=1, i=None) % Raf(raf=1))
Observable('R_obs', Raf(raf=None, i=None), match='species')
Observable('RR_obs', Raf(raf=1, i=None) % Raf(raf=1, i=None), match='species')
Observable('RRI_obs', Raf(raf=1, i=None) % Raf(raf=1, i=2) % I(raf=2), match='species')
Observable('IRRI_obs', I(raf=2) % Raf(raf=1, i=2) % Raf(raf=1, i=3) % I(raf=3), match='species')

Rule('GR', Ras(raf=None, state='GTP') + Raf(ras=None) | Ras(raf=1, state='GTP') % Raf(ras=1), kf_GR, kr_GR)
Rule('RR', Raf(raf=None) + Raf(raf=None) | Raf(raf=1) % Raf(raf=1), phi_RR, Ea0_RR, energy=True)
Rule('RI', Raf(i=None) + I(raf=None) | Raf(i=1) % I(raf=1), phi_RI, Ea0_RI, energy=True)

EnergyPattern('ep_RR', Raf(raf=1) % Raf(raf=1), Gf_RR)
EnergyPattern('ep_RI', Raf(i=1) % I(raf=1), Gf_RI)
EnergyPattern('ep_RRI', Raf(raf=1, i=None) % Raf(raf=1, i=2) % I(raf=2), Gf_RRI)
EnergyPattern('ep_IRRI', I(raf=2) % Raf(raf=1, i=2) % Raf(raf=1, i=3) % I(raf=3), Gf_IRRI)

Initial(Ras(raf=None, state='GTP') ** M, RasGTP_0)
Initial(Raf(ras=None, raf=None, i=None) ** C, Raf_0)
Initial(I(raf=None) ** C, I_0)

