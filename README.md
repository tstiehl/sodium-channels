# sodium-channel
Modified Hodgkin-Huxley model of CMi and Adelta fiber accounting for the sodium channels Nav1.1-Nav1.3 and Nav1.5-Nav1.9.

This is the Matlab code to: Sodium channels expressed in nociceptors contribute distinctly to action potential subthreshold phase, upstroke and shoulder 
by Phil Alexander Köster, Enrico Leipold, Jenny Tigerholm, Anna Maxion, Barbara Namer, Thomas Stiehl, Angelika Lampert.

## EXECUTABLE SCRIPTS
### Scripts for CMi fiber
- CMi_30.m: Simulation of CMi fiber with a current injection of 30 uA/cm². This script generates Fig. 9B, Fig. 9C, upper left panel of Fig. 9A, SFig. 8A and SFig. 8B.
- CMi_30_subtype_contribution.m: Simulation of CMi fiber with current injection of 30uA/cm² and with individual Nav subtypes removed during action potential. This script generates Fig. 9A and SFig. 6A.
- CMi_40.m: Simulation of CMi fiber with a current injection of 40 uA/cm². This script generates SFig. 4B, SFig. 4C and upper left panel of SFig. 4A.
- CMi_40_subtype_contribution.m: Simulation of CMi fiber with current injection of 40uA/cm² and with individual Nav subtypes removed during action potential. This script generates SFig 5A and SFig 6D.
- CMi_1_9.m: Simulation of CMi fiber with different expression levels ov Nav1.9 (current injection 30uA/cm²). This script generates SFig. 7.
- CMi_1_7_increased.m: Simulation of CMi fiber with 5-fold increased expression of Nav1.7. This script generates Fig. 11A.
- CMi_1_8_increased.m: Simulation of CMi fiber with 11-fold increased expression of Nav1.8. This script generates Fig. 11B.
- CMi_1_7_activation_shifted_0.m: Simulation of CMi fiber with different current injections. This script generates Fig. 11C.
- CMi_1_7_activation_shifted_5.m: Simulation of CMi fiber where Nav1.7 activation is shifted by 5mV. This script generates Fig. 11D.
- CMi_1_7_activation_shifted_8.m: Simulation of CMi fiber where Nav1.7 activation is shifted by 8mV. This script generates Fig 11E.
### Scripts for Adelta filber
- Ad_30.m: Simulation of Adelta fiber with a current injection of 30 uA/cm². This script generates Fig. 10B, Fig. 10C, upper left panel of Fig. 10A, SFig. 8C and SFig. 8D.
- Ad_30_subtype_contribution.m: Simulation of Adelta fiber with current injection of 30uA/cm² and with individual Nav subtypes removed during action potential. 
- Ad_40.m: Simulation of Adelta fiber with a current injection of 40 uA/cm². This script generates SFig. 5B, SFig. 5C, upper left panel of SFig. 5A.
- Ad_40_subtype_contribution.m: Simulation of Adelta fiber with current injection of 40uA/cm² and with individual Nav subtypes removed during action potential. 
- Ad_1_7_increased.m: Simulation of Adelta fiber with 5-fold increased expression of Nav1.7. This script generates Fig. 12A.
- Ad_1_8_increased.m: Simulation of Adelta fiber with 8-fold increased expression of Nav1.8. This script generates Fig. 12B.
- Ad_1_7_activation_shifted_0.m: Simulation of Adelta fiber with different current injections. This script generates Fig. 12C.
- Ad_1_7_activation_shifted_5.m: Simulation of Adelta fiber where Nav1.7 activation is shifted by 5mV. This script generates Fig. 12D.
- Ad_1_7_activation_shifted_8.m: Simulation of Adelta fiber where Nav1.7 activation is shifted by 8mV. This script generates Fig. 12E.


## FILES CONTAINING FUNCTIONS
- HH_model.m: Implementation of Hodgkin-Huxley ODE model with Nav1.1-Nav1.3, Nav1.4-Nav1.9
- HH_currents.m: Function to calculate currents based on solution obtained from the HH_model.m
- HH_conductances.m: Function to calculate conductances based on solution obtained from the HH_model.m
- HH_model_shifted_1_7.m: Implementation of Hodgkin-Huxley ODE model with Nav1.1-Nav1.3, Nav1.4-Nav1.9. The terms alpha_m and beta_m of Nav1.7 are shifted to model pathologies.
- chi.m: Implementation of function chi from equation (1) of the in silico supplement (parameterization of voltage dependent gating variables)
- alpha_n.m: Parametrization of voltage-dependency of gating variable alpha_n for generic potassium channel (equation (13) of the in silico supplement)
- beta_n.m: Parametrization of voltage-dependency of gating variable beta_n for generic potassium channel (equation (14) of the in silico supplement)


## PARAMETER FILES
- Nav_1_1_par.mat: contains gating parameters for Nav1.1 (as tabulated in section 4.1 of the in silico supplement)
- Nav_1_2_par.mat: contains gating parameters for Nav1.2 (as tabulated in section 4.2 of the in silico supplement)
- Nav_1_3_par.mat: contains gating parameters for Nav1.3 (as tabulated in section 4.3 of the in silico supplement)
- Nav_1_5_par.mat: contains gating parameters for Nav1.5 (as tabulated in section 4.4 of the in silico supplement)
- Nav_1_6_par.mat: contains gating parameters for Nav1.6 (as tabulated in section 4.5 of the in silico supplement)
- Nav_1_7_par.mat: contains gating parameters for Nav1.7 (as tabulated in section 4.6 of the in silico supplement)
- Nav_1_8_par.mat: contains gating parameters for Nav1.8 (as tabulated in section 4.7 of the in silico supplement)
- Nav_1_9_par.mat: contains gating parameters for Nav1.9 (as tabulated in section 4.8 of the in silico supplement)
- IC_CMi.mat: resting state of CMi fiber
- IC_CMi_1_9_100.mat: resting state of CMi fiber without Nav1.9 expression
- IC_CMi_1_9_75.mat: resting state of CMi fiber with Nav1.9 expression reduced by 75%
- IC_CMi_1_9_50.mat: : resting state of CMi fiber with Nav1.9 expression reduced by 50%
- IC_CMi_1_7.mat: resting state of CMi fiber with 5-fold increased Nav1.7 expression 
- IC_CMi_1_8.mat: resting state of CMi fiber with 11-fold increased Nav1.8 expression 
- IC_CMi_shift_1_7_5.mat: resting state of CMi fiber, where Nav1.7 activation is shifted by 5mV
- IC_CMi_shift_1_7_8.mat: resting state of CMi fiber, where Nav1.7 activation is shifted by 8mV
- IC_Ad.mat: resting state of Adelta fiber
- IC_Ad_1_7.mat: resting state of Adelta fiber with 5-fold increased Nav1.7 expression 
- IC_Ad_1_8.mat: resting state of Adelta fiber with 8-fold increased Nav1.8 expression 
- IC_Ad_shift_1_7_5.mat: resting state of Adelta fiber, where Nav1.7 activation is shifted by 5mV
- IC_Ad_shift_1_7_8.mat: resting state of Adelta fiber, where Nav1.7 activation is shifted by 8mV


## OUTPUT
- output.zip: zip file containing the output generated by the above scripts. The output was generated by running the code on MATLAB version '24.1.0.2537033 (R2024a)' on Ubuntu 24.04.1 LTS using 11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz.
