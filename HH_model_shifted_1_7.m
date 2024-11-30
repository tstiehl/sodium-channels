%% Reference: 
%% Sodium channels expressed in nociceptors contribute distinctly to action potential subthreshold phase, upstroke and shoulder 
%% Phil Alexander Köster, Enrico Leipold, Jenny Tigerholm, Anna Maxion, Barbara Namer, Thomas Stiehl, Angelika Lampert
%%
%% Implementation of the modified Hodgkin-Huxley ODE model
%% alpha_m and beta_m of Nav1.7 are shifted by shift1_7 to model pathologies
%% Equation numbers refer to the in silico supplement

function Y= HH_model_shifted_1_7(t,y)

global gK_max gL cap VK VNa Vl shift_K i0 ...
    alpha_m_1_8_par beta_m_1_8_par alpha_h_1_8_par beta_h_1_8_par delay_h_1_8_par ...
    alpha_m_1_7_par beta_m_1_7_par alpha_h_1_7_par beta_h_1_7_par delay_h_1_7_par ...
    alpha_m_1_6_par beta_m_1_6_par alpha_h_1_6_par beta_h_1_6_par delay_h_1_6_par ...
    alpha_m_1_5_par beta_m_1_5_par alpha_h_1_5_par beta_h_1_5_par delay_h_1_5_par ...
    alpha_m_1_3_par beta_m_1_3_par alpha_h_1_3_par beta_h_1_3_par delay_h_1_3_par ...
    alpha_m_1_2_par beta_m_1_2_par alpha_h_1_2_par beta_h_1_2_par delay_h_1_2_par ...
    alpha_m_1_1_par beta_m_1_1_par alpha_h_1_1_par beta_h_1_1_par delay_h_1_1_par ...
    alpha_m_1_9_par beta_m_1_9_par alpha_h_1_9_par beta_h_1_9_par  ...
    gNa_max_1_9 gNa_max_1_7 gNa_max_1_8 gNa_max_1_6 gNa_max_1_5 gNa_max_1_3 gNa_max_1_2 gNa_max_1_1 shift_1_7


%% model variables
V = y(1); % membrane potential [mV], defined in equation (3)
m_1_8 = y(2); % gating variable m for Nav1.8, defined in equation (7) 
n = y(3); % gating variable n for generic potassium channel, defined in eqaution (4)

h1_1_8 = y(4); % gating variable h1 for Nav1.8, defined in equation (8)
w_1_8 = y(5:9); % gating variabled hw1, ..., hw5 for Nav1.8, defined in equations (9)-(10)
h0_1_8 = y(10); % gating variable h0 for Nav1.8, defined in equation (11).

m_1_7 = y(11); % gating variable m for Nav1.7, defined in equation (7)
h1_1_7 = y(12); % gating variable h1 for Nav1.7, defined in equation (8)
w_1_7 = y(13:17); % gating variabled hw1, ..., hw5 for Nav1.7, defined in equations (9)-(10)
h0_1_7 = y(18); % gating variable h0 for Nav1.7, defined in equation (11).

m_1_6 = y(19); % gating variable m for Nav1.6, defined in equation (7)
h1_1_6 = y(20); % gating variable h1 for Nav1.6, defined in equation (8)
w_1_6 = y(21:25); % gating variabled hw1, ..., hw5 for Nav1.6, defined in equations (9)-(10)
h0_1_6 = y(26); % gating variable h0 for Nav1.6, defined in equation (11).

m_1_3 = y(27); % gating variable m for Nav1.3, defined in equation (7)
h1_1_3 = y(28); % gating variable h1 for Nav1.3, defined in equation (8)
w_1_3 = y(29:33); % gating variabled hw1, ..., hw5 for Nav1.3, defined in equations (9)-(10)
h0_1_3 = y(34); % gating variable h0 for Nav1.3, defined in equation (11).

m_1_2 = y(35); % gating variable m for Nav1.2, defined in equation (7)
h1_1_2 = y(36); % gating variable h1 for Nav1.2, defined in equation (8)
w_1_2 = y(37:41); % gating variabled hw1, ..., hw5 for Nav1.2, defined in equations (9)-(10)
h0_1_2 = y(42); % gating variable h0 for Nav1.2, defined in equation (11).

m_1_1 = y(43); % gating variable m for Nav1.1, defined in equation (7)
h1_1_1 = y(44); % gating variable h1 for Nav1.1, defined in equation (8)
w_1_1 = y(45:49); % gating variabled hw1, ..., hw5 for Nav1.1, defined in equations (9)-(10)
h0_1_1 = y(50); % gating variable h0 for Nav1.1 defined in equation (11).

m_1_5 = y(51); % gating variable m for Nav1.5, defined in equation (7)
h1_1_5 = y(52); % gating variable h1 for Nav1.5, defined in equation (8)
w_1_5 = y(53:57); % gating variabled hw1, ..., hw5 for Nav1.5, defined in equations (9)-(10)
h0_1_5 = y(58); % gating variable h0 for Nav1.5 defined in equation (11).

m_1_9 = y(59); % gating variable m for Nav1.9, defined in equation (5)
h_1_9 = y(60); % gating variable h for Nav1.9, defined in equation (6)


% functions to calculate values of gating variables in dependence of
% membrane potential: x is membrane potential in VOLTS

% Nav1.9
% The data used for parameter fitting of Nav1.9 were not corrected for liquid
% junction potential, therefore, a shift of 0.007 V is added to x
% data used for parameter fitting of all other Nav subtypes were corrected
% for liquid junction potential, therefore, no shift is added in case of
% other Nav subtypes
alpha_m_1_9 = @(x)max(0, chi(alpha_m_1_9_par,x+0.007)); % alpha_m for Nav 1.9, used in equation (5)
beta_m_1_9 = @(x)max(0, chi(beta_m_1_9_par,x+0.007)); % beta_m for Nav 1.9, used in equation (5)
alpha_h_1_9 = @(x)max(0.00022, chi(alpha_h_1_9_par,x+0.007)); % alpha_h for Nav 1.9, used in equation (6), the cutoff 0.00022 avoids alpha_h_1_9 to be zero within  the physiological range of membrane potentials. 
beta_h_1_9 = @(x)max(0.01, chi(beta_h_1_9_par,x+0.007)); % beta_h for Nav 1.9, used in equation (6), the cutoff 0.01 avoids alpha_h_1_9 to be zero within  the physiological range of membrane potentials. 

% Nav1.8
alpha_m_1_8 = @(x)max(0, chi(alpha_m_1_8_par,x)); % alpha_m for Nav 1.8, used in equation (7)
beta_m_1_8 = @(x)max(0, chi(beta_m_1_8_par,x));  % beta_m for Nav 1.8, equation (7)
alpha_h_1_8 = @(x)max(0, chi(alpha_h_1_8_par,x));  % alpha_h for Nav 1.8, used in equation (11)
beta_h_1_8 = @(x)max(0, chi(beta_h_1_8_par,x));  % beta_h for Nav 1.8, used in equation (8)
delay_h_1_8 = @(x)max(0, chi(delay_h_1_8_par,x)); % gamma_h for Nav1.8, used in euqations (9-11)

% Nav1.7
alpha_m_1_7 = @(x)max(0, chi(alpha_m_1_7_par,x+shift_1_7)); % alpha_m for Nav 1.7, used in equation (7), here shifted by shift1_7
beta_m_1_7 = @(x)max(0, chi(beta_m_1_7_par,x+shift_1_7));  % beta_m for Nav 1.7, used in equation (7), here shifted by shift1_7
alpha_h_1_7 = @(x)max(0, chi(alpha_h_1_7_par,x));  % alpha_h for Nav 1.7, used in equation (11)
beta_h_1_7 = @(x)max(0, chi(beta_h_1_7_par,x));  % beta_h for Nav 1.7, used in equation (8)
delay_h_1_7 = @(x)max(0, chi(delay_h_1_7_par,x)); % gamma_h for Nav1.7, used in euqations (9-11)

% Nav1.6
alpha_m_1_6 = @(x)max(0, chi(alpha_m_1_6_par,x)); % alpha_m for Nav 1.6, used in equation (7)
beta_m_1_6 = @(x)max(0, chi(beta_m_1_6_par,x));  % beta_m for Nav 1.6, used in equation (7)
alpha_h_1_6 = @(x)max(0, chi(alpha_h_1_6_par,x));  % alpha_h for Nav 1.6, used in equation (11)
beta_h_1_6 = @(x)max(0, chi(beta_h_1_6_par,x));  % beta_h for Nav 1.6, used in equation (8)
delay_h_1_6 = @(x)max(0, chi(delay_h_1_6_par,x)); % gamma_h for Nav1.6, euqations (9-11)

% Nav1.5
alpha_m_1_5 = @(x)max(0, chi(alpha_m_1_5_par,x)); % alpha_m for Nav 1.5, used in equation (7)
beta_m_1_5 = @(x)max(0, chi(beta_m_1_5_par,x));  % beta_m for Nav 1.5, used in equation (7)
alpha_h_1_5 = @(x)max(0, chi(alpha_h_1_5_par,x));  % alpha_h for Nav 1.5, used in equation (11)
beta_h_1_5 = @(x)max(0, chi(beta_h_1_5_par,x));  % beta_h for Nav 1.5, used in equation (8)
delay_h_1_5 = @(x)max(0, chi(delay_h_1_5_par,x));

% Nav1.3
alpha_m_1_3 = @(x)max(0, chi(alpha_m_1_3_par,x)); % alpha_m for Nav 1.3, used in equation (7)
beta_m_1_3 = @(x)max(0, chi(beta_m_1_3_par,x));  % beta_m for Nav 1.3, used in equation (7)
alpha_h_1_3 = @(x)max(0, chi(alpha_h_1_3_par,x));  % alpha_h for Nav 1.3, used in equation (11)
beta_h_1_3 = @(x)max(0, chi(beta_h_1_3_par,x));  % beta_h for Nav 1.3, used in equation (8)
delay_h_1_3 = @(x)max(0, chi(delay_h_1_3_par,x)); % gamma_h for Nav1.3, euqations (9-11)

% Nav1.2
alpha_m_1_2 = @(x)max(0, chi(alpha_m_1_2_par,x)); % alpha_m for Nav 1.2, used in equation (7)
beta_m_1_2 = @(x)max(0, chi(beta_m_1_2_par,x));  % beta_m for Nav 1.2, used in equation (7)
alpha_h_1_2 = @(x)max(0, chi(alpha_h_1_2_par,x));  % alpha_h for Nav 1.2, used in equation (11)
beta_h_1_2 = @(x)max(0, chi(beta_h_1_2_par,x));  % beta_h for Nav 1.2, used in equation (8)
delay_h_1_2 = @(x)max(0, chi(delay_h_1_2_par,x)); % gamma_h for Nav1.2, euqations (9-11)

% Nav1.1
alpha_m_1_1 = @(x)max(0, chi(alpha_m_1_1_par,x));  % alpha_m for Nav 1.1, used in equation (7)
beta_m_1_1 = @(x)max(0, chi(beta_m_1_1_par,x));  % beta_m for Nav 1.1, used in equation (7)
alpha_h_1_1 = @(x)max(0, chi(alpha_h_1_1_par,x));  % alpha_h for Nav 1.1, used in equation (11)
beta_h_1_1 = @(x)max(0, chi(beta_h_1_1_par,x));  % beta_h for Nav 1.1, used in equation (8)
delay_h_1_1 = @(x)max(0, chi(delay_h_1_1_par,x)); % gamma_h for Nav1.1, euqations (9-11)


%% ODEs for inactivating gates
% Nav1.8 (modified Hodgkin-Huxley model)
dh1_1_8 =  - beta_h_1_8(V/1000) * h1_1_8 + alpha_h_1_8(V/1000) * h0_1_8; % state h1, equation (8)
dw1_1_8 =  beta_h_1_8(V/1000) * h1_1_8 - delay_h_1_8(V/1000)*w_1_8(1); % state h_w1, equation (9)
dw_2_to_5_1_8 = delay_h_1_8(V/1000)*w_1_8(1:4) - delay_h_1_8(V/1000)*w_1_8(2:5); % states h_w2 to h_w5, equation (10)
dh0_1_8= -alpha_h_1_8(V/1000) * h0_1_8 + delay_h_1_8(V/1000)*w_1_8(5); % state h0, equation (11)

% Nav1.7 (modified Hodgkin-Huxley model)
dh1_1_7 =  - beta_h_1_7(V/1000) * h1_1_7 + alpha_h_1_7(V/1000) * h0_1_7; % state h1, equation (8)
dw1_1_7 =  beta_h_1_7(V/1000) * h1_1_7 - delay_h_1_7(V/1000)*w_1_7(1); % state h_w1, equation (9)
dw_2_to_5_1_7 = delay_h_1_7(V/1000)*w_1_7(1:4) - delay_h_1_7(V/1000)*w_1_7(2:5); % states h_w2 to h_w5, equation (10)
dh0_1_7= -alpha_h_1_7(V/1000) * h0_1_7 + delay_h_1_7(V/1000)*w_1_7(5); % state h0, equation (11)

% Nav1.5 (modified Hodgkin-Huxley model)
dh1_1_5 =  - beta_h_1_5(V/1000) * h1_1_5 + alpha_h_1_5(V/1000) * h0_1_5; % state h1, equation (8)
dw1_1_5 =  beta_h_1_5(V/1000) * h1_1_5 - delay_h_1_5(V/1000)*w_1_5(1); % state h_w1, equation (9)
dw_2_to_5_1_5 = delay_h_1_5(V/1000)*w_1_5(1:4) - delay_h_1_5(V/1000)*w_1_5(2:5); % states h_w2 to h_w5, equation (10)
dh0_1_5= -alpha_h_1_5(V/1000) * h0_1_5 + delay_h_1_5(V/1000)*w_1_5(5); % state h0, equation (11)

% Nav1.6 (modified Hodgkin-Huxley model)
dh1_1_6 =  - beta_h_1_6(V/1000) * h1_1_6 + alpha_h_1_6(V/1000) * h0_1_6; % state h1, equation (8)
dw1_1_6 =  beta_h_1_6(V/1000) * h1_1_6 - delay_h_1_6(V/1000)*w_1_6(1); % state h_w1, equation (9)
dw_2_to_5_1_6 = delay_h_1_6(V/1000)*w_1_6(1:4) - delay_h_1_6(V/1000)*w_1_6(2:5); % states h_w2 to h_w5, equation (10)
dh0_1_6= -alpha_h_1_6(V/1000) * h0_1_6 + delay_h_1_6(V/1000)*w_1_6(5); % state h0, equation (11)

% Nav1.3 (modified Hodgkin-Huxley model)
dh1_1_3 =  - beta_h_1_3(V/1000) * h1_1_3 + alpha_h_1_3(V/1000) * h0_1_3; % state h1, equation (8)
dw1_1_3 =  beta_h_1_3(V/1000) * h1_1_3 - delay_h_1_3(V/1000)*w_1_3(1); % state h_w1, equation (9)
dw_2_to_5_1_3 = delay_h_1_3(V/1000)*w_1_3(1:4) - delay_h_1_3(V/1000)*w_1_3(2:5); % states h_w2 to h_w5, equation (10)
dh0_1_3= -alpha_h_1_3(V/1000) * h0_1_3 + delay_h_1_3(V/1000)*w_1_3(5); % state h0, equation (11)

% Nav1.2 (modified Hodgkin-Huxley model)
dh1_1_2 =  - beta_h_1_2(V/1000) * h1_1_2 + alpha_h_1_2(V/1000) * h0_1_2; % state h1, equation (8)
dw1_1_2 =  beta_h_1_2(V/1000) * h1_1_2 - delay_h_1_2(V/1000)*w_1_2(1); % state h_w1, equation (9)
dw_2_to_5_1_2 = delay_h_1_2(V/1000)*w_1_2(1:4) - delay_h_1_2(V/1000)*w_1_2(2:5); % states h_w2 to h_w5, equation (10)
dh0_1_2= -alpha_h_1_2(V/1000) * h0_1_2 + delay_h_1_2(V/1000)*w_1_2(5); % state h0, equation (11)

% Nav1.1 (modified Hodgkin-Huxley model)
dh1_1_1 =  - beta_h_1_1(V/1000) * h1_1_1 + alpha_h_1_1(V/1000) * h0_1_1; % state h1, equation (8)
dw1_1_1 =  beta_h_1_1(V/1000) * h1_1_1 - delay_h_1_1(V/1000)*w_1_1(1); % state h_w1, equation (9)
dw_2_to_5_1_1 = delay_h_1_1(V/1000)*w_1_1(1:4) - delay_h_1_1(V/1000)*w_1_1(2:5); % states h_w2 to h_w5, equation (10)
dh0_1_1= -alpha_h_1_1(V/1000) * h0_1_1 + delay_h_1_1(V/1000)*w_1_1(5); % state h0, equation (11)

% Nav1.9 (ORIGINAL Hodgkin-Huxley model)
dh_1_9 = (alpha_h_1_9(V/1000) * (1.0 - h_1_9)) - (beta_h_1_9(V/1000) * h_1_9); % state h, equation (6)


%% ODEs for activating gates
dm_1_9 = (alpha_m_1_9(V/1000) * (1.0 - m_1_9)) - (beta_m_1_9(V/1000) * m_1_9); % state m for Nav1.9, equation (5)
dm_1_8 = (alpha_m_1_8(V/1000) * (1.0 - m_1_8)) - (beta_m_1_8(V/1000) * m_1_8); % state m for Nav1.8, equation (7)  
dm_1_7 = (alpha_m_1_7(V/1000) * (1.0 - m_1_7)) - (beta_m_1_7(V/1000) * m_1_7); % state m for Nav1.7, equation (7)
dm_1_6 = (alpha_m_1_6(V/1000) * (1.0 - m_1_6)) - (beta_m_1_6(V/1000) * m_1_6); % state m for Nav1.6, equation (7)
dm_1_5 = (alpha_m_1_5(V/1000) * (1.0 - m_1_5)) - (beta_m_1_5(V/1000) * m_1_5); % state m for Nav1.5, equation (7)
dm_1_3 = (alpha_m_1_3(V/1000) * (1.0 - m_1_3)) - (beta_m_1_3(V/1000) * m_1_3); % state m for Nav1.3, equation (7)
dm_1_2 = (alpha_m_1_2(V/1000) * (1.0 - m_1_2)) - (beta_m_1_2(V/1000) * m_1_2); % state m for Nav1.2, equation (7)
dm_1_1 = (alpha_m_1_1(V/1000) * (1.0 - m_1_1)) - (beta_m_1_1(V/1000) * m_1_1); % state m for Nav1.1, equation (7)
dn = (alpha_n(V,shift_K) * (1.0 - n)) - (beta_n(V,shift_K) * n); % state m for generic potassium channel, equation (4)
% gating variables for sodium channels are parameterized in volts, gating variables for potassium channel are parameterized in millivolts 

%% conductance-capacitance ratios for current time point
% conductances are calculated as specified in equation (3)
K = (gK_max / cap) * n^4; % generic potassium channel
Na_1_9 = (gNa_max_1_9  / cap) * m_1_9^3 * h_1_9; % Nav1.9 (ORIGINAL Hodgkin-Huxley)
Na_1_8 = (gNa_max_1_8  / cap) * m_1_8^3 * (1-h0_1_8); % Nav1.8 (modified Hodgkin-Huxley)
Na_1_7 = (gNa_max_1_7  / cap) * m_1_7^3 * (1-h0_1_7); % Nav1.7 (modified Hodgkin-Huxley)
Na_1_6 = (gNa_max_1_6  / cap) * m_1_6^3 * (1-h0_1_6); % Nav1.6 (modified Hodgkin-Huxley)
Na_1_5 = (gNa_max_1_5  / cap) * m_1_5^3 * (1-h0_1_5); % Nav1.5 (modified Hodgkin-Huxley)
Na_1_3 = (gNa_max_1_3  / cap) * m_1_3^3 * (1-h0_1_3); % Nav1.3 (modified Hodgkin-Huxley)
Na_1_2 = (gNa_max_1_2  / cap) * m_1_2^3 * (1-h0_1_2); % Nav1.2 (modified Hodgkin-Huxley)
Na_1_1 = (gNa_max_1_1  / cap) * m_1_1^3 * (1-h0_1_1); % Nav1.1 (modified Hodgkin-Huxley)
L = gL / cap; % ratio of leakeage conductance to membrane capacitance

I = i0; % injected current per area [uA/cm²]

%% ODE for membrane potential
% equation (3)
dV = (I / cap) - (K * (V - VK)) ...
    - (Na_1_9 * (V - VNa))...
    - (Na_1_8 * (V - VNa))...
    - (Na_1_7 * (V - VNa)) ...
    - (Na_1_6 * (V - VNa)) ...
    - (Na_1_5 * (V - VNa)) ...
    - (Na_1_3 * (V - VNa)) ...
    - (Na_1_2 * (V - VNa)) ...
    - (Na_1_1 * (V - VNa)) ...
    - (L * (V - Vl));


% formatting of return value for ODE solver
Y=[dV;dm_1_8;dn;dh1_1_8;dw1_1_8; dw_2_to_5_1_8; dh0_1_8; dm_1_7;dh1_1_7;dw1_1_7; dw_2_to_5_1_7;dh0_1_7;dm_1_6;dh1_1_6;dw1_1_6; dw_2_to_5_1_6;dh0_1_6;dm_1_3;dh1_1_3;dw1_1_3; dw_2_to_5_1_3;dh0_1_3;dm_1_2;dh1_1_2;dw1_1_2; dw_2_to_5_1_2;dh0_1_2;dm_1_1;dh1_1_1;dw1_1_1; dw_2_to_5_1_1;dh0_1_1;dm_1_5;dh1_1_5;dw1_1_5; dw_2_to_5_1_5;dh0_1_5;dm_1_9;dh_1_9];

end

