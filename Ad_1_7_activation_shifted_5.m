%% Reference: 
%% Sodium channels expressed in nociceptors contribute distinctly to action potential subthreshold phase, upstroke and shoulder 
%% Phil Alexander Köster, Enrico Leipoldc, Jenny Tigerholmb, Anna Maxiona, Barbara Namer, Thomas Stiehl, Angelika Lampert
%%
%% Simulation of Adelta fiber where activation of Nav1.7 is shifted by 5mV
%% This script generates Fig 12D

clear all
close all

%% setting of parameters
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

% load parameters for Nav subtypes
% these are the parameters from Table 2 of the in silico supplement
load('Nav_1_9_par.mat') % Nav1.9
load('Nav_1_8_par.mat') % Nav1.8
load('Nav_1_7_par.mat') % Nav1.7
load('Nav_1_6_par.mat') % Nav1.6
load('Nav_1_5_par.mat') % Nav1.5
load('Nav_1_3_par.mat') % Nav1.3
load('Nav_1_2_par.mat') % Nav1.2
load('Nav_1_1_par.mat') % Nav1.1

% set expression levels of Nav subtypes (taken from doi: 10.1126/scitranslmed.abj8186, rounded to two digits)
expr_1_1 = 0.09; % expression leven of Nav1.1 in Adelta
expr_1_2 = 0.09; % expression leven of Nav1.2 in Adelta
expr_1_3 = 0.01; % expression leven of Nav1.3 in Adelta
expr_1_5 = 0.03; % expression leven of Nav1.5 in Adelta
expr_1_6 = 0.40; % expression leven of Nav1.6 in Adelta
expr_1_7 = 2.02; % expression leven of Nav1.7 in Adelta
expr_1_8 = 1.18; % expression leven of Nav1.8 in Adelta
expr_1_9 = 0.92; % expression leven of Nav1.9 in Adelta

% total expression level (sum over all considered Nav subtypes)
expr_sum = expr_1_9 + expr_1_8 + expr_1_7 + expr_1_6 + expr_1_5 + expr_1_3 + expr_1_2 + expr_1_1;

% calculate relative expression for each subtype (corresponding to values in Table 4 of the in silico supplement)
f_1_9 = expr_1_9/expr_sum; % Nav1.9
f_1_8 = expr_1_8/expr_sum; % Nav1.8
f_1_7 = expr_1_7/expr_sum; % Nav1.7
f_1_6 = expr_1_6/expr_sum; % Nav1.6
f_1_5 = expr_1_5/expr_sum; % Nav1.5
f_1_3 = expr_1_3/expr_sum; % Nav1.3
f_1_2 = expr_1_2/expr_sum; % Nav1.2
f_1_1 = expr_1_1/expr_sum; % Nav1.1

% electric circuit parameters
gK_max = 72; % max. potassium conductance per area [mS/cm²]
gNa_max = 46; % max. sodium conductance per area [mS/cm²]
gL =1.0; % leak channel conductance per area [mS/cm²]
cap = 1.0; % Membrane capacitance per area [uF/cm²]
VK = 61*log(5.6/121.7)/log(10); % Potassium potential [mV], Nernst equation
VNa = 61*log(154/11.4)/log(10); % Sodium potential [mV], Nernst equation
Vl = -70.613; % Leak potential [mV]
shift_K = 65; % shift of K gating variables in relation to original Hodgkin-Huxley [mV]

% Sodium channel conductance per unit area [mS/cm^2] (corresponding to values in Table 4 of the in silico supplement)
gNa_max_1_9 = f_1_9*gNa_max; % Nav1.9
gNa_max_1_8 = f_1_8*gNa_max; % Nav1.8
gNa_max_1_7 = f_1_7*gNa_max; % Nav1.7 
gNa_max_1_6 = f_1_6*gNa_max; % Nav1.6
gNa_max_1_5 = f_1_5*gNa_max; % Nav1.5
gNa_max_1_3 = f_1_3*gNa_max; % Nav1.3
gNa_max_1_2 = f_1_2*gNa_max; % Nav1.2
gNa_max_1_1 = f_1_1*gNa_max; % Nav1.1

% load resting state of CMi fiber with increased Nav1.7 as initial condition
load("IC_Ad_shift_1_7_5.mat","IC")
shift_1_7 = 0.005; % shift of alpha_m and beta_m of Nav1.7 [V]

% check initial condition
figure()
i0=0; % injected current per area [uA/cm²] 
[T_AP,Y_AP] = ode23s(@HH_model_shifted_1_7,[0,10000],IC);
plot(T_AP, Y_AP(:,1))

%set simulation time
T_start = 0; % [ms]
T_end = 50; % [ms]

%% simulation
i0=8; % injected current per area [uA/cm²] 
[T_AP_8,Y_AP_8] = ode23s(@HH_model_shifted_1_7,[T_start,T_end],IC);
i0=10; % injected current per area [uA/cm²] 
[T_AP_10,Y_AP_10] = ode23s(@HH_model_shifted_1_7,[T_start,T_end],IC);


figure()
plot(T_AP_8,Y_AP_8(:,1),'k','LineWidth',2,'Color','#B7B7B7') % plot simulation for current injection of 8 uA/cm2
hold on
plot(T_AP_10,Y_AP_10(:,1),'k-','LineWidth',2) % plot simulation for current injection of 10 uA/cm2
yline(0,'k-','LineWidth',1.5)
xlim([0, 50])
ylim([-80,20])
xlabel('time [ms]')
ylabel ('U [mV]')
set(gca,'FontSize',20)
xticks([0,25,50])
legend('8 \muA/cm^2', '10 \muA/cm^2')
grid on
pbaspect([1 1 1])
exportgraphics(gca,'Fig12D.png','Resolution',300) 

