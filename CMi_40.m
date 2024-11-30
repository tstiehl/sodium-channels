
%% Reference: 
%% Sodium channels expressed in nociceptors contribute distinctly to action potential subthreshold phase, upstroke and shoulder 
%% Phil Alexander Köster, Enrico Leipold, Jenny Tigerholm, Anna Maxion, Barbara Namer, Thomas Stiehl, Angelika Lampert
%%
%% Simulation of CMi fiber with current injection of 40uA/cm²
%% This script generates SFig 4B, SFig 4C and upper left panel of SFig 4A

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
    gNa_max_1_9 gNa_max_1_7 gNa_max_1_8 gNa_max_1_6 gNa_max_1_5 gNa_max_1_3 gNa_max_1_2 gNa_max_1_1

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
expr_1_1 = 0.03; % expression leven of Nav1.1 in CMi
expr_1_2 = 0.07; % expression leven of Nav1.2 in CMi
expr_1_3 = 0.01; % expression leven of Nav1.3 in CMi
expr_1_5 = 0.03; % expression leven of Nav1.5 in CMi
expr_1_6 = 0.17; % expression leven of Nav1.6 in CMi
expr_1_7 = 2.31; % expression leven of Nav1.7 in CMi
expr_1_8 = 0.95; % expression leven of Nav1.8 in CMi
expr_1_9 = 2.64; % expression leven of Nav1.9 in CMi

% total expression level (sum over all considered Nav subtypes)
expr_sum = expr_1_9 + expr_1_8 + expr_1_7 + expr_1_6 + expr_1_5 + expr_1_3 + expr_1_2 + expr_1_1;

% calculate relative expression for each subtype (corresponding to values in Table 3 of the in silico supplement)
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
i0 = 0; % injected current per area [uA/cm²] 

% Sodium channel conductance per unit area [mS/cm^2] (corresponding to values in Table 3 of the in silico supplement)
gNa_max_1_9 = f_1_9*gNa_max; % Nav1.9
gNa_max_1_8 = f_1_8*gNa_max; % Nav1.8
gNa_max_1_7 = f_1_7*gNa_max; % Nav1.7
gNa_max_1_6 = f_1_6*gNa_max; % Nav1.6
gNa_max_1_5 = f_1_5*gNa_max; % Nav1.5
gNa_max_1_3 = f_1_3*gNa_max; % Nav1.3
gNa_max_1_2 = f_1_2*gNa_max; % Nav1.2
gNa_max_1_1 = f_1_1*gNa_max; % Nav1.1

% load resting state of CMi fiber as initial condition
load('IC_CMi.mat')

% check initial condition
figure()
[T_AP,Y_AP] = ode23s(@HH_model,[0,10000],IC);
plot(T_AP, Y_AP(:,1))

%% simulation
% set current injection to 40 uA/cm²
i0=40; 

% set simulation time
T_start = 0; % [ms]
T_end = 20; % [ms]

% solve HH model
[T_AP,Y_AP] = ode23s(@HH_model,[T_start,T_end],IC);

%% plotting of membrane potential
% plot AP shape
figure()
plot(T_AP,Y_AP(:,1),'k-','LineWidth',2)
hold on
yline(0,'k-','LineWidth',1.5)
xticks([0,10,20])
xlim([0, 20])
ylim([-80,20])
xlabel('time [ms]')
ylabel ('U [mV]')
set(gca,'FontSize',20)
grid on
pbaspect([1 1 1])
exportgraphics(gca,'FigS4_A_upper_left.png','Resolution',300) 

% plot zoom into AP shape
figure()
plot(T_AP,Y_AP(:,1),'k-','LineWidth',2)
hold on
yline(0,'k-','LineWidth',1.5)
xlim([0,3])
ylim([-80,20])
xlabel('time [ms]')
ylabel ('U [mV]')
set(gca,'FontSize',20)
grid on
pbaspect([1 1 1])
exportgraphics(gca,'FigS4_B_upper.png','Resolution',300) 

% plot zoom into AP shape
figure()
plot(T_AP,Y_AP(:,1),'k-','LineWidth',2)
hold on
yline(0,'k-','LineWidth',1.5)
xlim([0,1.5])
ylim([-80,20])
xlabel('time [ms]')
ylabel ('U [mV]')
set(gca,'FontSize',20)
grid on
pbaspect([1 1 1])
exportgraphics(gca,'FigS4_C_upper.png','Resolution',300) 

%% plotting of currents
% calculate currents
currents = HH_currents(Y_AP);
currents_abs = abs(currents); % take absolute values
total_currents_abs = sum(currents_abs,2); % sum over all absolute currents, including potassium
currents_abs_relative = currents_abs./total_currents_abs; % calculate relative contributions to total current

% scale action potential shape for overlay with current visualizations
AP = Y_AP(:,1);
AP_max = max(AP);
AP_min = min(AP);
AP_scaled = AP;
AP_scaled = (AP_scaled + abs(AP_min))/(AP_max-AP_min);

% plot relative contribution of all currents to total current
figure()
a=area(T_AP,currents_abs_relative,'LineStyle','none');
% assign colors (colors should match other figures)
a(1).FaceColor="#FF007F"; % K
a(2).FaceColor="#F0E442"; % Nav1.9
a(3).FaceColor="#D55E00"; % Nav1.8
a(4).FaceColor="#6A6A6A"; % Nav1.7
a(5).FaceColor="#CC79A7"; % Nav1.6
a(6).FaceColor="#0072B2"; % Nav1.5
a(7).FaceColor="#009E73"; % Nav1.3
a(8).FaceColor="#E69F00"; % Nav1.2
a(9).FaceColor="#56B4E9"; % Nav1.1
a(10).FaceColor="#AEF359"; %L
hold on
plot(T_AP, AP_scaled, '-k','LineWidth',2) % overlay action potential (scaled to match axis range)
xlim([0,3])
xline(1,'--k')
xline(2,'--k')
ylim([0,1])
xlabel('time [ms]')
ylabel ({'contribution to current'})
set(gca,'FontSize',20)
xticks([0,1,2,3])
grid on
pbaspect([1 1 1])
exportgraphics(gca,'FigS4_B_middle.png','Resolution',300) 

% plot zoom into previous plot 
figure()
a=area(T_AP,currents_abs_relative,'LineStyle','none');
% assign colors (colors should match other figures)
a(1).FaceColor="#FF007F"; % K
a(2).FaceColor="#F0E442"; % Nav1.9
a(3).FaceColor="#D55E00"; % Nav1.8
a(4).FaceColor="#6A6A6A"; % Nav1.7
a(5).FaceColor="#CC79A7"; % Nav1.6
a(6).FaceColor="#0072B2"; % Nav1.5
a(7).FaceColor="#009E73"; % Nav1.3
a(8).FaceColor="#E69F00"; % Nav1.2
a(9).FaceColor="#56B4E9"; % Nav1.1
a(10).FaceColor="#AEF359"; %L
hold on
plot(T_AP, AP_scaled, '-k','LineWidth',2) % overlay action potential (scaled to match axis range)
xlim([0,1.5])
xline(1,'--k')
xline(0.5,'--k')
ylim([0,1])
xlabel('time [ms]')
ylabel ({'contribution to current'})
set(gca,'FontSize',20)
xticks([0,0.5,1.0,1.5])
grid on
pbaspect([1 1 1])
exportgraphics(gca,'FigS4_C_middle.png','Resolution',300) 

% repeat plotting taking out potassium and leak current
sodium_currents_abs = abs(currents(:, 2:9)); % select sodium currents
total_sodium_currents_abs = sum(sodium_currents_abs,2);
sodium_currents_abs_relative = sodium_currents_abs./total_sodium_currents_abs; % calculate relative contributions

figure()
a=area(T_AP,sodium_currents_abs_relative,'LineStyle','none');
a(1).FaceColor="#F0E442"; % Nav1.9
a(2).FaceColor="#D55E00"; % Nav1.8
a(3).FaceColor="#6A6A6A"; % Nav1.7
a(4).FaceColor="#CC79A7"; % Nav1.6
a(5).FaceColor="#0072B2"; % Nav1.5
a(6).FaceColor="#009E73"; % Nav1.3
a(7).FaceColor="#E69F00"; % Nav1.2
a(8).FaceColor="#56B4E9"; % Nav1.1
xlim([0,3])
ylim([0.2,1])
xlabel('time [ms]')
ylabel ({'contribution to Na^+ current'})
set(gca,'FontSize',19)
xline(1,'--k')
xline(2,'--k')
hold on
plot(T_AP, AP_scaled*0.8+0.2, '-k','LineWidth',2) % overlay action potential shape, rescaled to fit y-axis range
grid on
pbaspect([1 1 1])
exportgraphics(gca,'FigS4_B_lower.png','Resolution',300) 

% zoom into previous plot
figure()
a=area(T_AP,sodium_currents_abs_relative,'LineStyle','none');
a(1).FaceColor="#F0E442"; % Nav1.9
a(2).FaceColor="#D55E00"; % Nav1.8
a(3).FaceColor="#6A6A6A"; % Nav1.7
a(4).FaceColor="#CC79A7"; % Nav1.6
a(5).FaceColor="#0072B2"; % Nav1.5
a(6).FaceColor="#009E73"; % Nav1.3
a(7).FaceColor="#E69F00"; % Nav1.2
a(8).FaceColor="#56B4E9"; % Nav1.1
xlim([0,1.5])
ylim([0.2,1])
xlabel('time [ms]')
ylabel ({'contribution to Na^+ current'})
set(gca,'FontSize',20)
xline(1,'--k')
xline(0.5,'--k')
hold on
plot(T_AP, AP_scaled*0.8+0.2, '-k','LineWidth',2) % overlay action potential shape, rescaled to fit y-axis range
grid on
pbaspect([1 1 1])
exportgraphics(gca,'FigS4_C_lower.png','Resolution',300) 

