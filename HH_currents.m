% calculate currents based on the solution of HH_model.m

function Y= HH_currents(y)
% y: solution from ODE solver

global gK_max gL VK VNa Vl gNa_max_1_9 gNa_max_1_7 gNa_max_1_8 gNa_max_1_6 gNa_max_1_5 gNa_max_1_3 gNa_max_1_2 gNa_max_1_1


%% model variables
V = y(:,1); % membrane potential [mV], defined in equation (3)
m_1_8 = y(:,2); % gating variable m for Nav1.8, defined in equation (7) 
n = y(:,3); % gating variable n for generic potassium channel, defined in eqaution (4)
h0_1_8 = y(:,10); % gating variable h0 for Nav1.8, defined in equation (11).
m_1_7 = y(:,11); % gating variable m for Nav1.7, defined in equation (7)
h0_1_7 = y(:,18); % gating variable h0 for Nav1.7, defined in equation (11).
m_1_6 = y(:,19); % gating variable m for Nav1.6, defined in equation (7)
h0_1_6 = y(:,26); % gating variable h0 for Nav1.6, defined in equation (11).
m_1_3 = y(:,27); % gating variable m for Nav1.3, defined in equation (7)
h0_1_3 = y(:,34); % gating variable h0 for Nav1.3, defined in equation (11).
m_1_2 = y(:,35); % gating variable m for Nav1.2, defined in equation (7)
h0_1_2 = y(:,42); % gating variable h0 for Nav1.2, defined in equation (11).
m_1_1 = y(:,43); % gating variable m for Nav1.1, defined in equation (7)
h0_1_1 = y(:,50); % gating variable h0 for Nav1.1 defined in equation (11).
m_1_5 = y(:,51); % gating variable m for Nav1.5, defined in equation (7)
h0_1_5 = y(:,58); % gating variable h0 for Nav1.5 defined in equation (11).
m_1_9 = y(:,59); % gating variable m for Nav1.9, defined in equation (5)
h_1_9 = y(:,60); % gating variable h for Nav1.9, defined in equation (6)


% calculate conductances
K = gK_max * n.^4; % potassium
Na_1_9 = gNa_max_1_9 * m_1_9.^3 .* h_1_9; % Nav1.9
Na_1_8 = gNa_max_1_8 * m_1_8.^3 .* (1-h0_1_8); % Nav1.8
Na_1_7 = gNa_max_1_7 * m_1_7.^3 .* (1-h0_1_7); % Nav1.7
Na_1_6 = gNa_max_1_6 * m_1_6.^3 .* (1-h0_1_6); % Nav1.6
Na_1_5 = gNa_max_1_5 * m_1_5.^3 .* (1-h0_1_5); % Nav1.5
Na_1_3 = gNa_max_1_3 * m_1_3.^3 .* (1-h0_1_3); % Nav1.3
Na_1_2 = gNa_max_1_2 * m_1_2.^3 .* (1-h0_1_2); % Nav1.2
Na_1_1 = gNa_max_1_1 * m_1_1.^3 .* (1-h0_1_1); % Nav1.1
L = gL; % leakage

% calculate currents
i_K = (K .* (V - VK)); % potassium
i_1_9 = (Na_1_9 .* (V - VNa)); % Nav1.9
i_1_8 = (Na_1_8 .* (V - VNa)); % Nav1.8
i_1_7 = (Na_1_7 .* (V - VNa)); % Nav1.7
i_1_6 = (Na_1_6 .* (V - VNa)); % Nav1.6
i_1_5 = (Na_1_5 .* (V - VNa)); % Nav1.5
i_1_3 = (Na_1_3 .* (V - VNa)); % Nav1.3
i_1_2 = (Na_1_2 .* (V - VNa)); % Nav1.2
i_1_1 = (Na_1_1 .* (V - VNa)); % Nav1.1
i_L = (L * (V - Vl)); % leakage

% format return values
Y=[i_K, i_1_9, i_1_8, i_1_7, i_1_6, i_1_5, i_1_3, i_1_2, i_1_1, i_L];
end

