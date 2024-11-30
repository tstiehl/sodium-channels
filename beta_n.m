%% parametrization of voltage-dependency of gating variable beta_n for generic potassium channel (equation (14))

function b_n = beta_n(Vm,s)
% Vm: membrane potential in millivolts
% s: voltage shift compared to parametrization from original work from Hodgkin-Huxley
    Vm = Vm +s;
    b_n = 0.125 * exp(-Vm ./ 80.0);
end

