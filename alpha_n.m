%% parametrization of voltage-dependency of gating variable alpha_n for generic potassium channel (equation (13))

function a_n = alpha_n(Vm, s)
% Vm: membrane potential in millivolts
% s: voltage shift compared to parametrization from original work from Hodgkin-Huxley
    Vm = Vm +s;
    a_n = (0.01 * (10.0 - Vm)) ./ (exp(1.0 - (0.1 * Vm)) - 1.0);
end

