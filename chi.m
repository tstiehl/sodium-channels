
%% parameterization of voltage dependent gating variables
%% function chi from equation (1) of in silico supplement

function res = chi(par, x)
% par: channel-specific parameters
% x: membrane potential in volts
	res = (par(1)+par(2)*x)./(par(3)+par(4)*exp( (x+par(5))/par(6) ))+par(7);
end
