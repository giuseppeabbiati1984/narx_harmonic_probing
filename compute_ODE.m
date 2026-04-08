function [yr, ur] = fun_ODEeval(odeHandle,tr, ur, SNRf, SNRy)
% This function evaluates an ODE in state-space 

% Input:
% odeHandle: handle to the ODE in state-space form
% tr: time axis
% fr: input loading
% SNRf: signal-to-noise ratio on input
% SNRy: signal-to-noise ratio on output

% Output:
% yr: output of the system
% fr: input signal of the system

% initial conditions
ui = 0 ;
vi = 0 ;
yi = [ui;vi] ;

% time integration
[~,yr] = ode45(odeHandle, tr, yi) ; yr = yr(:,1) ; % Displacement

% noise contamination
yr = yr + std(yr) / SNRy * randn(size(yr)) ; % Add noise to the output. Eq. 28. (pg. 6)
ur = ur + std(ur) / SNRf * randn(size(ur)) ; % Add noise to the output. Eq. 28. (pg. 6)

end