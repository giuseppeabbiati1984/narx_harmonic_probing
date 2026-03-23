function [yr,fr] = fun_ODEeval(odeHandle,tr, fr, SNRf, SNRy)
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
[~,yr] = ode45(odeHandle, tr, yi) ; yr = yr(:,1) ; % displacement

% noise contamination
yr = yr + std(yr) / SNRy * randn(size(yr)) ;
fr = fr + std(fr) / SNRf * randn(size(fr)) ;

end