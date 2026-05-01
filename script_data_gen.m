%% ######## GENERATE NUMERICAL DUFFING DATA  ########################
% This script generates synthetic data from a Duffing oscillator which is used for training of the NARX model 
% which is consequently probed with the numerical probing algorithm.

close all; clear; clc;
rng(10, 'twister') ;

%------------------------------------- Oscillator input: ----------------------------------------------------------------------------
odePars.m = 1 ; % mass [kg]
odePars.k1 = 1e4 ; % stiffness (linear) [N/m] 
odePars.k2 = 1e7 ; % stiffness (quadratic) [N/m2]
odePars.k3 = 5e9 ; % stiffness (cubic) [N/m3]
odePars.w0 = sqrt(odePars.k1/odePars.m) ; % modal frequency [rad/s]
odePars.f0 = odePars.w0/(2*pi) ; % modal frequency [Hz]
odePars.xi = 0.10 ; % damping ratio [  ]
odePars.c = 2*odePars.xi*sqrt(odePars.k1*odePars.m) ; % damping [Ns/m]

% Time Series Input:
fsr = 5000 ; % sampling frequency for ODE
dtr = 1/fsr ; % reference time step for ODE [s]
dts = 0.005 ; % subsampled time axis [s]  
tmax = 10; % time length [s];
tr = 0:dtr:tmax ; % high frequency time axis [s]
ts = 0:2*max(dts):tmax ; % sampling time axis [s]
us = normrnd(0, 5, [numel(ts), 1]) ; % input signal [N]
ur = interp1(ts, us, tr) ; % Nyquist-proof input

% ODE response and force with added noise
SNRf = 100 ; % signal-to-noise ratio for input
SNRy = 100 ; % signal-to-noise ratio for output
odeHandle = @(t, y) compute_SSduffing(t, y, odePars, tr, ur) ; % state-space definition of the duffing OED
														
[yr, ur] = compute_ODE(odeHandle, tr ,ur, SNRf, SNRy) ; % solving the OED

tr = tr(:);
ur = ur(:);
yr = yr(:);

save duff_train_data.mat tr ur yr

return
%% Synthetich data plot
PlotFontSize = 22;

figure('Renderer', 'painters', 'Position', [10, 100, 1000, 500]) ;

subplot('Position', [0.15, 0.62, 0.75, 0.32]);
plot(tr, ur, "r-") ;
% xlabel('$t$','Interpreter','latex')
ylabel('$u(t)$','Interpreter','latex')
set(gca, 'FontSize',  PlotFontSize) ;


subplot('Position', [0.15, 0.15, 0.75, 0.32]);
plot(tr, yr, "b-") ; 
xlabel('$t$','Interpreter','latex')
ylabel('$y(t)$','Interpreter','latex')
set(gca, 'FontSize',  PlotFontSize) ;


folder = "C:\Users\AU657332\OneDrive - Aarhus universitet\Giuseppe Abbiatis filer - david_stamenov\dissemination\MethodsX" ;  % your target folder
filename = strcat("synth_data.eps");  % filename
fullpath = fullfile(folder, filename);  % create full path

% print(gcf, '-depsc', fullpath);	% save as color EPS

