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
ry = 3 ; % max lagged sample on output (k-1, ..., k-ry)
ru = 3 ; % max lagged sample on input  (k-0, ..., k-ru )
ord = "3" ; % order of the poly-NARX
numSeg = 2; % number of segments ( >1 )
lenSeg = 100 ; % number of points in each segment
fsr = 5000 ; % sampling frequency for ODE
dtr = 1/fsr ; % reference time step for ODE [s]
dts = 0.005236 ; % subsampling time step(s) for NARX [s]   opt = 5.236e-3
tmax = lenSeg*numSeg*dts;
tr = 0:dtr:tmax ; % high frequency time axis [s]
ts = 0:2*max(dts):tmax ; % sampling time axis [s]
fs = normrnd(0,5,[numel(ts), 1]) ; % input signal [N]
ur = interp1(ts, fs, tr) ; % Nyquist-proof input
clear ts fs

% ODE response and force with added noise
SNRf = 100 ; % signal-to-noise ratio for input
SNRy = 100 ; % signal-to-noise ratio for output
odeHandle = @(t, y) compute_SSduffing(t, y, odePars, tr, ur) ; % state-space definition of the duffing OED
														
[yr, ur] = compute_ODE(odeHandle, tr ,ur, SNRf, SNRy) ; % solving the OED

%Compute scaling factors (input and output)
stdU = std(ur) ; % input std. dev. for scaling
stdY = std(yr) ; % output std. dev. for scaling

%Compute transfer function scaling factors (nonlinear scaling)
scale0 = stdY/stdU^0 ; 
scale1 = stdY/stdU^1 ; % LTF Scale
scale2 = stdY/stdU^2 ; % QTF Scale
scale3 = stdY/stdU^3 ; % CTF Scale

% Scale the signals by the std. dev.
ur_scl = ur/stdU ; % Scaled input 
yr_scl = yr/stdY ; % Scaled output

save duff_train_data.mat 

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

