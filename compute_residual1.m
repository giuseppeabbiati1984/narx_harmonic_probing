function [eps, Eps, f_ax] = compute_residual1(C1, H1w1, ry, ru, w1, dts, t)
%This function calcualtes the residual between a 1st-order NARX model and 1st-order Volterra expansion
% The computation is vectorized for efficiency

% INPUT description:
% ry: Lags on output signal
% ru: Lags on input signal
% r = ry + ru + 1: Total length of the input vector
% C1: Vector of coefficients contributing to the linear NARX term. [1 x r]
% H1w1: value of the linear GFRF at w1.
% w1: frequency value for probing [rad/s]
% dts: time-step [s]
% t: time-vector [s]

% OUTPUT:
% eps: time-domain residual, eps = Volterra - NARX. [1 x numel(t)]
% Eps: frequncy-domain residual [numel(t) x 1]
% f_ax: frequncy axis [Hz] [1 x numel(t)]

dt = (t(end) - t(1))/(length(t) - 1); % Obtain the time-step for computing the frequency axis of the DFT
fs = 1/dt;	% Obtain frequency discretization of the axis of the DFT


u = exp(1i*w1*(t.' - dts*(0:1:(ru) )))  ; % Monochromatic input at frequency w1. [numel(t) x (ru + 1)]
y = H1w1 * exp(1i*w1*(t.' - dts*(0:1:(ry) )))  ; % Volterra series expansion of order 1. [numel(t) x (ru + 1)]

X_vec = [y(:, 2:end), u]; % input vector. [numel(t) x r]

eps = y(:, 1).' - (C1 * X_vec.')  ; % Time-domain residual. Eq. 24. (pg. 5)


%  Double-sided spectrum:
Y = fftshift( fft(eps) ) ;
L = numel(t) ;

P2 = (Y/L) ; % Scale the amplitude

f_ax =  fs/L*(-(L)/2 : (L-1)/2) ; % Frequncy axis of the DFT in hertz = 1/s. [1 x numel(t)]
Eps = transpose(P2) ; % [numel(t) x 1]
end
