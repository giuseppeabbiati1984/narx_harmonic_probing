function [eps, Eps, f_ax] = compute_residual2(C1, C2, H1w1, H1w2, H2w1w2, ry, ru, w1, w2, dts, t)
%This function calcualtes the residual between a 2nd-order NARX model and 2nd order Volterra expansion
% The computation is vectorized for efficiency

% INPUT description:
% ry: Lags on output signal
% ru: Lags on input signal
% r = ry + ru + 1: Total length of the input vector
% C1: Vector of coefficients contributing to the linear NARX term. [1 x r]
% C2: Array of coefficients contributing to the quadratic NARX term. [r x r]
% H1w1: value of the linear GFRF at w1
% H1w2: value of the linear GFRF at w2
% H2w1w2: value of the quadratic GFRF at (w1,w2)
% w1: frequency value for probing [rad/s]
% w2: frequency value for probing [rad/s]
% dts: time-step [s]
% t: time-vector [s]

% OUTPUT:
% eps: time-domain residual, eps = Volterra - NARX. [1 x numel(t)]
% Eps: frequncy-domain residual [numel(t) x 1]
% f_ax: frequncy axis [Hz] [1 x numel(t)]



% fs = 1/(t(2) - t(1));
dt = (t(end) - t(1))/(length(t) - 1); % Obtain the time-step for computing the frequency axis of the DFT
fs = 1/dt; % Obtain frequency discretization of the axis of the DFT
w12 = w1 + w2; % sum frequency


u = exp(1i*w1*(t.' - dts*(0:1:(ru) ))) + ...
    exp(1i*w2*(t.' - dts*(0:1:(ru) ))) ; % Bi-chromatic input with freqs w1, w2. [numel(t) x (ru + 1)]

y = ...
    H1w1 * exp(1i*w1*(t.' - dts*(0:1:(ry) ))) + ... % Linear w1
    H1w2 * exp(1i*w2*(t.' - dts*(0:1:(ry) ))) + ... % Linear w2
    (H2w1w2 + (H2w1w2))  * exp(1i*w12*(t.' - dts*(0:1:(ry) )))  ; % Quad w1+w2. [numel(t) x (ru + 1)]


X_vec = [y(:, 2:end), u] ; % input vector. [numel(t) x r]


eps = y(:, 1).' - ((C1 * X_vec.')  +  1/2 * sum(X_vec * C2 .* X_vec, 2).');  % Time-domain residual. r = Volterra - NARX. [1 x numel(t)]. Eq. 24. (pg. 5)


%  Double-sided spectrum
Y = fftshift( fft(eps) ) ;
L = numel(t) ;

P2 = (Y/L) ; % Scale the amplitude

f_ax =  fs/L*(-(L)/2 : (L-1)/2) ; % Frequncy axis of the DFT in hertz = 1/s
Eps = transpose(P2) ;

end