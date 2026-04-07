function [eps, R, freq] = fun_residual3(C1, C2, C3, H1w1, H1w2, H1w3,...
    H2w1w1, H2w2w2, H2w3w3, H2w1w2, H2w1w3, H2w2w3, H3w1w2w3, ry, ru, w1, w2, w3, dts, t)
%This function calcualtes the residual between a 3rd-order NARX model and 3rd-order Volterra expansion
% The computation is vectorized for efficiency

% INPUT description:
% ry: Lags on output signal
% ru: Lags on input signal
% r = ry + ru + 1: Total length of the input vector
% C1: Vector of coefficients contributing to the linear NARX term. [1 x r]
% C2: Array of coefficients contributing to the quadratic NARX term. [r x r]
% C3: Tensor of coefficients contributing to the cubic NARX term. [r x r x r]
% H1w1: value of the linear GFRF at w1
% H1w2: value of the linear GFRF at w2
% H1w3: value of the linear GFRF at w3
% H2w1w2: value of the quadratic GFRF at (w1,w2)
% .
% .
% H3w1w2w3: value of the cubic GFRF at (w1,w2, w3) 
% w1: frequency value for probing [rad/s]
% w2: frequency value for probing [rad/s]
% w3: frequency value for probing [rad/s]
% dts: time-step [s]
% t: time-vector [s]


dt = (t(end) - t(1))/(length(t) - 1); % Obtain the time-step for computing the frequency axis of the DFT
fs = 1/dt;  % Obtain frequency discretization of the axis of the DFT

w12 = w1 + w2 ;
w13 = w1 + w3 ;
w23 = w2 + w3 ;
w123 = w1 + w2 + w3 ;

% Tri-chromatic input at w1, w2, w3 [numel(t) x (ru + 1)]:
u = exp(1i*w1*(t.' - dts*(0:1:(ru) ))) + ...
    exp(1i*w2*(t.' - dts*(0:1:(ru) ))) + ...
    exp(1i*w3*(t.' - dts*(0:1:(ru) ))) ; 

% Volterra output [numel(t) x (ru + 1)] :
y = ...
    H1w1 * exp(1i*w1*(t.' - dts*(0:1:(ry) ))) + ...	% linear w1
    H1w2 * exp(1i*w2*(t.' - dts*(0:1:(ry) ))) + ...	 % linear w2
    H1w3 * exp(1i*w3*(t.' - dts*(0:1:(ry) ))) + ...	 % linear w3
    H2w1w1 * exp(1i*2*w1*(t.' - dts*(0:1:(ry) )))  + ...	 % quad w12
    H2w2w2 * exp(1i*2*w2*(t.' - dts*(0:1:(ry) )))  + ... % quad w13
    H2w3w3 * exp(1i*2*w3*(t.' - dts*(0:1:(ry) )))  + ... % quad w12
    (H2w1w2+H2w1w2) * exp(1i*w12*(t.' - dts*(0:1:(ry) )))  + ... % quad w12
    (H2w1w3+H2w1w3) * exp(1i*w13*(t.' - dts*(0:1:(ry) )))  + ... % quad w13
    (H2w2w3+H2w2w3) * exp(1i*w23*(t.' - dts*(0:1:(ry) )))  + ... % quad w12
    6*H3w1w2w3 * exp(1i*w123*(t.' - dts*(0:1:(ry) )))  ; % cubic w123. 


r  = ry + ru + 1; % system size

X_vec  = [y(:, 2:end), u];	% Input vector. [numel(t) x r]

%========================
% Linear NARX term: C1*X'
%========================
lin = X_vec*C1.'; % [numel(t) x 1]

%========================
% Quadratic NARX term: 1/2 * X*C2*X'
%========================
quad = 0.5 * sum( (X_vec*C2) .* X_vec , 2 ); % [numel(t) x 1]

%========================
% Cubic NARX term: 1/6 * 1/2 * X*C3*X'*X'
%========================
C3_stack = reshape(C3, r*r, r); % stack the 3D tensor into a 2D array % [(r^2) x r]
vecS = C3_stack * X_vec.'; % compute 1 page for each time-step [(r^2) x numel(t)]
S  = reshape(vecS, r, r, numel(t)); % [r x r x numel(t)]  (pages are 3rd dim)

% Make X a row/col per page
Xp = permute(X_vec, [3 2 1]); % [1 x r x numel(t)  (rearanges the dimensions of the X matrix)
Xc = permute(X_vec.', [1 3 2]); % [r x 1 x numel(t)]  (rearanges the dimensions of the X matrix)

% cube(i) = (1/6) * X(i,:) * S(:,:,i) * X(i,:).'
cube = (1/6) * squeeze( pagemtimes( pagemtimes(Xp, S), Xc ) ); % [numel(t) x 1]
% pagemtimes: multiplies two 2D arrays over the 3rd dimension

%========================
% Residual vector [numel(t) x 1]
%========================
eps = y(:, 1) - (lin + quad + cube); % Time-domain residual. r = Volterra - NARX. [numel(t) x 1]. Eq. 24. (pg. 5)


%  Double-sided spectrum
Y = fftshift( fft(eps) ) ;
L = numel(t) ;

P2 = Y/(L) ; % Scale the amplitude

freq =  fs*(-L/2 : L/2-1)/L ; % Frequncy axis of the DFT in hertz = 1/s
R = transpose(P2) ;

end
