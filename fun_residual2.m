function [r, R, freq] = fun_residual2(J0, J, JJ, H0, H1w1, H1w2, H2w1w2, nf, nz, w1, w2, dts, t)
%this function calcualtes the residual due to a second order linearization
%and Volterra probing
     
             % fs = 1/(t(2) - t(1));
	   dt = (t(end) - t(1))/(length(t) - 1); % Updated after talking with Øyvind, Feb 27, 2024
	   fs = 1/dt;				   % Updated after talking with Øyvind, Feb 27, 2024
             w12 = w1 + w2;


	   % vectorized implementation
	   z = exp(1i*w1*(t.' - dts*(0:1:(nz) ))) + ...
	         exp(1i*w2*(t.' - dts*(0:1:(nz) ))) ;

	   f = H0  + ...													% Bias term
	       H1w1 * exp(1i*w1*(t.' - dts*(0:1:(nf) ))) + ...								% Linear w1
	       H1w2 * exp(1i*w2*(t.' - dts*(0:1:(nf) ))) + ...								% Linear w2
	       (H2w1w2 + (H2w1w2))  * exp(1i*w12*(t.' - dts*(0:1:(nf) )))  ;						% Quad w1+w2


	   X_vec = [f(:, 2:end), z] ;

	   r = f(:, 1).' - (J0 + (J * X_vec.')  +  1/2 * sum(X_vec * JJ .* X_vec, 2).');


%  Double-sided spectrum
Y = fftshift( fft(r) ) ;
L = numel(t) ;

P2 = (Y/L) ;															% distribute the amplitude over all the time-steps

freq =  fs/L*(-(L)/2 : (L-1)/2) ;												% in hertz = 1/s 
% freq_w = freq*2*pi ;
R = transpose(P2) ;

end

% slow implementation (looping)
% for i=1:numel(t)
%     % bi-chromatic excitation
%     z = exp(1i*w1*(t(i) - dts*(0:1:(nz) ))) + ...
%           exp(1i*w2*(t(i) - dts*(0:1:(nz) ))) ;
%
%     % response up to QTF
%     f = H0  + ...														% bias term
%         H1w1 * exp(1i*w1*(t(i) - dts*(0:1:(nf) ))) + ...									% linear w1
%         H1w2 * exp(1i*w2*(t(i) - dts*(0:1:(nf) ))) + ...									% linear w2
%  H2w1w1 * exp(1i*2*w1*(t(i) - dts*(0:1:(nf) )))  + ...								% quad 2w1
%          H2w2w2 * exp(1i*2*w2*(t(i) - dts*(0:1:(nf) )))  + ...								% quad 2w2
% (H2w1w2 + conj(H2w1w2))  * exp(1i*w12*(t(i) - dts*(0:1:(nf) )))  ;						% quad w1+w2
%
%
%     X = [f(2:end), z] ;
%     f  = f(1) ;
%
%     %Residual in time domain
%     r(i) = f - (J0 + J * X.' +  1/2 * X * JJ * X.') ;
% end