function [r, R, freq] = fun_residual1(J0, J, H0, H1w1, nf, nz, w1, dts, t)
%This function calcualtes the residual between a linear NARX model and 1st order Volterra expansion
     
	   dt = (t(end) - t(1))/(length(t) - 1); % Obtain the time-step for computing the frequency axis of the DFT
	   fs = 1/dt;	% Obtain frequency discretization of the axis of the DFT


	   % vectorized implementation
	   z = exp(1i*w1*(t.' - dts*(0:1:(nz) )))  ; % Monochromatic input at frequency w1

	   f = H0  + H1w1 * exp(1i*w1*(t.' - dts*(0:1:(nf) )))  ;	% Volterra series expansion of order 1


	   X_vec = [f(:, 2:end), z] ;

	   r = f(:, 1).' - J0 - (J * X_vec.')  ;	% Time-domain residual. r = NARX - Volterra


%  Double-sided spectrum:
Y = fftshift( fft(r) ) ;
L = numel(t) ;

P2 = (Y/L) ; % Scale the amplitude

freq =  fs/L*(-(L)/2 : (L-1)/2) ; % Frequncy axis of the DFT in hertz = 1/s 
R = transpose(P2) ;
end
