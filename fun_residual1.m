function [r, R, freq] = fun_residual1(J0, J, H0, H1w1, nf, nz, w1, dts, t)
%this function calcualtes the residual due to a second order linearization
%and Volterra probing
     
             % fs = 1/(t(2) - t(1));
	   dt = (t(end) - t(1))/(length(t) - 1);										% Updated after talking with Øyvind, Feb 27, 2024
	   fs = 1/dt;														% Updated after talking with Øyvind, Feb 27, 2024


	   % vectorized implementation
	   z = exp(1i*w1*(t.' - dts*(0:1:(nz) )))  ;

	   f = H0  + H1w1 * exp(1i*w1*(t.' - dts*(0:1:(nf) )))  ;								% Constant + Linear term


	   X_vec = [f(:, 2:end), z] ;

	   r = f(:, 1).' - J0 - (J * X_vec.')  ;											% time-domain residual


%  Double-sided spectrum
Y = fftshift( fft(r) ) ;
L = numel(t) ;

P2 = (Y/L) ;															% distribute the amplitude over all the time-steps

freq =  fs/L*(-(L)/2 : (L-1)/2) ;												% in hertz = 1/s 
R = transpose(P2) ;
end
