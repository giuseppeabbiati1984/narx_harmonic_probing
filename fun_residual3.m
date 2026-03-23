function [r, R, freq] = fun_residual3(J0, J, JJ, JJJ, H0, H1w1, H1w2, H1w3,...
	      H2w1w1, H2w2w2, H2w3w3, H2w1w2, H2w1w3, H2w2w3, H3, nf, nz, w1, w2, w3, dts, t) 

% This function calcualtes the residual of a third order system and a Volterra expansion

	   dt = (t(end) - t(1))/(length(t) - 1); 
	   fs = 1/dt; 

             w12 = w1 + w2 ;
             w13 = w1 + w3 ;
             w23 = w2 + w3 ;
             w123 = w1 + w2 + w3 ;

	 
	   z = exp(1i*w1*(t.' - dts*(0:1:(nz) ))) + ...
		exp(1i*w2*(t.' - dts*(0:1:(nz) ))) + ...
		exp(1i*w3*(t.' - dts*(0:1:(nz) ))) ;

	   %
	   f = H0  + ...   % bias term
	       H1w1 * exp(1i*w1*(t.' - dts*(0:1:(nf) ))) + ...								% linear w1
	       H1w2 * exp(1i*w2*(t.' - dts*(0:1:(nf) ))) + ...								% linear w2
	       H1w3 * exp(1i*w3*(t.' - dts*(0:1:(nf) ))) + ...								% linear w3
	       H2w1w1 * exp(1i*2*w1*(t.' - dts*(0:1:(nf) )))  + ...							% quad w12
	       H2w2w2 * exp(1i*2*w2*(t.' - dts*(0:1:(nf) )))  + ...							% quad w13
	       H2w3w3 * exp(1i*2*w3*(t.' - dts*(0:1:(nf) )))  + ...							% quad w12
	       (H2w1w2+H2w1w2) * exp(1i*w12*(t.' - dts*(0:1:(nf) )))  + ...						% quad w12
	       (H2w1w3+H2w1w3) * exp(1i*w13*(t.' - dts*(0:1:(nf) )))  + ...						% quad w13
	       (H2w2w3+H2w2w3) * exp(1i*w23*(t.' - dts*(0:1:(nf) )))  + ...						% quad w12
	       6*H3 * exp(1i*w123*(t.' - dts*(0:1:(nf) )))  ;								% cubic w123


n  = numel(J);														% system size
Nt = numel(t);														% time series length 


X  = [f(:, 2:end), z];														% Nt x n
f0 =  f(:, 1) ;															% Nt x 1


%========================
% Linear term: J*X'
%========================
lin = X*J.';															% Nt x 1

%========================
% Quadratic term: 1/2 * X*JJ*X'
%========================
quad = 0.5 * sum( (X*JJ) .* X , 2 ); 											 % Nt x 1

%========================
% Cubic term: 1/6 * sum_{a,b,c} JJJ(a,b,c) x_a x_b x_c
%========================

JJJ2 = reshape(JJJ, n*n, n); 												% (n^2) x n
vecS = JJJ2 * X.';														% (n^2) x Nt
S    = reshape(vecS, n, n, Nt);												% n x n x Nt  (pages are 3rd dim)

% Make X a row/col per page 
Xp = permute(X, [3 2 1]); 													% 1 x n x Nt  (rearanges the dimensions of the X matrix)
Xc = permute(X.', [1 3 2]);													% n x 1 x Nt  (rearanges the dimensions of the X matrix)

% cube(i) = (1/6) * X(i,:) * S(:,:,i) * X(i,:).'
cube = (1/6) * squeeze( pagemtimes( pagemtimes(Xp, S), Xc ) );							% Nt x 1


%========================
% Residual vector r (Nt x 1)
%========================
r = f0 - (J0 + lin + quad + cube);


% %DOUBLE SIDED SPECTRUM
Y = fftshift( fft(r) ) ;
L = numel(t) ;

P2 = Y/(L) ; % distribute the amplitude over all the time-steps

freq =  fs*(-L/2 : L/2-1)/L ; % in hertz = 1/s 
R = transpose(P2) ;

end
