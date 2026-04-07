function H2_sol = fun_probing_H2(J0, J, JJ, nf, nz, H1, w1_ax, dts)
%This function performs Harmonic Probing numerically by minimizing
% different points on the frequency spectrum.

% H1= H1./ (a/b); %scale to correspond to normalized data

%Point at which to evaluate the residual z = H2
z1 = 1; % Point 1 for eval of residual
z2 = 0; % Point 2 for eval of residual

%% QTF 
H2_sol = ones( [numel(w1_ax), numel(w1_ax)] )./1e5 +  1i*ones( [numel(w1_ax), numel(w1_ax)] )./1e5 ;	% Array for initial guesses
% 

% QTF Optimization 
dw = abs(w1_ax(2) - w1_ax(1)) ;

% parfor i = 1 : numel(w1_ax)
for j = 1: numel(w1_ax)
    for i = 1: numel(w1_ax)
        w1 = w1_ax(i) ; % in rad/s 1st harmonic
        w2 = w1_ax(j) ; % in rad/s 2nd harmonic


         % Check for frequency mixing:
        if round(w1, 5)==0 || round(w2, 5)==0 ...
	      || w1+w2 == 2*w1 ...
	      || w1+w2 == 2*w2 ...
	      || w1 == w2 

	  H2_sol(i, j) = NaN ; % If freq cannot be obtained set NaN for later interp
	  continue
        end

        % Interp LTF at the correct w
        H1w1 = interp1(w1_ax, H1, w1, "spline", "extrap")  ;
        H1w2 = interp1(w1_ax, H1, w2, "spline", "extrap")  ;

        [~, t] = fun_freq_finder([abs(w1), abs(w2)], dw) ; % find t that gives no leakagein the FFT


        [~, R, freq] = fun_residual2(J0, J, JJ, 0, H1w1, H1w2, z1, nf, nz, w1, w2, dts, t) ; % Compute the residual for Point 1
        [~, index] = min(abs(freq - (w1+w2)/(2*pi) )) ; % identify the correct index for the optimization
        delta_1 =  R(index);


        [~, R, freq] = fun_residual2(J0, J, JJ, 0, H1w1, H1w2, z2, nf, nz, w1, w2, dts, t) ; % Compute the residual for Point 2
        [~, index] = min(abs(freq - (w1+w2)/(2*pi) )) ; % identify the correct index for the optimization
        delta_0 =  R(index);

        H2_sol(i, j) =  delta_0/(delta_0 - delta_1);

        if rem(j, 25) == 0
	  disp(strcat('2nd order HP, w1 =', num2str(w1), ', w2 =', num2str(w2), ', w1+w2 =', num2str(w1+w2)))
	  disp(' ')
        end
    end
end

% interpolate the missing values (NaN)
for i=1:4 % ensure all values are filled
    % mov_wind = floor(sqrt(numel(w1_ax)))+2 ;
    % H2_sol = fillmissing2(H2_sol, "movmean",  {mov_wind, mov_wind}) ;
    H2_sol = fillmissing2(H2_sol, "cubic") ;
end

H2_sol = fillmissing2(H2_sol, "nearest") ;

end