function H3_sol = fun_probing_H3(J0, J, JJ, JJJ, nf, nz, H1, H2, w1_ax, dts)

%This function performs Harmonic Probing numerically by minimizing
% different points on the frequency spectrum.
% profile on

%initial guesses
H0 = 0;

%Point at which to evaluate the residual z = H2
z1 = 1 ; 
z2 = 0 ; 

%% CTF (ONLY DIAGONAL FOR NOW)
H3_sol = ones( [ numel(w1_ax), numel(w1_ax),  numel(w1_ax)] )./1e5 ...
    +  1i*ones( [ numel(w1_ax), numel(w1_ax),  numel(w1_ax)]  )./1e5 ;


% CTF Optimization (w3 = fixed)
dw = w1_ax(2) - w1_ax(1) ; % all axes discretized with the same step

for k = 1 : numel(w1_ax)
    for j = 1 : numel(w1_ax)
        for i = 1 : numel(w1_ax)
	  w1 = w1_ax(i) ;     % in rad/s 1st harmonic
	  w2 = w1_ax(j) ;     % in rad/s 2nd harmonic
	  w3 = w1_ax(k) ;     % in rad/s 3rd harmonic

	   % Interp LTF at the correct w
	   H1w1 = interp1(w1_ax, H1, w1, "linear", "extrap")  ;
	   H1w2 = interp1(w1_ax, H1, w2, "linear", "extrap")  ;
	   H1w3 = interp1(w1_ax, H1, w3, "linear", "extrap")  ;


	   H2w1w1 = interp2(w1_ax, w1_ax, H2, w1, w1, 'spline', 0) ;
	   H2w2w2 = interp2(w1_ax, w1_ax, H2, w2, w2, 'spline', 0) ;
	   H2w3w3 = interp2(w1_ax, w1_ax, H2, w3, w3, 'spline', 0) ;


	   H2w1w2 = interp2(w1_ax, w1_ax, H2, w1, w2, 'spline', 0) ;
	   H2w1w3 = interp2(w1_ax, w1_ax, H2, w1, w3, 'spline', 0) ;
	   H2w2w3 = interp2(w1_ax, w1_ax, H2, w2, w3, 'spline', 0) ;


	  % Check for frequency mixing:
        if round(w1, 5)==0 || round(w2, 5)==0 || round(w3, 5)==0 || ...
	      w1+w2+w3 == 3*w1 || ...
	      w1+w2+w3 == 3*w2 || ...
	      w1+w2+w3 == 3*w3 || ...
	      w1+w2+w3 == 2*w1 + w3 || ...
	      w1+w2+w3 == 2*w1 + w2 || ...
	      w1+w2+w3 == 2*w2 + w1 || ...
	      w1+w2+w3 == 2*w1 + w3 || ...
	      w1+w2+w3 == 2*w3 + w1 || ...
	      w1+w2+w3 == 2*w3 + w2 

		H3_sol(i, j, k) = NaN ;
	  continue
        end


	  [~, t] = fun_freq_finder([abs(w1), abs(w2), abs(w3)], dw) ;						% find t that gives no leakage

	  [~, R, freq] = fun_residual3(J0, J, JJ, JJJ, H0, H1w1, H1w2, H1w3,...
	      H2w1w1, H2w2w2, H2w3w3, H2w1w2, H2w1w3, H2w2w3,...
	      z1, nf, nz, w1, w2, w3, dts, t) ;
	  [~, index] = min(abs(freq - (w1+w2+w3)/(2*pi) )) ;								% identify the correct index for the optimization
	  delta_1 =  R(index);


	  [~, R, freq] = fun_residual3(J0, J, JJ, JJJ, H0, H1w1, H1w2, H1w3,...
	      H2w1w1, H2w2w2, H2w3w3, H2w1w2, H2w1w3, H2w2w3,...
	      z2, nf, nz, w1, w2, w3, dts, t) ;
	  [~, index] = min(abs(freq - (w1+w2+w3)/(2*pi) )) ; % identify the correct index for the optimization
	  delta_0 =  R(index);

	  H3_sol(i, j, k) =   delta_0/(delta_0 - delta_1) ; % estimate the slope of the linear residual equation

        end
    end


    if rem(k, 5) ==0
        disp(strcat('3rd order HP, frequency: ', num2str(k), " of ", num2str(numel(w1_ax)))) ;
        % disp(' ')
    end

end

% interpolate the missing values (NaN) in plane w1-w2
for k = 1 : numel(w1_ax)
    for i=1:4  % Try several times to ensure all values are filled
        H3_sol(:,:,k) = fillmissing2(H3_sol(:,:,k), "cubic") ;
    end
end

% interpolate the missing values (NaN) in plane w1-w3
for k = 1 : numel(w1_ax)
	H3_sol(:,k,:) = fillmissing2(squeeze(H3_sol(:,k,:)), "cubic") ;
end


end
