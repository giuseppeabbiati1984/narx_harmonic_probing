function H3_sol = compute_probing_H3(C1, C2, C3, ry, ru, H1, H2, w1_ax, dts)
% This function estimates the cubic transfer function based on a polynomial NARX model.
% This function performs Harmonic Probing numerically by minimizing different points on the frequency spectrum.
% The equation references in this function refer to the paper:
% Stamenov et al. (2025) "Numerical estimation of generalized frequency response functions from time series data using NARX"
% https://doi.org/10.1016/j.ymssp.2025.113278

% INPUT:
% ry: Lags on output signal
% ru: Lags on input signal
% r = ry + ru + 1: Total length of the input vector
% C1: Vector of coefficients contributing to the linear NARX term. [1 x r]
% C2: Array of coefficients contributing to the quadratic NARX term. [r x r]
% C3: Tensor of coefficients contributing to the cubic NARX term. [r x r x r]
% w1_ax: frequency axis for probing [rad/s]
% dts: time-step [s]
% H1: linear GFRF
% H2: quadratic GFRF

% OUTPUT:
% H3_sol: estimated complex values for the cubic GFRF at frequency-triplets defined by the grid (w1_ax X w1_ax X w1_ax)

% Point at which to evaluate the residual. That is, z represents trial values for the GFRF (H3)
z1 = 1; % Point 1 for eval of residual
z2 = 0; % Point 2 for eval of residual

H3_sol = ones( [ numel(w1_ax), numel(w1_ax),  numel(w1_ax)] )./1e5 ...
    +  1i*ones( [ numel(w1_ax), numel(w1_ax),  numel(w1_ax)]  )./1e5 ; % Declare an empty array

dw = w1_ax(2) - w1_ax(1) ; % Frequency grid discretization
% Loop through all freqeuncy axes
for k = 1 : numel(w1_ax)
    for j = 1 : numel(w1_ax)
        for i = 1 : numel(w1_ax)
	  w1 = w1_ax(i) ; % 1st harmonic  in rad/s 
	  w2 = w1_ax(j) ; % 2nd harmonic  in rad/s 
	  w3 = w1_ax(k) ; % 3rd harmonic  in rad/s 

	   % Interp LTF at the correct w
	   H1w1 = interp1(w1_ax, H1, w1, "linear", "extrap")  ; % Ensure H1 value is available at w1
	   H1w2 = interp1(w1_ax, H1, w2, "linear", "extrap")  ; % Ensure H1 value is available at w2
	   H1w3 = interp1(w1_ax, H1, w3, "linear", "extrap")  ; % Ensure H1 value is available at w3

	   H2w1w1 = interp2(w1_ax, w1_ax, H2, w1, w1, 'spline', 0) ; % Ensure H2 value is available at w1,w1
	   H2w2w2 = interp2(w1_ax, w1_ax, H2, w2, w2, 'spline', 0) ; % Ensure H2 value is available at w2,w2
	   H2w3w3 = interp2(w1_ax, w1_ax, H2, w3, w3, 'spline', 0) ; % Ensure H2 value is available at w3,w3

	   H2w1w2 = interp2(w1_ax, w1_ax, H2, w1, w2, 'spline', 0) ; % Ensure H2 value is available at w1,w2
	   H2w1w3 = interp2(w1_ax, w1_ax, H2, w1, w3, 'spline', 0) ; % Ensure H2 value is available at w1,w3
	   H2w2w3 = interp2(w1_ax, w1_ax, H2, w2, w3, 'spline', 0) ; % Ensure H2 value is available at w2,w3


	   % Check for frequency mixing. Line 3 of Algorithm 1. (pg. 6):
	   w_check = compute_freq_check(3, [w1, w2, w3]); % compute the unidentifiable frequencies
	   if sum(w1+w2+w3 == w_check) >0
	       H3_sol(i, j, k) = NaN ; % If frequency cannot be probed due to mixing, set NaN for later interp
	       continue
	   end




	  [~, t] = find_freq([abs(w1), abs(w2), abs(w3)], dw) ; % find t that gives no leakage

	  [~, R, freq] = compute_residual3(C1, C2, C3, H1w1, H1w2, H1w3,...
	      H2w1w1, H2w2w2, H2w3w3, H2w1w2, H2w1w3, H2w2w3,...
	      z1, ry, ru, w1, w2, w3, dts, t) ; % Compute the residual for Point 1. Line 20 (pg. 6)
	  [~, index] = min(abs(freq - (w1+w2+w3)/(2*pi) )) ; % Identify the correct index (freqeuncy bin)
	  delta_1 =  R(index);


	  [~, R, freq] = compute_residual3(C1, C2, C3, H1w1, H1w2, H1w3,...
	      H2w1w1, H2w2w2, H2w3w3, H2w1w2, H2w1w3, H2w2w3,...
	      z2, ry, ru, w1, w2, w3, dts, t) ; % Compute the residual for Point 1. Line 19 (pg. 6)
	  [~, index] = min(abs(freq - (w1+w2+w3)/(2*pi) )) ; % Identify the correct index (freqeuncy bin)
	  delta_0 =  R(index);

	  H3_sol(i, j, k) =   delta_0/(delta_0 - delta_1) ; % Equation on Line 21 from Algorithm 1 (pg. 6)

        end
    end

    if rem(k, 5) ==0 % Print probing progress
        disp(strcat('3rd order HP, frequency: ', num2str(k), " of ", num2str(numel(w1_ax)))) ;
        disp(' ')
    end
end

% Interpolate the missing values (NaN) in plane w1-w2
for k = 1 : numel(w1_ax)
    for i=1:4  % Try several times to ensure all values are filled
        H3_sol(:,:,k) = fillmissing2(H3_sol(:,:,k), "cubic") ;
    end
end

% Interpolate the missing values (NaN) in plane w1-w3
for k = 1 : numel(w1_ax)
	H3_sol(:,k,:) = fillmissing2(squeeze(H3_sol(:,k,:)), "cubic") ;
end

end
