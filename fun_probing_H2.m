function H2_sol = fun_probing_H2(J0, J, JJ, nf, nz, H1, w1_ax, dts)
% This function estimates the quadratic transfer function based on a polynomial NARX model.
% This function performs Harmonic Probing numerically by minimizing different points on the frequency spectrum.
% The equation references in this function refer to the paper:
% Stamenov et al. (2025) "Numerical estimation of generalized frequency response functions from time series data using NARX"
% https://doi.org/10.1016/j.ymssp.2025.113278

% H1= H1./ (a/b); %scale to correspond to normalized data

% Point at which to evaluate the residual. That is, z represents trial values for the GFRF (H2)
z1 = 1; % Point 1 for eval of residual
z2 = 0; % Point 2 for eval of residual

%% QTF 
H2_sol = ones( [numel(w1_ax), numel(w1_ax)] )./1e5 +  1i*ones( [numel(w1_ax), numel(w1_ax)] )./1e5 ; % Declare an empty array
% 


dw = abs(w1_ax(2) - w1_ax(1)) ; % Frequency grid discretization

for j = 1: numel(w1_ax)
    for i = 1: numel(w1_ax)
        w1 = w1_ax(i) ; % 1st harmonic in rad/s 
        w2 = w1_ax(j) ; % 2nd harmonic in rad/s 


        % Check for frequency mixing. Line 3 of Algorithm1. (pg. 6):
        if round(w1, 5)==0 || round(w2, 5)==0 ...
	      || w1+w2 == 2*w1 ...
	      || w1+w2 == 2*w2 ...
	      || w1 == w2

	  H2_sol(i, j) = NaN ; % If frequency cannot be probed due to mixing, set NaN for later interp
	  continue
        end

        % Interp LTF at the correct w
        H1w1 = interp1(w1_ax, H1, w1, "spline", "extrap")  ; % Ensure H1 value is available at w1
        H1w2 = interp1(w1_ax, H1, w2, "spline", "extrap")  ; % Ensure H1 value is available at w2

        [~, t] = fun_freq_finder([abs(w1), abs(w2)], dw) ; % find t that gives no leakage in the DFT


        [~, R, freq] = fun_residual2(J0, J, JJ, 0, H1w1, H1w2, z1, nf, nz, w1, w2, dts, t) ; % Compute the residual for Point 1
        [~, index] = min(abs(freq - (w1+w2)/(2*pi) )) ; % identify the correct index for the optimization
        delta_1 =  R(index);


        [~, R, freq] = fun_residual2(J0, J, JJ, 0, H1w1, H1w2, z2, nf, nz, w1, w2, dts, t) ; % Compute the residual for Point 2
        [~, index] = min(abs(freq - (w1+w2)/(2*pi) )) ; % identify the correct index for the optimization
        delta_0 =  R(index);

        H2_sol(i, j) =  delta_0/(delta_0 - delta_1); % Equation on Line 21 from Algorithm 1 (pg. 6)

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