function H1_sol = compute_probing_H1(C1, ry, ru, w1_ax, dts)
% This function estimates the linear transfer function based on a polynomial NARX model.
% This function performs Harmonic Probing numerically by minimizing different points on the frequency spectrum.
% The equation/Line references in this function refer to the paper:
% Stamenov et al. (2025) "Numerical estimation of generalized frequency response functions from time series data using NARX"
% https://doi.org/10.1016/j.ymssp.2025.113278

% INPUT:
% ry: Lags on output signal
% ru: Lags on input signal
% r = ry + ru + 1: Total length of the input vector
% C1: Vector of coefficients contributing to the linear NARX term. [1 x r]
% w1_ax: frequency axis for probing [rad/s]
% dts: time-step [s]

% OUTPUT:
% H1_sol: estimated complex values for the linear GFRF at frequencies defined in w1_ax

% Point at which to evaluate the residual. That is, z represents trial values for the GFRF (H1)
z1 = 1; % Point 1 for eval of residual
z2 = 0; % Point 2 for eval of residual

H1_sol = ones( [1, numel(w1_ax)] ).*1e-5 +  1i*ones( [1, numel(w1_ax)] ).*1e-5 ; % Declare an empty array
dw = abs(w1_ax(2) - w1_ax(1)) ; % Frequency grid discretization

for i=1:1:numel(w1_ax)
    w1 = w1_ax(i) ; % rad/s

    % Esnure the time vector does not cause leakage in the DFT computation (factor 5 introduced)
    dt = (pi)/(1*(sum(2*w1))); % Upper bound computation m*sum(\Omega_p). Line 7. Algorithm 1. (pg. 6)
    r =  single((2*pi/dw)/dt); % Length of the time-vector (must be an integer). Line 8+9. Algorithm 1. (pg. 6)
    t = 0:dt:(r-1)*dt ; % Time vector, last element removed for periodicity

    [~, R1, freq] = compute_residual1(C1, z1,  ry, ru, w1, dts, t) ; % Compute the residual for Point 1. Line 20 (pg. 6)
    [~, index1] = min(abs(freq - (w1)/(2*pi) )) ; % Identify the correct index for the optimization
    delta_1 =  R1(index1);

    [~, R0, freq] = compute_residual1(C1, z2,  ry, ru, w1, dts, t) ; % Compute the residual for Point 2. Line 19 (pg. 6)
    [~, index0] = min(abs(freq - (w1)/(2*pi) )) ; % Identify the correct index for the optimization
    delta_0 =  R0(index0);

    H1_sol(i) =  delta_0/(delta_0 - delta_1); % Equation on Line 21 from Algorithm 1 (pg. 6)
end


end