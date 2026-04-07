function H1_sol = fun_probing_H1(J0, J, nf, nz, w, dts)
% This function performs Harmonic Probing numerically by minimizing
% different points on the frequency spectrum.
% a, b - scaling factors on the output, input

%Point at which to evaluate the residual z = H2
z1 = 1; % Point 1 for eval of residual
z2 = 0; % Point 2 for eval of residual

%% LTF Optimization
H1_sol = ones( [1, numel(w)] ).*1e-5 +  1i*ones( [1, numel(w)] ).*1e-5 ;


for i=1:1:numel(w)
    w1 = w(i) ; % rad/s

    [~, t] = fun_freq_finder(5*w1, 1); % Find the length of the signal for the probing. Use 5*w for stability

    [~, R, freq] = fun_residual1(J0, J, 0, z1,  nf, nz, w1, dts, t) ; % Compute the residual for Point 1
    [~, index] = min(abs(freq - (w1)/(2*pi) )) ; % identify the correct index for the optimization
    delta_1 =  R(index);


    [~, R, freq] = fun_residual1(J0, J, 0, z2,  nf, nz, w1, dts, t) ; % Compute the residual for Point 2
    [~, index] = min(abs(freq - (w1)/(2*pi) )) ; % identify the correct index for the optimization
    delta_0 =  R(index);

    H1_sol(i) =  delta_0/(delta_0 - delta_1); 
end


end