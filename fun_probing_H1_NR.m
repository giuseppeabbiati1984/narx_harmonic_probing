function H1_sol = fun_probing_H1_NR(J0, J, nf, nz, w, dts)
% This function performs Harmonic Probing numerically by minimizing a time-domain residual
% It uses a Newton-Raphson solving scheme for minimizing the residual in 
% a, b - scaling factors on the output, input

dt = dts ;
t = 1000:dt:1000+dt*50;		% This parameter should cancel out. Make sure it's far from 0 to avoid infinities in the residual
H0 = 0;

%% LTF Optimization
H1_sol = ones( [1, numel(w)] ).*1e-5 +  1i*ones( [1, numel(w)] ).*1e-5 ;
tol = 1e-9;    % solver tolerance


% Initial guesses
H1w1 = 1e4 + 1e4i ; dH1 = 1e-6 ;

for i=1:1:numel(w)
    w1 = w(i) ;

    if i~=1
        H1w1 = H1_sol(i-1) ; % Initial guess
    end

    % Nonlinear Solver
    for k = 1:1:5 % Newton-Raphson Iterations
        [r, ~, ~] = fun_residual1(J0, J, H0, H1w1, nf, nz, w1, dt, t) ;

        if sum(abs(r))/numel(r) < tol
	  break
        end

        [r_back, ~, ~] = fun_residual1(J0, J, H0, H1w1-dH1,  nf, nz, w1, dt, t) ;
        [r_forw, ~, ~] = fun_residual1(J0, J, H0, H1w1+dH1, nf, nz, w1, dt, t) ;

        H1w1_new = H1w1 - r/((r_forw - r_back)/(2*dH1)) ;

        H1w1 = H1w1_new;

    end
    
    H1_sol(i) = H1w1;
end
 
end